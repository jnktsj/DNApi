#! /usr/bin/env python3

import sys
import os.path
import re
import uuid
import signal
import fileinput
import subprocess
from argparse import ArgumentParser

import dnapilib
from dnapilib.io_utils import get_file_obj
from dnapilib.apred import adapter_prediction
from dnapilib.apred import iterative_adapter_prediction
from dnapilib.exhaust import rm_temp_dir
from dnapilib.exhaust import fastq_input_prep
from dnapilib.exhaust import map_clean_reads
from dnapilib.exhaust import make_stats_report


TEMP_DIR = None
MAP_TO_GENOME = False
SAMPLE_NUM = 50000


def convert_interval(s_in, s_op, func):
    """Return range of kmers or filtering ratios.

    """
    msg = "bad {}: {} {}"
    try:
        s = list(map(func, s_in.split(":")))
    except:
        raise Exception(msg.format("value", s_op, s_in))
    if len(s) == 1:
        return s
    if len(s) == 3:
        beg, end, interval = s
        values = []
        while beg < end:
            values.append(beg)
            beg += interval
        values.append(end)
        return values
    else:
        raise Exception(msg.format("interval", s_op, s_in))


def parse_args():
    """Return options and required arguments.

    """
    parser = ArgumentParser(
                 usage="%(prog)s [options] FASTQ",
                 description="Predict or evaluate 3'adapter sequence(s)",
                 epilog="Report bug to: Junko Tsuji <jnktsj@gmail.com>")

    parser.add_argument("FASTQ",
        type=str,
        help="including stdin or compressed file {zip,gz,tar,bz}")
    parser.add_argument("--version", action="version",
        version="%(prog)s {}".format(dnapilib.__version__))

    predop = parser.add_argument_group("adapter prediction parameters")
    predop.add_argument("-k",
        metavar="[KMER_BEG:KMER_END:INCREMENT | KMER_LEN]",
        default="9:11:2",
        help="range of kmers or a single kmer to predict 3'adapters "
             "(default: %(default)s)")
    predop.add_argument("-r",
        metavar="[RATIO_BEG:RATIO_END:INTCREMENT | RATIO]",
        default="1.2:1.4:0.1",
        help="range of ratios or a single ratio to filter less abundant kmers"
             " (default: %(default)s)")
    predop.add_argument("--show-all",
        action="store_true",
        help="show other candidates if any")

    exhaop = parser.add_argument_group("exhaustive adapter search")
    exhaop.add_argument("--map-command",
        metavar="COMMAND",
        default=None,
        help="read mapping command to be tested")
    exhaop.add_argument("--subsample-rate",
        metavar="FLOAT",
        default=1.0, type=float,
        help="subsampling fraction of reads (default: %(default)s)")
    exhaop.add_argument("--output-dir",
        metavar="DIRECTORY",
        default="./dnapi_out",
        help="output directory to write report and cleansed reads"
             " (default: ./dnapi_out)")
    exhaop.add_argument("--no-output-files",
        action="store_true",
        help="only display report and suppress output files")
    exhaop.add_argument("--temp-path",
        metavar="DIRECTORY",
        default="/tmp",
        help="place to make temporary directory (default: %(default)s)")

    evalop = parser.add_argument_group("evaluation of candidate adapters")
    evalop.add_argument("--adapter-seq",
        dest="seq", nargs="+",
        default=None,
        help="list of 3'adapters for evaluation")

    adrmop = parser.add_argument_group("adapter removal parameters")
    adrmop.add_argument("--prefix-match",
        metavar="LENGTH",
        default=7, type=int,
        help="3'adapter match length to trim (default: %(default)s)")
    adrmop.add_argument("--min-len",
        metavar="LENGTH",
        default=16, type=int,
        help="minimum read length to keep for mapping (default: %(default)s)")
    adrmop.add_argument("--max-len",
        metavar="LENGTH",
        default=36, type=int,
        help="maximum read length to keep for mapping (default: %(default)s)")
    adrmop.add_argument("--trim-5p",
        metavar="LENGTH",
        default=0, type=int,
        help="trim specified number of bases from 5'ends after adapter removal"
             " (default: %(default)s)")
    adrmop.add_argument("--trim-3p",
        metavar="LENGTH",
        default=0, type=int,
        help="trim specified number of bases from 3'ends after adapter removal"
             " (default: %(default)s)")

    args = parser.parse_args()

    if args.map_command:
        err_find = "can't find {}"
        soft = os.path.expanduser(args.map_command.split()[0])
        if os.path.dirname(soft):
            if not os.path.exists(soft):
                raise Exception(err_find.format(soft))
        else:
            try:
                subprocess.call("which {}".format(soft).split())
            except OSError:
                raise Exception(err_find.format(soft))
        if not re.findall("@in", args.map_command):
            raise Exception("can't locate input argument: @in")
        if not re.findall("@out", args.map_command):
            raise Exception("can't locate output argument: @out")
        if args.prefix_match <= 0:
            raise Exception("bad value: --prefix-match")
        if args.min_len <= 0:
            raise Exception("bad value: --min-len")
        if args.max_len <= 0:
            raise Exception("bad value: --max-len")
        if args.trim_5p < 0:
            raise Exception("bad value: --trim-5p")
        if args.trim_3p < 0:
            raise Exception("bad value: --trim-3p")
        if args.subsample_rate <= 0 or 1 < args.subsample_rate:
            raise Exception("bad subsampling rate")
        global MAP_TO_GENOME
        MAP_TO_GENOME = True

    return args


def main():
    args = parse_args()
    fastq = args.FASTQ

    Ks = convert_interval(args.k, "-k", int)
    Rs = convert_interval(args.r, "-r", float)

    if not MAP_TO_GENOME:
        if len(Ks) > 1 or len(Rs) > 1:
            adapts = iterative_adapter_prediction(fastq, Rs, Ks, SAMPLE_NUM)
        else:
            adapts = adapter_prediction(fastq, Rs[0], Ks[0], SAMPLE_NUM)
        if args.show_all:
            for x in adapts:
                print("{}\tscore={:.2f}".format(*x))
        else:
            print(adapts[0][0])

    else:
        global TEMP_DIR
        TEMP_DIR = "{}/DNApi_tmp_{}".format(
            args.temp_path, str(uuid.uuid4()))
        subprocess.call(("mkdir {}".format(TEMP_DIR)).split())

        original_fastq = fastq
        fastq, total_read, sd = fastq_input_prep(
            fastq, args.subsample_rate, TEMP_DIR)

        if args.seq:
            adapts = set(args.seq)
            setstr = ["user-input" for i in range(len(adapts))]
        else:
            msg = "warning: predicted adapter is too short (<{0}): '{1}'\n" \
                + "warning: '{1}' will not be further investigated\n"
            params = {}
            for k in Ks:
                for r in Rs:
                    aout = adapter_prediction(fastq, r, k, SAMPLE_NUM)[0][0]
                    if len(aout) < args.prefix_match:
                        sys.stderr.write(msg.format(l, s))
                        continue
                    aseq = aout[: args.prefix_match+5]
                    params.setdefault(aseq,[]).append("{}:{:.1f}".format(k,r))
            adapts = list(params.keys())
            setstr = [';'.join(s) for s in params.values()]
            adapts.append("RAW_INPUT")
            setstr.append("NO_TREATMENT")

        if not adapts:
            raise Exception("no valid adapters to further process")

        table = []
        for i, aseq in enumerate(adapts):
            cnts = map_clean_reads(
                       fastq, aseq[:args.prefix_match], args.trim_5p,
                       args.trim_3p, args.min_len, args.max_len,
                       args.map_command, TEMP_DIR)
            read_stats = [c / total_read * 100 for c in cnts]
            table.append([aseq, cnts[0], read_stats[0],
                          cnts[1], read_stats[1], setstr[i]])
        make_stats_report(
            table, total_read, args.subsample_rate, args.prefix_match,
            sd, original_fastq, args.output_dir, TEMP_DIR, args.no_output_files)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    if sys.version_info.major <= 2:
        raise ValueError("DNApi requires python version 3 or higher")
    try:
        main()
    except KeyboardInterrupt:
        rm_temp_dir(TEMP_DIR)
    except Exception as e:
        prog = os.path.basename(sys.argv[0])
        rm_temp_dir(TEMP_DIR)
        sys.exit("{}: error: {}".format(prog, str(e)))
    finally:
        rm_temp_dir(TEMP_DIR)
