#!/usr/bin/env python3

"""Perform quality trimming with the same algorithm as
   bwa_trim_read() in bwaseqio.c, BWA. For Solexa quliaty,
   the scores are converted to Phred quality for trimming.

   Formula to convert Solexa quality to Phred quality is:
   Phred = 10 * log_10(1 + 10 ** (Solexa / 10.0))

   Formulas to calculate Phred and Solexa quality scores
   from sequencing error probability are:
   Phred  = -10 * log_10(P)
   Solexa = -10 * log_10(P / (1 - P))

"""

import sys
import re
import os.path
import signal
import math
from argparse import ArgumentParser


cur = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(cur))
from dnapilib.io_utils import get_file_obj, fastq_record


def solexa_to_phred(x):
    return int(round(10 * math.log10(1+10**(x/10.0))))


def illumina_64B(x):
    if x == 2:
        return 0
    else:
        return x


def illumina_33(x):
    return x


def calc_qual_score(p, solexa):
    if solexa:
        d = 1.0 - p if not p else 1
        return solexa_to_phred(int(-10 * math.log10(p/d)))
    else:
        return int(-10 * math.log10(p))


def qual_trim(args):
    if args.solexa:
        args.b = 64
        func = solexa_to_phred
    elif args.illumina5:
        func = illumina_64B
    else:
        func = illumina_33

    if args.b not in (33, 64):
        raise Exception("wrong quality score base")
    if args.l < 1:
        raise Exception("specify longer read length")
    if args.p < 0 or args.p > 1:
        raise Exception("bad error probability cutoff")
    if not args.solexa and args.q < 0:
        raise Exception("bad quality score cutoff")

    if args.q:
        cutoff = args.q
    else:
        cutoff = calc_qual_score(args.p, args.solexa)

    base = args.b
    minlen = args.l
    ns = re.compile('N', re.IGNORECASE)
    fastqs = fastq_record(get_file_obj(args.FASTQ))
    for read in fastqs:
        read = read.rstrip().split("\n")
        qual = read[3]
        s, max_s = 0, 0
        max_i = len(read[3])
        if minlen > max_i:
            continue
        for i in reversed(range(max_i)):
            q = func(ord(qual[i]) - base)
            s += cutoff - q
            if s < 0:
                break
            if s > max_s:
                max_s, max_i = s, i
        read[1] = read[1][:max_i]
        read[3] = read[3][:max_i]
        n_num = len(ns.findall(read[1]))
        if n_num < len(read[1]) and len(read[1]) >= minlen:
            print("\n".join(read))


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    prog = os.path.basename(sys.argv[0])
    if sys.version_info.major <= 2:
        raise ValueError("{} requires python version 3 or higher".format(prog))

    parser = ArgumentParser(
                 description="Perform quality trimming for single-end reads.")
    parser.add_argument("FASTQ",
        type=str,
        help="including stdin or compressed file {zip,gz,tar,bz}")
    parser.add_argument("-b",
        metavar="BASE",
        type=int, default=33,
        help="ASCII-encoded quality offset, e.g. 33 or 64 (default: %(default)s)")
    parser.add_argument("-p",
        metavar="PROB",
        type=float, default=0.1,
        help="error probability cutoff (default: %(default)s)")
    parser.add_argument("-q",
        metavar="SCORE",
        type=int, default=0,
        help="quality score cutoff (default: '-p 0.1')")
    parser.add_argument("-l",
        type=int, default=16,
        metavar="BP",
        help="minimum read length in bp (default: %(default)s)")
    parser.add_argument("--illumina5",
        action="store_true",
        help="Illumina 1.5+ encoding marked with 'B'")
    parser.add_argument("--solexa",
        action="store_true",
        help="Solexa encoding")

    args = parser.parse_args()

    try:
        qual_trim(args)
    except KeyboardInterrupt: pass
    except Exception as e:
        sys.exit(prog + ": error: " + str(e))
