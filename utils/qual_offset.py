#!/usr/bin/env python3

"""Guess fastq quality encoding offset by checking the range
    of the ASCII-encoded quality scores in FASTQ.

"""

import sys
import os.path
import signal
import subprocess
from argparse import ArgumentParser

cur = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(cur))
from dnapilib.io_utils import get_file_obj, fastq_quality


def guess_qual_offset(args):
    platform = [ ('Sanger/Illumina-1.8+', 33,  76, 33),
                 ('Illumina-1.5+', 67, 104, 64),
                 ('Illumina-1.3+', 54, 104, 64),
                 ('Solexa', 59, 104, 64) ]
    q_chars = set()
    fastqs = fastq_quality(get_file_obj(args.FASTQ))
    sample_num = 50000
    for i, quality in enumerate(fastqs):
        if i == sample_num:
            break
        q_chars = q_chars.union(quality)

    q_int = sorted(list(map(ord, q_chars)))
    if len(q_int) <= 1:
        raise Exception("unknown quality encoding")
    for pl in platform:
        if pl[1] <= q_int[1] and q_int[-2] <= pl[2]:
            return "{}:base={}".format(pl[0], pl[3])

    raise Exception("unknown quality encoding")


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    prog = os.path.basename(sys.argv[0])
    if sys.version_info.major <= 2:
        raise ValueError("{} requires python version 3 or higher".format(prog))

    parser = ArgumentParser(
                 description="Estimate quality score encoding")
    parser.add_argument("FASTQ",
        type=str,
        help="including stdin or compressed file {zip,gz,tar,bz}")
    args = parser.parse_args()

    try:
        print(guess_qual_offset(args))
    except KeyboardInterrupt: pass
    except Exception as e:
        sys.exit(prog + ": error: " + str(e))
