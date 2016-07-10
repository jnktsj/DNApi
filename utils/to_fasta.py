#!/usr/bin/env python3

"""Clip off adapters from reads in FASTQ by searching kmer
   exact prefix matches. There is an option '-s' that tries
   to find prefix matches in length of (k + 1) nucleotides
   with 1 mismatch if perfect prefix matches are not found.
   After adapter removal, the program tallies clean reads
   and writes non-redundant reads with the counts in FASTA.

"""


import sys
import re
import os.path
import signal
from operator import itemgetter
from argparse import ArgumentParser


cur = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(cur))
from dnapilib.io_utils import get_file_obj, fastq_sequence


def make_regex(a_seq, a_len, is_3prime, sensitive):
    if not a_seq:
        return None, None

    cutoff = a_len
    if sensitive:
        cutoff = a_len + 1

    if len(a_seq) < cutoff:
        message = "input adapters longer than {} nt"
        raise Exception(message.format(cutoff-1))

    a_seq = a_seq.upper().replace('U', 'T')
    if is_3prime:
        pat = "(.*)"
        p_seq, m_seq = a_seq[:a_len], a_seq[:a_len+1]
        beg, end = 0, a_len
    else:
        pat = "(.*?)"
        p_seq, m_seq = a_seq[-a_len:], a_seq[-(a_len+1):]
        beg, end = 1, a_len+1

    pp = re.compile(pat+p_seq, re.IGNORECASE)
    if not sensitive:
        return pp, []

    mp = []
    seed = list(m_seq)
    for i in range(beg, end):
        ps = seed[:]
        ps[i] = '.'
        x = re.compile(pat+''.join(ps), re.IGNORECASE)
        mp.append(x)
    return pp, mp


def match_adapters(seq, pp, mps, sensitive, pi=0, mi=0, l=0):
    if not pp:
        return l, '*'

    p = pp.search(seq)
    if p:
        return (p.end()-pi), '0'

    if not sensitive:
        return l, '*'
    else:
        for mp in mps:
            m = mp.search(seq)
            if m:
                return (m.end()-mi), '1'
    return l, '*'


def to_fasta(args):
    if args.m <= 0:
        raise Exception("bad value: -m")
    if args.x <= 0:
        raise Exception("bad value: -x")
    if args.m == args.x:
        raise Exception("bad read length cutoff range")
    if not args.f and not args.b:
        raise Exception("input adapter sequence")
    if args.f == args.b:
        raise Exception("5' and 3' adapters are same sequences")
    if args.seed_5p <= 0:
        raise Exception("bad value: --seed-5p")
    if args.seed_3p <= 0:
        raise Exception("bad value: --seed-3p")
    if args.trim_3p < 0:
        raise Exception("input positive value for 3'trimming")
    if args.trim_5p < 0:
        raise Exception("input positive value for 5'trimming")

    f_seq, f_len = args.f, args.seed_5p
    b_seq, b_len = args.b, args.seed_3p
    if not f_seq and b_seq:
        req = lambda x, y: y != '*' or args.a
    elif f_seq and b_seq:
        if args.B:
            req = lambda x, y: x != '*' and y != '*' or args.a
        else:
            req = lambda x, y: x != '*' or y != '*' or args.a
    elif f_seq and not b_seq:
        req = lambda x, y: x != '*' or args.a

    f_pp, f_mp = make_regex(f_seq, f_len, False, args.s)
    b_pp, b_mp = make_regex(b_seq, b_len, True,  args.s)

    fas = {}
    for seq in fastq_sequence(get_file_obj(args.FASTQ)):
        seq_len = len(seq)
        f_i, f_mm = match_adapters(seq, f_pp, f_mp, args.s)
        b_i, b_mm = match_adapters(seq, b_pp, b_mp, args.s,
                                   b_len, b_len+1, seq_len)
        ins = seq[f_i+args.trim_5p : b_i-args.trim_3p]
        ins_len = len(ins)
        if req(f_mm, b_mm) and (args.m <= ins_len and ins_len <= args.x):
            fas[ins] = fas.get(ins, 0) + 1

    fas = sorted(fas.items(), key=itemgetter(0))
    for seq, cnt in fas:
        print(">{0}_{1}\n{0}".format(seq, cnt))


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    prog = os.path.basename(sys.argv[0])
    if sys.version_info.major <= 2:
        raise ValueError("{} requires python version 3 or higher".format(prog))

    parser = ArgumentParser(
        description="Remove adapters and collapse reads from FASTQ to FASTA")
    parser.add_argument("FASTQ",
        type=str,
        help="including stdin or compressed file {zip,gz,tar,bz}")
    parser.add_argument("-3",
        metavar="SEQ",
        dest="b",
        help="3'adapter sequence")
    parser.add_argument("-5",
        metavar="SEQ",
        dest="f",
        help="5'adapter sequence")
    parser.add_argument("--trim-5p",
        metavar="BP",
        type=int, default=0,
        help="trim specified number of bases from 5'ends")
    parser.add_argument("--trim-3p",
        metavar="BP",
        type=int, default=0,
        help="trim specified number of bases from 3'ends")
    parser.add_argument("--seed-5p",
        metavar="BP",
        type=int, default=7,
        help="5' adapter match length in bp (default: %(default)s)")
    parser.add_argument("--seed-3p",
        metavar="BP",
        type=int,default=7,
        help="3' adapter match length in bp (default: %(default)s)")
    parser.add_argument("-m",
        metavar="BP",
        type=int, default=16,
        help="minimum read length in bp (default: %(default)s)")
    parser.add_argument("-x",
        metavar="BP",
        type=int, default=36,
        help="maximum read length in bp (default: %(default)s)")
    parser.add_argument("-s",
        action="store_true",
        help="sensitive adapter search with 1 mismatch (default: off)")
    parser.add_argument("-B",
        action="store_true",
        help="only print the reads with both 5' and 3' adapter matches")
    parser.add_argument("-a",
        action="store_true",
        help="print all reads with and without adapter matches if the "
             "reads are in the range specified with '-m' and '-x'")

    args = parser.parse_args()

    try:
        to_fasta(args)
    except KeyboardInterrupt: pass
    except Exception as e:
        sys.exit(prog + ": error: " + str(e))
