"""Functions to handle k-mer related proceses.

"""

import sys
import re
from operator import itemgetter


def _calc_overlap(x, y, seed):
    """Return an overlapping position between a pair of k-mers.

    """
    if not x or not y:
        return 0
    for m in re.finditer(y[:seed], x):
        overlap = seed
        p = m.end()
        if len(x) == p:
            return overlap
        tail = re.search(x[p:], y)
        if not tail:
            continue
        if tail.start() == seed:
            return tail.end()
    return 0


def filter_kmers(kmers, kmer_len, rate):
    """Return a clean set of k-mers in tuple.

       Filter low-complexity and low-frequency kmers.
    """
    low_comp = [re.compile(base * (kmer_len//2)) for base in "ACGTN"]
    i, x = -1, -1
    while x != len(low_comp):
        i += 1
        x = sum([not p.findall(kmers[i][0]) for p in low_comp])
    max_hits = kmers[i][1]

    clean = []
    total = 0
    for s, n in kmers[i:]:
        if sum([not p.findall(s) for p in low_comp]) != len(low_comp):
            continue
        if float(max_hits)/n > rate:
            break
        clean.append((s, n))
        total += n
    return [(s, round(float(n)/total*100, 4)) for s, n in clean]


def assemble_kmers(kmers, seed):
    """Return assembled k-mers and the frequency in tuple.

       Assemble given k-mers by checking suffix-prefix matches.
    """
    pre_l, new_l = 0, len(kmers)
    while pre_l != new_l:
        pre_l = len(kmers)
        for i in range(pre_l):
            kmer, hits = kmers[i]
            if not hits:
                continue
            max_o, max_j = 0, 0
            for j in range(pre_l):
                if i == j:
                    continue
                if kmers[j][0] in kmer:
                    hits += kmers[j][1]
                    kmers[i] = (kmer, hits)
                    kmers[j] = ('', 0)
                    continue
                overlap = _calc_overlap(kmer, kmers[j][0], seed)
                if overlap > max_o:
                    max_o, max_j = overlap, j
            if max_o > 0:
                kmer += kmers[max_j][0][max_o:]
                hits += kmers[max_j][1]
                kmers[i] = (kmer, hits)
                kmers[max_j] = ('', 0)
        kmers = [k for k in kmers if k != ('', 0)]
        new_l = len(kmers)
    return kmers


def count_kmers(seq_list, kmer_len, sample_num):
    """Return sorted k-mer frequency.

    """
    freq = {}
    for cnt, seq in enumerate(seq_list):
        if cnt == sample_num:
            break
        interval = len(seq) - kmer_len + 1
        for i in range(interval):
            kmer = seq[i : i+kmer_len]
            freq[kmer] = freq.get(kmer, 0) + 1
    return sorted(freq.items(), key=itemgetter(1), reverse=True)
