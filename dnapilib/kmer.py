#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015, 2016 Junko Tsuji

# This module contains a set of functions that handle
# K-mer related processes.

import sys, re
from operator import itemgetter

# check overlap between kmers
def _calcOverlap(x, y, seed):
    if not len(x) or not len(y):
        return 0
    for m in re.finditer(y[:seed], x):
        overlap = seed
        p = m.end()
        if len(x) == p: return overlap
        tail = re.search(x[p:], y)
        if not tail: continue
        if tail.start() == seed:
            return tail.end()
    return 0

# filter low-complexity and low-frequency kmers
def filterKmers(kmers, kmer_len, rate):
    cutoff = kmer_len/2
    pa = [re.compile(c * cutoff) for c in "ACGTN"]
    x, i = -1, -1
    while x != len(pa):
        i += 1
        x = sum([not p.findall(kmers[i][0]) for p in pa])
    max_hits = float(kmers[i][1])

    clean = []
    for s, n in kmers[i:]:
        x = sum([not p.findall(s) for p in pa])
        if x != len(pa): continue
        if max_hits/n > rate: break
        clean.append((s, n))
    return clean

# assemble kmers with exhausitve way
def assembleKmers(kmers, seed):
    pre_l, new_l = 0, len(kmers)
    while pre_l != new_l:
        pre_l = len(kmers)
        for i in xrange(pre_l):
            kmer, hits = kmers[i]
            if not hits: continue
            max_o, max_j = 0, 0
            for j in xrange(pre_l):
                if i == j: continue
                o = _calcOverlap(kmer, kmers[j][0], seed)
                if o > max_o: max_o, max_j = o, j
            if max_o > 0:
                kmer += kmers[max_j][0][max_o:]
                hits += kmers[max_j][1]
                kmers[i] = (kmer, hits)
                kmers[max_j] = ('', 0)
        kmers = [k for k in kmers if k != ('', 0)]
        new_l = len(kmers)
    return kmers

# count kmers
def countKmers(fIn, kmer_len, sample_num):
    freq = {}
    for cnt, seq in enumerate(fIn):
        if cnt == sample_num: break
        l = len(seq)
        for i in xrange(l - kmer_len + 1):
            kmer = seq[i:i+kmer_len]
            freq[kmer] = freq.get(kmer, 0) + 1
    freq = sorted(freq.items(), key=itemgetter(1), reverse=True)
    return freq
