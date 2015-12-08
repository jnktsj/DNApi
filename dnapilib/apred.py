#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015, 2016 Junko Tsuji

# This module contains the core function of adapter prediction.

from operator import itemgetter

from io import fileObject, fastqSequence
from kmer import countKmers, filterKmers, assembleKmers

def adapterPrediction(fq, ratio, kmer_len, sample_num):
    fqObj = fileObject(fq)
    fqIn = fastqSequence(fqObj)
    freq = countKmers(fqIn, kmer_len, sample_num)
    clean = filterKmers(freq, kmer_len, ratio)
    assembl = sorted(assembleKmers(clean, kmer_len/2),
                     key=itemgetter(1), reverse=True)
    fqObj.close()
    return assembl
