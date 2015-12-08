#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015, 2016 Junko Tsuji

import sys, re
from io import fileObject, fastqSequence

def clipAdapter(fqObj, aseed, tm5, tm3, clipMin, clipMax):
    seedLen = len(aseed)
    pp = re.compile("(.*)" + aseed, re.IGNORECASE)
    for seq in fastqSequence(fqObj):
        match = pp.search(seq)
        if not match: continue
        end = match.end() - seedLen
        clippedSeq = seq[tm5 : end-tm3]
        l = len(clippedSeq)
        if clipMin <= l and l <= clipMax:
            yield clippedSeq

def toFasta(fq, fa, aseed, tm5, tm3, clipMin, clipMax):
    fqObj = fileObject(fq)
    if aseed == "RAW_INPUT":
        iterator = fastqSequence(fqObj)
    else:
        iterator = clipAdapter(fqObj, aseed, tm5, tm3, clipMin, clipMax)
    fas = {}
    cleanReadCount = 0
    for seq in iterator:
        fas[seq] = fas.get(seq, 0) + 1  # tally sequences
    faObj = open(fa, "w")
    for seq, cnt in fas.items():
        cleanReadCount += cnt
        faObj.write(">%s_%d\n%s\n" % (seq, cnt, seq))
    faObj.close()
    fqObj.close()
    return cleanReadCount
