#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015, 2016 Junko Tsuji

# This module contains mainly for IO processes.

import os.path, fileinput, tarfile, gzip, zipfile, bz2

# return a file object
def fileObject(f):
    if not os.path.exists(f) and f != "-":
        raise Exception("can't open %s" % f)
    fIn = None
    if f.find(".tar") > 0:
        if f.endswith(".tar"):
            fIn = tarfile.open(f, "r")
        elif f.endswith(".tar.gz"):
            fIn = tarfile.open(f, "r:gz")
        elif f.endswith(".tar.bz2"):
            fIn = tarfile.open(f, "r:bz2")
    elif f.endswith(".gz"):
        fIn = gzip.open(f, "r")
    elif f.endswith(".zip"):
        zo = zipfile.ZipFile(f)
        fIn = zo.open(zo.namelist()[0], "r")
    elif f.endswith(".bz") or f.endswith(".bz2"):
        fIn = bz2.BZ2File(f, "r")
    else:
        fIn = fileinput.input(f)
    return fIn

# return a sequence in FASTQ
def fastqSequence(fobj):
    for i, x in enumerate(fobj):
        if i % 4 == 1: yield x.rstrip()

# return a quality score for a sequence in FASTQ
def fastqQuality(fobj):
    for i, x in enumerate(fobj):
        if i % 4 == 3: yield x.rstrip()

# return a record in FASTQ
def fastqRecord(fobj):
    record = ""
    for i, x in enumerate(fobj):
        record += x
        if i % 4 == 3:
            yield record
            record = ""
