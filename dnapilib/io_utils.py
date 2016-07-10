"""Functions for IO processes.

"""

import io
import bz2
import gzip
import zipfile
import tarfile
import os.path
import fileinput


def get_file_obj(in_file):
    """Return a file object from an input file.

    """
    if not os.path.exists(in_file) and in_file != "-":
        raise Exception("can't open {}".format(in_file))

    if in_file.find(".tar") > 0:
        if in_file.endswith(".tar.gz"):
            tp = tarfile.open(in_file, "r:gz")
        elif in_file.endswith(".tar"):
            tp = tarfile.open(in_file, "r")
        elif in_file.endswith(".tar.bz2"):
            tp = tarfile.open(in_file, "r:bz2")
        return io.TextIOWrapper(tp)
    elif in_file.endswith(".gz"):
        return gzip.open(in_file, "rt")
    elif in_file.endswith(".zip"):
        zobj = zipfile.ZipFile(in_file)
        zp = zobj.open(zobj.namelist()[0], "r")
        return io.TextIOWrapper(zp)
    elif in_file.endswith(".bz") or in_file.endswith(".bz2"):
        return bz2.BZ2File(in_file, "rt")
    else:
        return fileinput.input(in_file)


def fastq_sequence(fobj):
    """Return sequence lines in FASTQ.

    """
    for i, x in enumerate(fobj):
        if i % 4 == 1: yield x.rstrip()


def fastq_quality(fobj):
    """Return quality score lines in FASTQ.

    """
    for i, x in enumerate(fobj):
        if i % 4 == 3: yield x.rstrip()


def fastq_record(fobj):
    """Return sets of read records in FASTQ.

    """
    record = ""
    for i, x in enumerate(fobj):
        record += x
        if i % 4 == 3:
            yield record
            record = ""
