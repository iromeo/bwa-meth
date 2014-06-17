from __future__ import print_function

__all__ = ["convert_fasta", "convert_reads"]

import re
import sys
import os
from itertools import repeat
from future_builtins import zip

from toolshed import nopen, is_newer_b


def fasta_write(fh, bytes header, bytes seq, int width=100):
    cdef:
        int i
        int offset
    fh.write(">")
    fh.write(header)
    fh.write(os.linesep)
    for offset in range(0, len(seq), width):
        for i in range(width):
            fh.write(seq[i + offset])
        fh.write(os.linesep)


def fasta_iterread(fh):
    cdef bytes header = None
    cdef list chunks = []
    for line in fh:
        if line[0] == ">":
            if header is not None:
                yield header, "".join(chunks).upper()
                del chunks[:]

            header = line.strip(" \t\r\n>")
        else:
            chunks.add(line)


def convert_fasta(bytes ref_fasta, just_name=False):
    out_fasta = ref_fasta + ".bwameth.c2t"
    if just_name:
        return out_fasta

    msg = "c2t in %s to %s" % (ref_fasta, out_fasta)
    if is_newer_b(ref_fasta, out_fasta):
        sys.stderr.write("already converted: %s\n" % msg)
        return out_fasta
    sys.stderr.write("converting %s\n" % msg)

    ih = nopen(ref_fasta, "rb")
    oh = open(out_fasta, "wb")
    try:
        for header, seq in fasta_iterread(ih):
            fasta_write(oh, "r" + header, seq.replace("G", "A"))  # Reverse.
            fasta_write(oh, "f" + header, seq.replace("C", "T"))  # Forward.
    finally:
        ih.close()
        oh.close()
        os.unlink(out_fasta)

    return out_fasta


re_mate_tag = re.compile(r"(?:_R/)[12]")

cdef int _process_read(fh, out, int mate):
    cdef bytes name = fh.readline().split(maxsplit=2)[0]
    name = re_mate_tag.sub("", name)
    cdef bytes seq = fh.readline().upper().rstrip()
    fh.readline()
    cdef bytes qual = fh.readline()

    cdef char a
    cdef char b
    if mate == 1:
        a, b = "C", "T"
    else:
        a, b = "G", "A"

    out.write(os.linesep.join((name + " YS:Z:" + seq + "\tYC:Z:" + a + b,
                               seq.replace(a, b),
                               "+", qual)))
    return len(seq) < 80


cdef int is_eof(fh):
    return fh.tell() == os.fstat(fh.fileno()).st_size


def convert_reads(list fq1s, list fq2s, out=sys.stdout):
    cdef long long lt80 = 0
    for fq1, fq2 in zip(fq1s.split(","), fq2s.split(",")):
        sys.stderr.write("converting reads in %s,%s\n" % (fq1, fq2))
        fq1 = nopen(fq1)
        fq2 = nopen(fq2) if fq2 != "NA" else None

        while not is_eof(fq1):  # This is enough.
            lt80 += _process_read(fq1, out, 1)
            if fq2 is not None:
                lt80 += _process_read(fq2, out, 2)

    out.flush()
    out.close()
    if lt80 > 50:
        sys.stderr.write("WARNING: %i reads with length < 80\n" % lt80)
        sys.stderr.write("       : this program is designed for long reads\n")
