__all__ = ["convert_fasta"]

import sys
import os
from itertools import groupby

from toolshed import nopen, is_newer_b


def fasta_write(fh, bytes text, int width=100):
    cdef int i
    cdef int offset
    for offset in range(0, len(text), width):
        for i in range(width):
            fh.write(text[i + offset])
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


def convert_fasta(ref_fasta, just_name=False):
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
            oh.write(">r%s\n" % header)  # Reverse
            fasta_write(oh, seq.replace("G", "A"))
            oh.write(">f%s\n" % header)  # Forward
            fasta_write(oh, seq.replace("C", "T"))
    finally:
        ih.close()
        oh.close()
        os.unlink(out_fasta)

    return out_fasta
