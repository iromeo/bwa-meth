__all__ = ["convert_fasta"]

import sys
import os
from itertools import groupby

from toolshed import nopen, is_newer_b


def wrap(text, width=100): # much faster than textwrap
    for i in range(0, len(text), width):
        yield text[i:i + width]


def fasta_iter(fasta_name):
    fh = nopen(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        yield header, "".join(s.strip() for s in next(faiter)).upper()


def convert_fasta(ref_fasta, just_name=False):
    out_fasta = ref_fasta + ".bwameth.c2t"
    if just_name:
        return out_fasta

    msg = "c2t in %s to %s" % (ref_fasta, out_fasta)
    if is_newer_b(ref_fasta, out_fasta):
        sys.stderr.write("already converted: %s\n" % msg)
        return out_fasta
    sys.stderr.write("converting %s\n" % msg)
    try:
        fh = open(out_fasta, "w")
        for header, seq in fasta_iter(ref_fasta):
            ########### Reverse ######################
            fh.write(">r%s\n" % header)

            ########### Forward ######################
            fh.write(">f%s\n" % header)
            for line in wrap(seq.replace("C", "T")):
                fh.write(line + '\n')
        fh.close()
    finally:
        fh.close()
        os.unlink(out_fasta)

    return out_fasta
