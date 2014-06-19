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
        for i in range(min(len(seq) - offset, width)):
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
            chunks.append(line)


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
    except Exception:
        os.unlink(out_fasta)
        raise
    finally:
        ih.close()
        oh.close()

    return out_fasta


re_mate_tag = re.compile(r"(?:_R/)[12]")

def _process_read(fh, out, int mate):
    cdef bytes name = next(fh).split(None, 2)[0]
    name = re_mate_tag.sub("", name)
    cdef bytes seq = next(fh).upper().rstrip()
    _ = next(fh)  # Skip.
    cdef bytes qual = next(fh)

    cdef bytes a
    cdef bytes b
    if mate == 1:
        a, b = b"C", b"T"
    else:
        a, b = b"G", b"A"

    out.write(os.linesep.join((name + " YS:Z:" + seq + "\tYC:Z:" + a + b,
                               seq.replace(a, b),
                               "+", qual)))
    return len(seq) < 80


def convert_reads(bytes fq1s, bytes fq2s, out=sys.stdout):
    cdef long long lt80 = 0
    cdef list files = list(zip(fq1s.split(","), fq2s.split(",")))
    for i, (fq1, fq2) in enumerate(files, 1):
        msg = ("converting reads in {0},{1} ({2:02d} files left)"
               .format(fq1, fq2, len(files) - i))
        print(msg, file=sys.stderr)
        fq1 = nopen(fq1)
        fq2 = nopen(fq2) if fq2 != "NA" else None

        try:
            while True:
                lt80 += _process_read(fq1, out, 1)
                if fq2 is not None:
                    lt80 += _process_read(fq2, out, 2)
        except (StopIteration, IOError):
            pass

    out.flush()
    out.close()
    if lt80 > 50:
        print("WARNING: {0} reads with length < 80".format(lt80),
              file=sys.stderr)
        print("       : this program is designed for long reads",
              file=sys.stderr)
