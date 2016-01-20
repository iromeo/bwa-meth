import re

try:
    from future_builtins import zip
except ImportError: # python3
    from builtins import zip


class Bam(object):
    __slots__ = 'read flag chrom pos mapq cigar chrom_mate pos_mate tlen \
            seq qual other'.split()

    def __init__(self, args):
        for attr, value in zip(self.__slots__[:11], args):
            setattr(self, attr, value)
        self.other = args[11:]
        self.flag = int(self.flag)
        self.pos = int(self.pos)
        self.tlen = int(float(self.tlen))

    def __repr__(self):
        return "Bam({chr}:{start}:{read}".format(chr=self.chrom,
                                                 start=self.pos,
                                                 read=self.read)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__[:11]) \
                         + "\t" + "\t".join(self.other)

    def is_first_read(self):
        return bool(self.flag & 0x40)

    def is_second_read(self):
        return bool(self.flag & 0x80)

    def is_plus_read(self):
        return not (self.flag & 0x10)

    def is_minus_read(self):
        return bool(self.flag & 0x10)

    def is_mapped(self):
        return not (self.flag & 0x4)

    def cigs(self):
        cdef str cigar = self.cigar
        if cigar == "*":
            return []

        cdef int length = len(cigar)
        cdef int i = 0
        cdef int offset
        cdef int n
        cdef list acc = []
        while i < length:
            offset = i
            while i < length and cigar[i].isdigit():
                i += 1

            n = int(cigar[offset:i])

            offset = i
            while i < length and not cigar[i].isdigit():
                i += 1

            acc.append((n, cigar[offset:i]))
        return acc

    def left_right_shift(self):
        cdef list cigs = self.cigs()
        cdef int left = self._shift(cigs)
        cigs.reverse()
        cdef int right = self._shift(cigs)
        return left, -right or None

    def _shift(self, list cigs):
        cdef int acc = 0
        cdef int n
        cdef str cig
        for n, cig in cigs:
            if cig == "H":
                acc += n
            elif cig == "M":
                break
        return acc

    @property
    def original_seq(self):
        cdef str tag
        for tag in self.other:
            if tag.startswith("YS:Z"):
                return tag[5:]

        # XXX we must always find it.

    def longest_match(self, patt=re.compile("(\d+)M")):
        return max(int(match.group(1))
                   for match in patt.finditer(self.cigar))
