"""Tiny exact multi-string index double shared by incumbent unit tests."""

from __future__ import annotations

import numpy as np

from evitea.incumbent.msbwt_queries import MsbwtSample


class TinyBWT:
    """Suffix-sorted stand-in for the msbwt-modern3 public query API."""

    def __init__(self, reads: list[str]):
        self.reads = reads
        self.suffixes = sorted((read[i:] + "$", index) for index, read in enumerate(reads) for i in range(len(read) + 1))

    def findIndicesOfStr(self, pattern: str) -> tuple[int, int]:
        keys = [suffix for suffix, _ in self.suffixes]
        lo = next((i for i, key in enumerate(keys) if key >= pattern), len(keys))
        hi = lo
        while hi < len(keys) and keys[hi].startswith(pattern):
            hi += 1
        return lo, hi

    def getSequenceDollarID(self, row: int) -> int:
        return self.suffixes[row][1]

    def recoverString(self, dollar_id: int) -> bytes:
        return ("$" + self.reads[dollar_id]).encode("ascii")


def paired_sample(pairs: list[tuple[str, str]]) -> MsbwtSample:
    """Dollar IDs interleave mates in a shuffled order to exercise about.npy."""

    records = [(0, pid, m1) for pid, (m1, _) in enumerate(pairs)] + [(1, pid, m2) for pid, (_, m2) in enumerate(pairs)]
    records = records[::-1]
    about = np.array([(f, p) for f, p, _ in records], dtype=[("f0", "u1"), ("f1", "<u8")])
    return MsbwtSample("unused", bwt=TinyBWT([seq for _, _, seq in records]), about=about)
