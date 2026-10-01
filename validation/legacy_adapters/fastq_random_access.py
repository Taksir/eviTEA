"""Validation-only random access to SRA-style FASTQ records by ordinal.

Requires headers of the form ``@<run>.<n> ...`` with ``n = ordinal + 1`` (as in
the AKA017 SRA FASTQs).  Used only to fetch small read subsets for stage-level
differential tests without streaming a 12-GB file.
"""

from __future__ import annotations

import os
from pathlib import Path


class SraFastqIndex:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.size = self.path.stat().st_size
        self.handle = self.path.open("rb")

    def _record_at(self, offset: int) -> tuple[int, int, bytes]:
        """First complete record at or after ``offset``: (ordinal, start, sequence)."""

        self.handle.seek(max(0, offset - 1))
        if offset > 0:
            self.handle.readline()
        while True:
            start = self.handle.tell()
            line = self.handle.readline()
            if not line:
                raise EOFError
            if line.startswith(b"@"):
                seq = self.handle.readline()
                plus = self.handle.readline()
                qual = self.handle.readline()
                if plus.startswith(b"+") and len(qual.rstrip(b"\r\n")) == len(seq.rstrip(b"\r\n")):
                    token = line[1:].split()[0].decode()
                    return int(token.rsplit(".", 1)[1]) - 1, start, seq.strip()
                self.handle.seek(start)
                self.handle.readline()

    def sequence(self, ordinal: int) -> str:
        lo, hi = 0, self.size
        while hi - lo > 4096:
            mid = (lo + hi) // 2
            found, _, _ = self._record_at(mid)
            if found <= ordinal:
                lo = mid
            else:
                hi = mid
        self.handle.seek(lo)
        if lo:
            self.handle.readline()
        offset = self.handle.tell()
        while True:
            found, start, seq = self._record_at(offset)
            if found == ordinal:
                return seq.decode().upper()
            if found > ordinal:
                raise KeyError(ordinal)
            offset = self.handle.tell()


def pair_reads(read1: str | Path, read2: str | Path, pair_ids) -> dict[int, tuple[str, str]]:
    a, b = SraFastqIndex(read1), SraFastqIndex(read2)
    return {pid: (a.sequence(pid), b.sequence(pid)) for pid in sorted(set(pair_ids))}
