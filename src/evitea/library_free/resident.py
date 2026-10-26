"""The resident-copy test.

Non-reference insertions land inside old TE copies at about the genome's TE fraction. Differences between old copies of a
family look like junctions in short reads. A library-free call inside an annotated reference TE copy is therefore not reported
when its family's places lie inside such copies more than twice as often as the genome's TE fraction (``family``). The family
condition switches the rule off where TEs cover most of the genome (a per-call rule alone would drop many real
insertions there). ``per-call`` drops every call inside a copy. The test requires the reference's TE annotation. Without one it
is off.

The annotation is read as a RepeatMasker ``.out`` table (lines of >= 15 fields whose first field is a number) or as a BED file
(contig, 0-based start, end, ...). Copies are merged per contig. A position is inside a copy when begin <= position <= end.
"""

from __future__ import annotations

import bisect
from collections import defaultdict
from pathlib import Path

from ..contigs import match_contig_name
from ..errors import InputContractError

__all__ = ["CopyIndex", "FAMILY_ENRICHMENT", "family_inside", "resident_flags"]

FAMILY_ENRICHMENT = 2.0


class CopyIndex:
    def __init__(self, intervals: dict[str, list[tuple[int, int]]]):
        self.starts: dict[str, list[int]] = {}
        self.ends: dict[str, list[int]] = {}
        self.bp = 0
        for c, iv in intervals.items():
            merged: list[list[int]] = []
            for s, e in sorted(iv):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.starts[c] = [m[0] for m in merged]
            self.ends[c] = [m[1] for m in merged]
            self.bp += sum(m[1] - m[0] for m in merged)

    @classmethod
    def from_file(cls, path, contigs=None) -> "CopyIndex":
        """Read a RepeatMasker .out or a BED file. ``contigs`` (names of the reference) lets a bare name such as ``2L`` match ``chr2L``."""
        names = set(contigs or ())
        iv: dict[str, list[tuple[int, int]]] = defaultdict(list)
        try:
            text = Path(path).read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            raise InputContractError(f"cannot read repeat annotation {path}: {exc}") from exc
        for line in text.splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            try:
                if len(f) >= 15 and f[0].isdigit():                  # RepeatMasker .out
                    chrom, s, e = f[4], int(f[5]), int(f[6])
                else:                                               # BED
                    tab = line.split("\t")
                    chrom, s, e = tab[0], int(tab[1]), int(tab[2])
            except (ValueError, IndexError):
                continue
            iv[match_contig_name(chrom, names)].append((s, e))
        if not iv:
            raise InputContractError(f"repeat annotation {path} holds no interval")
        return cls(iv)

    def inside(self, chrom: str, x: int) -> bool:
        st = self.starts.get(chrom)
        if not st:
            return False
        j = bisect.bisect_right(st, x) - 1
        return j >= 0 and self.ends[chrom][j] >= x

    def share(self, total_bp: int) -> float:
        """Annotated bases (merged) over ``total_bp`` bases of the reference."""
        return self.bp / total_bp if total_bp else 0.0


def family_inside(rows, copies: CopyIndex) -> dict[str, float]:
    """cluster -> fraction of the family's family-call rows (any status) whose break lies inside an annotated copy."""
    n: dict[str, int] = defaultdict(int)
    k: dict[str, int] = defaultdict(int)
    for r in rows:
        n[r["cluster"]] += 1
        k[r["cluster"]] += copies.inside(r["chrom"], int(r["bL"]))
    return {c: k[c] / n[c] for c in n}


def resident_flags(rows, copies: CopyIndex, share: float, mode: str):
    """For each family-call row (text fields): (inside, family inside fraction, dropped)."""
    fam = family_inside(rows, copies)
    out = []
    for r in rows:
        ins = copies.inside(r["chrom"], int(r["bL"]))
        f = fam[r["cluster"]]
        drop = ins and (mode == "per-call" or (mode == "family" and f > FAMILY_ENRICHMENT * share))
        out.append((bool(ins), f, bool(drop)))
    return out
