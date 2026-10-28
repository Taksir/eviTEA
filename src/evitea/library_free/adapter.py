"""Adapter test: a place whose inserted segment holds a sequencing-adapter core is not a junction.

A read that runs past the end of a short insert continues into the adapter. At a unique place this is new sequence at
single-copy depth, which is what the one-site check accepts as a junction, and the few reads that carry it look like a
second allele at a low share. No TE family contains an adapter core. The test runs on the scanned places of every pass
before they are grouped, so these places enter no family, no family count and no family record. Their stage becomes
``S4_adapter``.

Cores: AGATCGGAAGAGC (the start of both TruSeq adapters) and CTGTCTCTTATACACATCT (Nextera), with their reverse complements,
anywhere in the left or right segment. The walk follows the reference while the next adapter letters match it by chance, so
the junction can sit up to SHIFT letters inside the adapter. The junction end of a segment (the start of the left segment,
the start of the reverse complement of the right one) that begins with a core less its first 1 to SHIFT letters is therefore
an adapter segment too. A shift of k letters happens with probability 4^-k; 4 letters cover 99.6%.

The test reads only the place records, not the scan, so saved scan pieces stay valid with it on or off.
"""

from __future__ import annotations

from collections.abc import Iterable

__all__ = ["CORES", "SHIFT", "STAGE", "TESTED_STAGES", "adapter_segment", "mark_adapter_places"]

CORES = ("AGATCGGAAGAGC", "CTGTCTCTTATACACATCT")
SHIFT = 4
STAGE = "S4_adapter"
TESTED_STAGES = ("S4_pass", "S4_reject")

_COMP = str.maketrans("ACGTN", "TGCAN")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


_ANYWHERE = CORES + tuple(_rc(c) for c in CORES)
_SHIFTED = tuple(c[k:] for c in CORES for k in range(1, SHIFT + 1))


def adapter_segment(side: str, seg: str) -> bool:
    """Does the segment of ``side`` (L or R, forward orientation: the left one starts at the junction, the right one ends
    there) hold an adapter core anywhere, or a shifted core at its junction end?"""
    s = (seg or "").upper()
    if not s:
        return False
    if any(c in s for c in _ANYWHERE):
        return True
    return (s if side == "L" else _rc(s)).startswith(_SHIFTED)


def mark_adapter_places(places: Iterable[dict]) -> int:
    """Set the stage of every S4_pass or S4_reject place with an adapter segment on either side to ``S4_adapter``; return
    how many."""
    n = 0
    for p in places:
        if p.get("stage") in TESTED_STAGES and (adapter_segment("L", p.get("lseg") or "")
                                                 or adapter_segment("R", p.get("rseg") or "")):
            p["stage"] = STAGE
            n += 1
    return n
