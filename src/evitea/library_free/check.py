"""The one-site check of a library-free call.

A real insertion junction is new sequence read at single-copy depth, whatever the zygosity. A library-free call has no second
line of evidence (no library, no read pair), so each junction side it rests on is measured:

  depth    T / D <= tolerance. T = the reads that continue the walk's context at the break (the sum of the four letter counts after
           the 30-letter context; the site total when the walk follows the other letter), D = the place's local depth. A junction
           between two reference copies, or inside a collapsed repeat, has many times single-copy depth.
  novelty  the junction (30 host letters + the first 30 inserted letters) must be new to the reference: at most 10 of its
           31-letter windows (offsets 0..20, both strands) occur in the reference, and the junction (up to 60 letters) has no
           ungapped alignment at >= 95% identity to the reference on either strand. A difference between old copies of a
           repeat reads like a junction in short reads but its words are in the reference.

A call passes when it has at least one side and every side it has passes. The sides are those the detection recorded
(``place.process_place``); no read is read again here, only the reference.

Tolerance: the p95 or p99 of the side depth ratios of the library-confirmed family calls of the first pass (calls with an EviTEA
call within 500 bp: the library layer has already confirmed them as insertions, so their ratio shows what a real junction of
this sample looks like at this depth). With no library-confirmed side, or with ``own-calls``, the layer's own two-sided calls
that the reference-path filter keeps and that no rule drops are used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

K = 31
OFFSETS = tuple(range(21))
IDENT = 0.95
WIN_IN_REF_MAX = 10

_COMP = str.maketrans("ACGTN", "TGCAN")
_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _c in enumerate("ACGT"):
    _LUT[ord(_c)] = _i

__all__ = ["RefIndex", "SideInput", "Tolerance", "call_pass", "novelty", "side_inputs", "side_pass", "tolerance"]


def rc(s):
    return s.translate(_COMP)[::-1]


def acgt(s):
    return len(s) > 0 and all(c in "ACGT" for c in s)


def enc(s):
    v = 0
    for ch in s:
        v = (v << 2) | "ACGT".index(ch)
    return v


def to_arr(s):
    return _LUT[np.frombuffer(s.encode("ascii"), dtype=np.uint8)]


def kmer_codes(b, k):
    """2-bit codes of every k-mer of the uint8 array b (A0 C1 G2 T3, other 255): (codes uint64, ok bool)."""
    n = len(b) - k + 1
    if n <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, bool)
    inv = np.concatenate([[0], np.cumsum(b > 3, dtype=np.int64)])
    ok = (inv[k:] - inv[:-k]) == 0
    c = np.zeros(n, dtype=np.uint64)
    two = np.uint64(2)
    for j in range(k):
        c = (c << two) | (b[j:j + n] & 3).astype(np.uint64)
    return c, ok


def _in_sorted(sorted_q, x):
    if len(sorted_q) == 0:
        return np.zeros(len(x), bool)
    idx = np.searchsorted(sorted_q, x)
    idx[idx == len(sorted_q)] = 0
    return sorted_q[idx] == x


class RefIndex:
    """Reference as uint8 arrays; exact 31-mer presence and ungapped >= 95% alignment queries."""

    def __init__(self, genome):
        self.arr = {name: to_arr(s.upper()) for name, s in genome.items()}

    def present_codes(self, codes, k=K):
        """Subset of ``codes`` (ints) that occur as k-mers in the reference forward strand."""
        q = np.unique(np.array(sorted(set(codes)), dtype=np.uint64))
        found = set()
        for b in self.arr.values():
            c, ok = kmer_codes(b, k)
            hit = ok & _in_sorted(q, c)
            if hit.any():
                found.update(int(x) for x in np.unique(c[hit]))
        return found

    def aligned(self, queries, ident=IDENT):
        """queries: qid -> ACGT sequence. {qid: best mismatch count} for the qids with an ungapped alignment of at most
        floor(L * (1 - ident)) mismatches on either strand (reference windows with a non-ACGT letter count it as a mismatch).
        Seed-and-extend with exact non-overlapping pieces, so every such alignment is found."""
        best = {}
        jobs = {}
        info = {}
        for qid, seq in queries.items():
            L = len(seq)
            m = int(L * (1 - ident) + 1e-9)
            plen = L // (m + 1)
            info[qid] = m
            if plen < 8:
                continue
            for strand in (seq, rc(seq)):
                if not acgt(strand):
                    continue
                qarr = to_arr(strand)
                for i in range(m + 1):
                    o = i * plen
                    jobs.setdefault(plen, {}).setdefault(enc(strand[o:o + plen]), []).append((qid, o, qarr))
        ar_cache = {}
        for b in self.arr.values():
            for plen, seedmap in jobs.items():
                S = np.array(sorted(seedmap), dtype=np.uint64)
                pc, ok = kmer_codes(b, plen)
                hp = np.flatnonzero(ok & _in_sorted(S, pc))
                if len(hp) == 0:
                    continue
                hc = pc[hp]
                order = np.argsort(hc, kind="stable")
                hp, hc = hp[order], hc[order]
                for code, items in seedmap.items():
                    lo, hi = np.searchsorted(hc, np.uint64(code), "left"), np.searchsorted(hc, np.uint64(code), "right")
                    if hi <= lo:
                        continue
                    pos = hp[lo:hi]
                    for qid, o, qarr in items:
                        L = len(qarr)
                        ar = ar_cache.setdefault(L, np.arange(L))
                        starts = pos - o
                        starts = starts[(starts >= 0) & (starts + L <= len(b))]
                        for s0 in range(0, len(starts), 200000):
                            st = starts[s0:s0 + 200000]
                            mm = (b[st[:, None] + ar[None, :]] != qarr[None, :]).sum(axis=1)
                            mn = int(mm.min()) if len(mm) else None
                            if mn is not None and (qid not in best or mn < best[qid]):
                                best[qid] = mn
        return {qid: mm for qid, mm in best.items() if mm <= info[qid]}


def junction_windows(junction):
    return [junction[j:j + K] for j in OFFSETS if j + K <= len(junction)]


def novelty(junctions, ref: RefIndex):
    """junctions: id -> junction string. id -> dict(n_win, n_win_in_ref, b2_novel, aligned_mm)."""
    allcodes, wins = set(), {}
    for jid, jn in junctions.items():
        ws = [w for w in junction_windows(jn) if acgt(w)]
        wins[jid] = ws
        for w in ws:
            allcodes.add(enc(w))
            allcodes.add(enc(rc(w)))
    present = ref.present_codes(allcodes) if allcodes else set()
    q = {jid: jn for jid, jn in junctions.items() if acgt(jn) and len(jn) >= 40}
    al = ref.aligned(q) if q else {}
    out = {}
    for jid, jn in junctions.items():
        ws = wins[jid]
        in_ref = [w for w in ws if enc(w) in present or enc(rc(w)) in present]
        out[jid] = dict(n_win=len(ws), n_win_in_ref=len(in_ref), b2_novel=(jid in q and jid not in al),
                        aligned_mm=al.get(jid, ""))
    return out


@dataclass(frozen=True)
class SideInput:
    side: str          # "L" or "R"
    T: int
    ratio: float       # T / D
    junction: str


def side_inputs(row) -> list[SideInput]:
    """The sides the detection recorded for a place row (text fields); empty when it recorded none."""
    out = []
    try:
        D = float(row["D"])
    except (KeyError, TypeError, ValueError):
        return out
    for s in ("L", "R"):
        t = row.get(f"chk{s}_T", "")
        j = row.get(f"chk{s}_junction", "")
        if t not in ("", None) and j:
            T = int(float(t))
            out.append(SideInput(s, T, T / D if D else math.inf, j))
    return out


def side_pass(side: SideInput, nov: dict, tol: float) -> int:
    return int(side.ratio <= tol and nov["n_win"] > 0 and nov["n_win_in_ref"] <= WIN_IN_REF_MAX and nov["b2_novel"])


def call_pass(sides, novs, tol: float) -> int:
    """1 when the call has a side and every side passes, else 0."""
    return int(len(sides) > 0 and all(side_pass(s, novs[s.side], tol) == 1 for s in sides))


@dataclass(frozen=True)
class Tolerance:
    p95: float
    p99: float
    source: str        # "library" | "own-calls"
    n_sides: int

    def value(self, which: str) -> float:
        return self.p99 if which == "p99" else self.p95


def tolerance(library_ratios, own_ratios, source: str) -> Tolerance | None:
    """The tolerance from the side depth ratios of the library-confirmed calls (or, with source ``own-calls`` or no
    library-confirmed side, of the layer's own two-sided calls). None when neither has a side."""
    use, name = (library_ratios, "library") if source == "library" and len(library_ratios) else (own_ratios, "own-calls")
    if not len(use):
        return None
    a = np.asarray(use, dtype=float)
    return Tolerance(p95=float(np.percentile(a, 95)), p99=float(np.percentile(a, 99)), source=name, n_sides=int(len(a)))
