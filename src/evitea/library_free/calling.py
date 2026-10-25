"""Family calls of the library-free detector: the reference-path filter, the family-identity filter and the merge with EviTEA's calls.

Everything here reads the places table as text fields, exactly as the research
scorers read ``P4_PLACES.tsv`` (``p4_library_free/v2/code/v2_lib.py``,
``precision_diag/code/precb_run.py`` with ``v2_diag/code/v2diag_lib.py``).
"""

from __future__ import annotations

import bisect
import csv
import io
import math
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import astuple, dataclass, fields

import numpy as np

from . import NEAR_EVITEA_BP, RULE_C_MAX, T_E

__all__ = [
    "ADDED_CALL_COLUMNS",
    "CallIndex",
    "FAMILY_CALL_COLUMNS",
    "FamilyCall",
    "added_calls_tsv_bytes",
    "cluster_members",
    "e_value",
    "family_calls",
    "family_calls_tsv_bytes",
    "family_similarity",
    "feat_c",
    "is_family_call",
    "pair_identity",
    "rule_c_dropped",
    "sw",
    "v2_dropped",
]

# family-identity filter alignment (v2diag_lib): Smith-Waterman, match +2, mismatch -3, linear gap 5; N never matches
MATCH, MISM, GAP = 2, -3, 5
SEG_MIN = 60
C_MINCOLS, C_MAXPAIRS, C_SEED = 50, 3000, 20261003
_RC = str.maketrans("ACGTN", "TGCAN")


# ------------------------------------------------------------------ reference-path filter (v2_lib)
def _num(value) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        v = float(value)
    else:
        text = str(value).strip()
        if text == "":
            return None
        try:
            v = float(text)
        except ValueError:
            return None
    return v if math.isfinite(v) else None


def e_value(fracL, fracR, rstat) -> float:
    """Reference-path fraction at the break: min(fracL, fracR) if the right walk broke and fracR is a number, else fracL."""

    fl, fr = _num(fracL), _num(fracR)
    if rstat == "break" and fr is not None:
        return float("nan") if fl is None else min(fl, fr)
    return fl if fl is not None else float("nan")


def v2_dropped(e, t: float = T_E) -> bool:
    """Return True when the reference-path filter drops a call (e >= T_E). A value that is not a number is never dropped."""

    if e is None or (isinstance(e, float) and math.isnan(e)):
        return False
    return e >= t


def is_family_call(p: Mapping[str, str]) -> bool:
    """Return True for a family call of a valid-TSD family: stage S4_pass, not a duplicate, in a cluster, fam_vtsd 1."""

    return p["stage"] == "S4_pass" and p["dup"] != "1" and p["cluster"] != "" and p["fam_vtsd"] == "1"


def rule_c_dropped(c: float) -> bool:
    """Return True when the family-identity filter drops a call whose family similarity c <= RULE_C_MAX. NaN is never dropped."""

    return not math.isnan(c) and c <= RULE_C_MAX


# ------------------------------------------------------------------ family-identity filter: within-family similarity (v2diag_lib)
def rc(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _clean(q: str, t: str) -> tuple[str, str]:
    q = "".join(c if c in "ACGT" else "N" for c in q.upper())
    t = "".join(c if c in "ACGT" else "X" for c in t.upper())
    return q, t


def _stats(ops, q_beg, t_beg, qlen, score) -> dict:
    m = sum(n for n, o in ops if o == "=")
    cols = sum(n for n, o in ops)
    qspan = sum(n for n, o in ops if o in "=XI")
    tspan = sum(n for n, o in ops if o in "=XD")
    return dict(score=score, qcov=qspan / qlen if qlen else 0.0, ident=m / cols if cols else 0.0, cols=cols, matches=m,
                q_beg=q_beg, q_end=q_beg + qspan, t_beg=t_beg, t_end=t_beg + tspan)


_NULL = dict(score=0, qcov=0.0, ident=0.0, cols=0, matches=0, q_beg=0, q_end=0, t_beg=0, t_end=0)
_PARASAIL = None


def _parasail():
    global _PARASAIL
    if _PARASAIL is None:
        import parasail

        _PARASAIL = (parasail, parasail.matrix_create("ACGTNX", MATCH, MISM))
    return _PARASAIL


def sw(q: str, t: str) -> dict:
    """Return a parasail local alignment (striped, 16-bit), summarized as in the research code."""

    if not q or not t:
        return dict(_NULL)
    q, t = _clean(q, t)
    ps, mat = _parasail()
    r = ps.sw_trace_striped_16(q, t, GAP, GAP, mat)
    if r.score <= 0:
        return dict(_NULL)
    cig = r.cigar
    raw = cig.decode.decode() if isinstance(cig.decode, bytes) else cig.decode
    ops, num = [], ""
    for ch in raw:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num), ch))
            num = ""
    qb, tb = cig.beg_query, cig.beg_ref
    while ops and ops[0][1] in "ID":  # parasail reports the start offset as leading gaps
        n, o = ops.pop(0)
        if o == "I":
            qb += n
        else:
            tb += n
    while ops and ops[-1][1] in "ID":
        ops.pop()
    return _stats(ops, qb, tb, len(q), r.score)


def pair_identity(segs_a, segs_b, aligner: Callable | None = None):
    """Return the identity of the best-scoring local alignment over all segment pairs and both strands; None below 50 columns."""

    aligner = aligner or sw
    best = None
    for a in segs_a:
        for b in segs_b:
            for bb in (b, rc(b)):
                r = aligner(a, bb)
                if best is None or r["score"] > best["score"]:
                    best = r
    if best is None or best["cols"] < C_MINCOLS:
        return None
    return best["ident"]


def feat_c(member_segs, aligner: Callable | None = None, rng_seed: int = C_SEED):
    """Return (c, pairs, used): the median pair identity over a family's members (a seeded sample of 3,000 pairs above that)."""

    n = len(member_segs)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if len(pairs) > C_MAXPAIRS:
        rng = np.random.default_rng(rng_seed)
        sel = rng.choice(len(pairs), C_MAXPAIRS, replace=False)
        pairs = [pairs[k] for k in sorted(sel.tolist())]
    vals = []
    for i, j in pairs:
        if member_segs[i] and member_segs[j]:
            v = pair_identity(member_segs[i], member_segs[j], aligner)
            if v is not None:
                vals.append(v)
    return (float(np.median(vals)) if vals else float("nan")), len(pairs), len(vals)


def cluster_members(places: Sequence[Mapping[str, str]], segments: Mapping) -> dict[str, list[list[str]]]:
    """Return cluster -> member segment lists (segments of >= 60 letters), members in table order (grouped places only)."""

    out: dict[str, list[list[str]]] = defaultdict(list)
    for p in places:
        if not (p["stage"] == "S4_pass" and p["dup"] != "1" and p["cluster"] != ""):
            continue
        d = segments.get((p["chrom"], p["bL"], p["cluster"]), {})
        if len(d.get("lseg", "")) != int(p["lseg_len"]) or len(d.get("rseg", "")) != int(p["rseg_len"]):
            raise ValueError(f"segments do not match the places table at {p['chrom']}:{p['bL']}")
        out[p["cluster"]].append([d[s] for s in ("lseg", "rseg") if len(d.get(s, "")) >= SEG_MIN])
    return dict(out)


def family_similarity(places, segments, clusters: Iterable[str] | None = None, aligner: Callable | None = None):
    """Return cluster -> (c, pairs, used) for the given clusters (default: clusters holding a valid-TSD family call)."""

    members = cluster_members(places, segments)
    wanted = sorted({p["cluster"] for p in places if is_family_call(p)} if clusters is None else set(clusters), key=int)
    return {c: feat_c(members.get(c, []), aligner) for c in wanted}


# ------------------------------------------------------------------ merge with EviTEA's calls
class CallIndex:
    """EviTEA calls per contig. ``near`` is true when some call [start, end] lies within ``tol`` bases of x (0 inside)."""

    def __init__(self, calls: Mapping[str, Iterable[tuple[int, int]]], tol: int = NEAR_EVITEA_BP):
        self.tol = tol
        self.by_chrom = {ch: sorted((int(s), int(e)) for s, e in v) for ch, v in calls.items()}
        self.max_len = {ch: max((e - s for s, e in v), default=0) for ch, v in self.by_chrom.items()}
        self.starts = {ch: [s for s, _ in v] for ch, v in self.by_chrom.items()}

    @classmethod
    def from_native(cls, rows: Iterable[Mapping[str, str]], tol: int = NEAR_EVITEA_BP) -> "CallIndex":
        calls: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for r in rows:
            calls[r["chr"]].append((int(r["start"]), int(r["end"])))
        return cls(calls, tol)

    def near(self, chrom: str, x: int) -> bool:
        v = self.by_chrom.get(chrom)
        if not v:
            return False
        hi = bisect.bisect_right(self.starts[chrom], x + self.tol)
        lo = bisect.bisect_left(self.starts[chrom], x - self.tol - self.max_len[chrom])
        for s, e in v[lo:hi]:
            d = 0 if s <= x <= e else (s - x if x < s else x - e)
            if d <= self.tol:
                return True
        return False


@dataclass(frozen=True, slots=True)
class FamilyCall:
    place_idx: int
    chrom: str
    bL: int
    bR: str
    tsd: str
    cluster: str
    fracL: str
    fracR: str
    rstat: str
    e: float
    c: float
    c_pairs: int
    c_used: int
    v2_dropped: bool
    rule_c_dropped: bool
    near_evitea: bool
    status: str


FAMILY_CALL_COLUMNS = tuple(f.name for f in fields(FamilyCall))
ADDED_CALL_COLUMNS = ("event_id", "chr", "start", "end", "family", "cluster", "bL", "bR", "tsd", "e", "c")


def family_calls(places: Sequence[Mapping[str, str]], evitea: CallIndex, c_by_cluster: Mapping[str, tuple],
                 *, rule_c: bool) -> tuple[FamilyCall, ...]:
    """Return one row per valid-TSD family call, in table order. The status is DROPPED_V2 (reference-path filter), then
    DROPPED_RULE_C (family-identity filter, when on), then NEAR_EVITEA (an EviTEA call within 500 bp of bL), else ADDED."""

    out = []
    for i, p in enumerate(places):
        if not is_family_call(p):
            continue
        e = e_value(p["fracL"], p["fracR"], p["rstat"])
        c, pairs, used = c_by_cluster.get(p["cluster"], (float("nan"), 0, 0))
        dv, dc = v2_dropped(e), rule_c_dropped(c)
        x = int(p["bL"])
        near = evitea.near(p["chrom"], x)
        if dv:
            status = "DROPPED_V2"
        elif rule_c and dc:
            status = "DROPPED_RULE_C"
        elif near:
            status = "NEAR_EVITEA"
        else:
            status = "ADDED"
        out.append(FamilyCall(i, p["chrom"], x, p["bR"], p["tsd"], p["cluster"], p["fracL"], p["fracR"], p["rstat"],
                              e, c, pairs, used, dv, dc, near, status))
    return tuple(out)


def _tsv(header, rows) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return handle.getvalue().encode("utf-8")


def _cell(value):
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, float):
        return "" if math.isnan(value) else repr(value)
    return value


def family_calls_tsv_bytes(rows: Iterable[FamilyCall]) -> bytes:
    return _tsv(FAMILY_CALL_COLUMNS, ([_cell(v) for v in astuple(r)] for r in rows))


def added_calls_tsv_bytes(rows: Iterable[FamilyCall], *, sample_id: str = "") -> bytes:
    """Return the added calls. The call interval is [bL, bL + 1), as in the scored union of the research code."""

    del sample_id  # reserved: ids are made unique per contig and position
    seen: dict[str, int] = defaultdict(int)
    body = []
    for r in rows:
        if r.status != "ADDED":
            continue
        base = f"LF_{r.chrom}_{r.bL}"
        seen[base] += 1
        eid = base if seen[base] == 1 else f"{base}_{seen[base]}"
        body.append([eid, r.chrom, r.bL, r.bL + 1, f"LF_cluster_{r.cluster}", r.cluster, r.bL, r.bR, r.tsd,
                     _cell(r.e), _cell(r.c)])
    return _tsv(ADDED_CALL_COLUMNS, body)
