"""Reference-guided k-mer walk and J1 junction placement over a sample msBWT.

Port of the frozen S3 research modules (``oracle_walk``, ``junction``, ``j1``,
``s3a_outcome``, ``s3_rules``; see ``tests/reference_s3``).  Every constant and
tie rule is kept; the differential tests require identical results.

All positions are 0-based.  A side is walked on the forward reference (left
side) or on its reverse complement (right side); ``site`` converts both to
forward coordinates.  The read counter and the segment aligner are injected:

* ``KmerCounter(count)``: ``count(seq)`` returns the forward-strand occurrence
  count of ``seq`` in the reads (msBWT ``countOccurrencesOfSeq``);
* ``SegmentClassifier(aligner, ...)``: ``aligner.map(seq)`` yields hits with
  minimap2/mappy attributes ``ctg, q_st, r_st, r_en, mlen, blen``.
"""

from __future__ import annotations

import bisect
import statistics
from collections.abc import Callable, Iterable

# --- reference-continuity walk (research ``oracle_walk``) -------------------
K = 31
W = 20
PRE = 20  # the walk starts W + PRE bases before the window
REF_OK, LOW, SUPP, HETMAX = 0.15, 0.65, 0.15, 0.85  # branch support (x D); reference-fraction cut-offs
from ..depth_window import D_MIN, upper_bound as _upper_bound  # noqa: E402

D_MAX = _upper_bound()  # no absolute upper bound unless --depth-cap / EVITEA_DEPTH_CAP sets one
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


class KmerCounter:
    """Both-strand read count of a k-mer, cached; a k-mer with N counts 0."""

    def __init__(self, count: Callable[[str], int]) -> None:
        self.count, self.cache, self.n = count, {}, 0

    def cnt(self, seq: str) -> int:
        value = self.cache.get(seq)
        if value is None:
            if "N" in seq:
                value = 0
            else:
                self.n += 2
                value = self.count(seq) + self.count(rc(seq))
            self.cache[seq] = value
        return value


def classify(C: KmerCounter, ref: str, a: int, e: int, w: int = W) -> dict:
    """Reference-continuity class of the window [a, e] walked left to right."""

    lo, hi = a - w - PRE, e + w + PRE
    if lo - 300 < 0 or hi + 300 > len(ref):
        return {"cls": "UNINFORMATIVE", "why": "edge"}
    fl = [C.cnt(ref[s:s + K]) for s in list(range(lo - 250, lo - 50 - K, 5)) + list(range(hi + 50, hi + 250, 5))
          if "N" not in ref[s:s + K]]
    if len(fl) < 20:
        return {"cls": "UNINFORMATIVE", "why": "flank_N"}
    D = statistics.median(fl)
    out = {"D": D}
    if not (D_MIN <= D <= D_MAX):
        out.update(cls="UNINFORMATIVE", why="depth")
        return out
    if ref[lo - K + 1:hi + K].count("N") > 0.2 * (hi - lo):
        out.update(cls="UNINFORMATIVE", why="window_N")
        return out
    path = ref[lo - K + 1:lo]  # 30-base reference anchor
    if C.cnt(ref[lo - K:lo]) < REF_OK * D:
        out.update(cls="UNINFORMATIVE", why="anchor_absent")
        return out
    minfrac, at, n_snp, n_div = 1.0, None, 0, 0
    i = lo
    while i < hi:
        ctx = path[-(K - 1):]
        n = {b: C.cnt(ctx + b) for b in "ACGT"}
        T = sum(n.values())
        rb = ref[i]
        nref = n.get(rb, 0)
        if T < REF_OK * D:  # dead end on every branch
            frac = 0.0 if nref == 0 else nref / max(T, 1)
            if frac < minfrac:
                minfrac, at = frac, i - a
            break
        frac = nref / T
        alt = max((b for b in "ACGT" if b != rb), key=lambda b: n[b])
        rejoin = False
        if n[alt] >= REF_OK * D:
            # supported alternative: substitution if it rejoins the reference within 3 bases, else divergence
            rejoin = all(C.cnt((ctx + alt + ref[i + 1:i + 1 + j])[-K:]) >= REF_OK * D for j in (1, 2, 3))
            if rejoin:
                n_snp += 1
            else:
                n_div += 1
                if frac < minfrac:
                    minfrac, at = frac, i - a
        if nref >= REF_OK * D:
            path += rb
        elif rejoin:
            path += alt
        else:
            break  # reference path lost: divergence
        i += 1
    out.update(R=round(minfrac, 4), at_offset=at, n_snp_steps=n_snp, n_div_steps=n_div, walked=i - lo, span=hi - lo)
    out["cls"] = "SUPPORTED" if minfrac <= SUPP else ("HET_LIKE" if minfrac <= HETMAX else "CONTRADICTED")
    return out


# --- divergence events and segment labels (research ``junction``) -----------
EXT = 150
MIN_SEG = 40
MIN_BLOCK = 30


def walk_events(C: KmerCounter, ref: str, a: int, e: int, w: int = W):
    """Divergence events of the walk over [a, e]; ``(None, reason)`` when the side is uninformative."""

    lo, hi = a - w - PRE, e + w + PRE
    if lo - 300 < 0 or hi + 300 > len(ref):
        return None, "edge"
    fl = [C.cnt(ref[s:s + K]) for s in list(range(lo - 250, lo - 50 - K, 5)) + list(range(hi + 50, hi + 250, 5))
          if "N" not in ref[s:s + K]]
    if len(fl) < 20:
        return None, "flank_N"
    D = statistics.median(fl)
    if not (D_MIN <= D <= D_MAX):
        return None, "depth"
    if C.cnt(ref[lo - K:lo]) < REF_OK * D:
        return None, "anchor_absent"
    path = ref[lo - K + 1:lo]
    ev = []
    i = lo
    while i < hi:
        ctx = path[-(K - 1):]
        n = {b: C.cnt(ctx + b) for b in "ACGT"}
        T = sum(n.values())
        rb = ref[i]
        nref = n.get(rb, 0)
        if T < REF_OK * D:
            break
        alt = max((b for b in "ACGT" if b != rb), key=lambda b: n[b])
        rejoin = False
        if n[alt] >= REF_OK * D:
            rejoin = all(C.cnt((ctx + alt + ref[i + 1:i + 1 + j])[-K:]) >= REF_OK * D for j in (1, 2, 3))
            if not rejoin:
                ev.append(dict(i=i, frac=nref / T, ctx=ctx, alt=alt, n_alt=n[alt], D=D))
        if nref >= REF_OK * D:
            path += rb
        elif rejoin:
            path += alt
        else:
            break
        i += 1
    return ev, D


def build_rm_index(records: Iterable[tuple[str, int, int, str]]) -> dict:
    """chrom -> (starts, ends, te_ids, max length), sorted by (start, end, id): the reference TE annotation."""

    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for chrom, start, end, te_id in records:
        by_chrom.setdefault(chrom, []).append((int(start), int(end), te_id))
    out = {}
    for chrom, values in by_chrom.items():
        ordered = sorted(values)
        out[chrom] = ([x[0] for x in ordered], [x[1] for x in ordered], [x[2] for x in ordered],
                      max(x[1] - x[0] for x in values))
    return out


class SegmentClassifier:
    """Labels a divergent segment as TE consensus, reference TE copy, or non-TE reference sequence."""

    def __init__(self, aligner, cons_names, fam_of, rm) -> None:
        self.al, self.cons, self.fam, self.rm = aligner, cons_names, fam_of, rm

    def rm_family(self, ch: str, s: int, t: int):
        if ch not in self.rm:
            return None
        st, en, ids, maxlen = self.rm[ch]
        j = bisect.bisect_right(st, t) - 1
        while j >= 0 and st[j] >= s - maxlen:
            if min(en[j], t) - max(st[j], s) >= 0.5 * (t - s):
                return ids[j]
            j -= 1
        return None


def label_hit(cl: SegmentClassifier, seg: str):
    """(label, family, matching bases, chosen hit) of a segment."""

    if len(seg) < MIN_SEG:
        return "SHORT", None, 0, None
    hits = [h for h in cl.al.map(seg) if h.blen >= MIN_BLOCK]
    if not hits:
        return "NO_HIT", None, 0, None
    best = max(h.mlen for h in hits)
    cons = [h for h in hits if h.ctg in cl.cons and h.mlen >= 0.9 * best]
    if cons:
        h = max(cons, key=lambda h: h.mlen)
        return "TE_CONSENSUS", cl.fam.get(h.ctg, h.ctg), h.mlen, h
    h = max(hits, key=lambda h: h.mlen)
    f = cl.rm_family(h.ctg, h.r_st, h.r_en)
    if f:
        return "TE_GENOMIC", cl.fam.get(f, f), h.mlen, h
    return "NON_TE_REFERENCE", f"{h.ctg}:{h.r_st}", h.mlen, h


# --- J1 junction placement (research ``j1``) --------------------------------
MATCH, MISM, GAPO, GAPE = 1, -4, 6, 1  # a gap of length L costs GAPO + GAPE * L
BAND = 12  # max |net indel| inside the host prefix
HOST_MIN = 12  # minimum extension score for a host prefix
SUFFIX_MIN = 100  # extend until this many bases follow the host prefix
EXT0, EXT_MAX = EXT, 450
NEG = -10**9
QST_MAX, ROW_SLACK = 40, 12
MAX_RESUME = 3
SYNC = 12


def _sub(x: str, y: str) -> int:
    return MATCH if (x == y and x in "ACGT") else MISM


def dp(q: str, r: str, band: int = BAND):
    """Gotoh DP, both sequences anchored at 0, cells |i-j| <= band; H[i][j] = best score of q[:i] vs r[:j]."""

    n, m = len(q), len(r)
    H = [[NEG] * (m + 1) for _ in range(n + 1)]
    E = [[NEG] * (m + 1) for _ in range(n + 1)]
    F = [[NEG] * (m + 1) for _ in range(n + 1)]
    H[0][0] = 0
    for j in range(1, min(m, band) + 1):
        E[0][j] = H[0][j] = -(GAPO + GAPE * j)
    for i in range(1, n + 1):
        if i <= band:
            F[i][0] = H[i][0] = -(GAPO + GAPE * i)
        for j in range(max(1, i - band), min(m, i + band) + 1):
            e_ = max(H[i][j - 1] - GAPO - GAPE, E[i][j - 1] - GAPE)
            f_ = max(H[i - 1][j] - GAPO - GAPE, F[i - 1][j] - GAPE)
            E[i][j], F[i][j], H[i][j] = e_, f_, max(H[i - 1][j - 1] + _sub(q[i - 1], r[j - 1]), e_, f_)
    return H


def ext_align(q: str, r: str, band: int = BAND):
    """Best prefix-vs-prefix alignment (free end): (score, qe, re); ties -> smallest qe, then re; (0, 0, 0) if none > 0."""

    H = dp(q, r, band)
    best = (0, 0, 0)
    for i in range(len(H)):
        for j in range(len(H[0])):
            if H[i][j] > best[0]:
                best = (H[i][j], i, j)
    return best


def row_best(H, Q: int):
    """(score, j): best alignment of the whole q[:Q] against some prefix r[:j]; ties -> smallest j."""

    row = H[Q]
    s = max(row)
    return s, row.index(s)


def extend(C: KmerCounter, cur: str, seq: str, n_alt: int, maxlen: int):
    """Greedy extension of the divergent branch, resumable. Returns (seq, cur, dead)."""

    stop = max(3, 0.10 * n_alt)
    while len(seq) < maxlen:
        c = cur[-(K - 1):]
        n = {b: C.cnt(c + b) for b in "ACGT"}
        b = max(n, key=n.get)
        if n[b] < stop:
            return seq, cur, True
        seq += b
        cur += b
    return seq, cur, False


def refine(C: KmerCounter, ref: str, i: int, ctx: str, alt: str, n_alt: int, cl: SegmentClassifier) -> dict:
    """J1 for one divergence at walk position i.

    The divergent branch is extended; a host prefix (banded affine alignment to
    the reference from i, score >= HOST_MIN) absorbs host variants; the rest is
    labelled.  For a TE label whose hit starts q_st <= QST_MAX bases into the
    suffix, the junction moves to the TE start when the bases before it still
    align as host (within ROW_SLACK); otherwise the junction is the end of the
    host prefix.
    """

    seq, cur, dead = extend(C, ctx + alt, alt, n_alt, EXT0)
    while True:
        r_ = ref[i:i + len(seq) + BAND]
        H = dp(seq, r_)
        sc, qe, re_ = 0, 0, 0
        for a_ in range(len(H)):
            for b_ in range(len(H[0])):
                if H[a_][b_] > sc:
                    sc, qe, re_ = H[a_][b_], a_, b_
        hs = sc
        if sc < HOST_MIN:
            qe = re_ = 0
        if dead or len(seq) - qe >= SUFFIX_MIN or len(seq) >= EXT_MAX:
            break
        seq, cur, dead = extend(C, cur, seq, n_alt, min(EXT_MAX, qe + SUFFIX_MIN))
    suf = seq[qe:]
    lab, fam, mlen, h = label_hit(cl, suf)
    qst = h.q_st if (h is not None and lab.startswith("TE_")) else 0
    pos_r, used = re_, 0
    if qst > 0 and qst <= QST_MAX:
        rs, rj = row_best(H, qe + qst)
        if rs >= max(hs if qe > 0 else 0, 0) - ROW_SLACK:
            pos_r, used = rj, 1
    if qe > 0 and len(suf) < MIN_SEG:
        lab = "HOST_VARIANT"  # divergence fully explained by a host variant
    return dict(pos=i + pos_r, pos_host=i + re_, lab=lab, fam=fam, mlen=mlen, host_q=qe, host_r=re_, host_score=hs,
                qst=qst, qst_used=used, seg=seq, dead=dead)


def walk(C: KmerCounter, ref: str, i: int, path: str, hi: int, D):
    """The divergence walk from position i with context ``path`` to hi -> (events, i_end, broke)."""

    ev = []
    while i < hi:
        ctx = path[-(K - 1):]
        n = {b: C.cnt(ctx + b) for b in "ACGT"}
        T = sum(n.values())
        rb = ref[i]
        nref = n.get(rb, 0)
        if T < REF_OK * D:
            return ev, i, False
        alt = max((b for b in "ACGT" if b != rb), key=lambda b: n[b])
        rejoin = False
        if n[alt] >= REF_OK * D:
            rejoin = all(C.cnt((ctx + alt + ref[i + 1:i + 1 + j])[-K:]) >= REF_OK * D for j in (1, 2, 3))
            if not rejoin:
                ev.append(dict(i=i, frac=nref / T, ctx=ctx, alt=alt, n_alt=n[alt], D=D))
        if nref >= REF_OK * D:
            path += rb
        elif rejoin:
            path += alt
        else:
            return ev, i, True
        i += 1
    return ev, i, False


def resume_point(ref: str, i: int, r: dict) -> int:
    """First segment offset after the host variant where segment and reference agree for SYNC bases."""

    seg, d = r["seg"], r["host_r"] - r["host_q"]
    for qv in range(1, r["host_q"] - SYNC + 1):
        if seg[qv:qv + SYNC] == ref[i + qv + d:i + qv + d + SYNC]:
            return qv
    return r["host_q"]


def choose(C: KmerCounter, ref: str, a: int, e: int, cl: SegmentClassifier):
    """Side junction with host-variant skipping -> (status, event, refine result, n host variants)."""

    lo, hi = a - W - PRE, e + W + PRE
    ev0, D = walk_events(C, ref, a, e)
    if ev0 is None:
        return D, None, None, 0
    ev, iend, broke = walk(C, ref, lo, ref[lo - K + 1:lo], hi, D)
    if [x["i"] for x in ev] != [x["i"] for x in ev0]:
        raise RuntimeError("breakpoint walk is not deterministic")
    stops = {iend} if broke else set()
    cand = sorted(ev, key=lambda x: (x["frac"], x["i"]))
    nhv = 0
    first_hv = None
    res = 0
    while cand:
        m = cand.pop(0)
        r = refine(C, ref, m["i"], m["ctx"], m["alt"], m["n_alt"], cl)
        if r["lab"] != "HOST_VARIANT":
            return "ok", m, r, nhv
        nhv += 1
        first_hv = first_hv or (m, r)
        if m["i"] in stops and res < MAX_RESUME:
            res += 1
            qv = resume_point(ref, m["i"], r)
            j = m["i"] + qv + (r["host_r"] - r["host_q"])
            if j < hi:
                new, iend, broke = walk(C, ref, j, m["ctx"] + r["seg"][:qv], hi, D)
                if broke:
                    stops.add(iend)
                cand = sorted(cand + new, key=lambda x: (x["frac"], x["i"]))
    if first_hv:
        return "ok", first_hv[0], first_hv[1], nhv
    return "no_divergence", None, None, 0


# --- site outcome (research ``s3a_outcome``) --------------------------------
TSD_MIN, TSD_MAX = -3, 30
OUTCOMES = ("TE_JUNCTION_TWO_SIDED", "TE_JUNCTION_ONE_SIDED_TSD", "REFERENCE_CONTINUOUS", "NON_TE_DISCONTINUITY",
            "UNRESOLVED")
TE_JUNCTION = frozenset(OUTCOMES[:2])
ANCHOR_SHIFTS = (0, 32, 64, 96, 128)


def side(C: KmerCounter, seq: str, a: int, e: int, cl: SegmentClassifier) -> dict:
    c = classify(C, seq, a, e)
    st, m, r, nhv = choose(C, seq, a, e, cl)
    if st != "ok":
        s = {"status": st}
    else:
        s = {"status": "ok", "frac": round(m["frac"], 3), "div_offset": m["i"] - a, "junction_offset": r["pos"] - a,
             "n_host_var": nhv}
        s.update({k: r[k] for k in ("lab", "fam", "mlen", "host_q", "host_r", "host_score", "qst", "qst_used", "seg",
                                    "dead")})
    s.update(cls=c["cls"], walked=c.get("walked"), span=c.get("span"), why=c.get("why"), D=c.get("D"))
    return s


def side_shift(C: KmerCounter, seq: str, a: int, e: int, cl: SegmentClassifier) -> dict:
    for sh in ANCHOR_SHIFTS:
        s = side(C, seq, a - sh, e, cl)
        s["shift"] = sh
        s["a_used"] = a - sh
        if s["status"] != "anchor_absent":
            return s
    return s


def rm_overlaps(rmi: dict, ch, s: int, t: int) -> bool:
    if ch not in rmi or t <= s:
        return False
    st, en, ids, maxlen = rmi[ch]
    j = bisect.bisect_left(st, t) - 1
    while j >= 0 and st[j] >= s - maxlen:
        if en[j] > s and st[j] < t:
            return True
        j -= 1
    return False


def disc(s: dict) -> bool:
    return s["status"] == "ok" and s["frac"] <= HETMAX and s.get("lab") != "HOST_VARIANT"


def is_te(s: dict) -> bool:
    return disc(s) and str(s.get("lab", "")).startswith("TE_") and s.get("host_anchored", True)


def decide(L: dict, R: dict, a: int, e: int):
    for s in (L, R):
        if s["status"] not in ("ok", "no_divergence"):
            return "UNRESOLVED", "UNINFORMATIVE_" + s["status"], None
        if s["cls"] == "UNINFORMATIVE":
            return "UNRESOLVED", "UNINFORMATIVE_" + str(s["why"]), None
    tsd = None
    if disc(L) and disc(R):
        tsd = L["pos"] - R["pos"] - 1
    geom = tsd is not None and TSD_MIN <= tsd <= TSD_MAX
    nte = is_te(L) + is_te(R)
    if nte == 2:
        return ("TE_JUNCTION_TWO_SIDED", "", tsd) if geom else ("UNRESOLVED", "TE_GEOMETRY_INCONSISTENT", tsd)
    if nte == 1:
        return ("TE_JUNCTION_ONE_SIDED_TSD", "", tsd) if geom else ("UNRESOLVED", "TE_ONE_SIDED_NO_TSD", tsd)
    if not disc(L) and not disc(R):
        full = all(s["cls"] == "CONTRADICTED" and s["walked"] == s["span"] for s in (L, R))
        return ("REFERENCE_CONTINUOUS", "", None) if full else ("UNRESOLVED", "DEAD_END_OR_PARTIAL_WALK", None)
    labs = sorted({str(s.get("lab")) for s in (L, R) if disc(s)})
    if "NON_TE_REFERENCE" in labs:
        return "NON_TE_DISCONTINUITY", "", tsd
    if any(lab.startswith("UNANCHORED") for lab in labs):
        return "UNRESOLVED", "SIDE_NOT_HOST_ANCHORED", tsd
    return "UNRESOLVED", "DISC_" + "+".join(labs), tsd


def site(C: KmerCounter, fw: str, rv: str, a: int, e: int, cl: SegmentClassifier, ch=None, rmi=None):
    """Outcome of the window [a, e]: (outcome, reason, tsd, left side, right side)."""

    Lc = len(fw)
    L = side_shift(C, fw, a, e, cl)
    R = side_shift(C, rv, Lc - 1 - e, Lc - 1 - a, cl)
    # forward coordinates: pL = first non-reference base walking right; pR = first non-reference base walking left
    if L["status"] == "ok":
        L["pos"] = L["a_used"] + L["junction_offset"]
        L["div_pos"] = L["a_used"] + L["div_offset"]
    if R["status"] == "ok":
        R["pos"] = Lc - 1 - (R["a_used"] + R["junction_offset"])
        R["div_pos"] = Lc - 1 - (R["a_used"] + R["div_offset"])
    if L["status"] == "ok":
        L["host_anchored"] = not rm_overlaps(rmi or {}, ch, L["a_used"] - W - PRE - K, L["pos"])
    if R["status"] == "ok":
        R["host_anchored"] = not rm_overlaps(rmi or {}, ch, R["pos"] + 1, Lc - (R["a_used"] - W - PRE - K))
    for s_ in (L, R):
        if s_.get("host_anchored") is False:
            s_["lab"] = "UNANCHORED:" + str(s_.get("lab"))
    o, why, tsd = decide(L, R, a, e)
    return o, why, tsd, L, R


def breakpoint(o: str, L: dict, R: dict):
    """Insertion interval [lo, hi] (forward, inclusive): two-sided -> [pR+1, pL] (TSD span); one-sided -> the TE side."""

    if o == "TE_JUNCTION_TWO_SIDED":
        x, y = R["pos"] + 1, L["pos"]
        return (min(x, y), max(x, y))
    if o == "TE_JUNCTION_ONE_SIDED_TSD":
        if is_te(L):
            return (L["pos"], L["pos"])
        return (R["pos"] + 1, R["pos"] + 1)
    return None


# --- refinement rule (research ``s3_rules``) --------------------------------
def refine_decision(s0: int, e0: int, outcome: str, bp_lo, bp_hi):
    """-> (decision, start, end). REFINED only for a TE junction whose interval [bp_lo, bp_hi] lies inside [s0, e0]."""

    if outcome not in TE_JUNCTION:
        return "KEEP_NO_TE_JUNCTION", s0, e0
    lo, hi = int(bp_lo), int(bp_hi)
    if not (s0 <= lo <= hi <= e0):
        return "KEEP_OUTSIDE_WINDOW", s0, e0
    return "REFINED", lo, hi + 1
