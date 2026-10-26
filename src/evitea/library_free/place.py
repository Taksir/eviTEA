"""One candidate run -> one place (walk and decision), for the first pass and the haplotype scan.

``process_place`` is ``_p4_lib.process`` with these optional parts, each off by default:

* the walk start of the left walk ``left_strides`` words before the run, and the right anchor after the left break
  (``right_after_bL``, ``safe_right``): the haplotype scan (``scan.py``);
* the walk functions (``walk_fn``): the haplotype-aware walk (``walks.walk_alt``);
* the extension (``extender``): anchored extension, poly-A crossing (``extension.Extender``);
* the walk-out (``walkout``): when the only reasons to abstain are that the word where the left or right walk starts is
  repeated in the reference or has more than ``ANCHOR_RC`` x D reads, the start moves outward one ``STRIDE`` at a time, at most
  ``REP_FLANK`` (500 bp, the flank over which the layer judges repeat content), to the first word that passes the same test,
  and the walk starts there. A moved walk uses ``walks.walk_repeat``: a branch that would end the walk is passed when its word
  occurs in the reference, because another reference copy explains it. The repeat test is then applied at the new break. If no
  start passes, the new walk does not break or the repeat test fails, the plain abstention record is kept unchanged;
* the record of the one-site check (``record_check``): at each side that exists, the walk's read total T at the break and the
  junction (30 host letters + 30 inserted letters), read here because the reads are at hand.

With every option off the place equals ``_p4_lib.process``'s (tests/test_place.py).
"""

from __future__ import annotations

from . import _p4_lib as P
from .walks import left_walk_at, right_walk_at, walk_repeat

__all__ = ["PlaceContext", "WALK_MAX", "process_place", "right_only_place"]

WALK_MAX = P.REP_FLANK // P.STRIDE
INS_USE = 30


class PlaceContext:
    """What one contig's places need: the counter, the contig, the repeat words, the firewall, the extension."""

    def __init__(self, C, ref, rm, fw, chrom, extender=None, words=None, record_check=False, d_min=None):
        self.C, self.ref, self.rm, self.fw, self.chrom = C, ref, rm, fw, chrom
        self.d_min = d_min          # lower bound of the local-depth window; None = the fixed ``_p4_lib.D_MIN``
        self.extender = extender or P.extend
        self.words = words
        self.record_check = record_check
        self.walk_log: list[dict] = []
        self.stats = dict(walkout_tried=0, walkout_moved=0, walkout_no_start=0, walkout_no_break=0, walkout_rep=0)
        self._paralog = walk_repeat(words) if words is not None else None

    def extend(self, cand, side, *args, **kw):
        where = getattr(self.extender, "where", None)
        if where is not None:
            where.clear()
            where.update(chrom=self.chrom, first=cand["first"], last=cand["last"], side=side)
        return self.extender(*args, **kw)


def _anchor_ok(pc, p, D):
    """The plain walk-start test of ``_p4_lib.process`` (True = the start word may seed a walk)."""
    if p < 0:
        return False
    w = pc.ref[p:p + P.K]
    return not (len(w) < P.K or not P.is_acgt(w) or pc.rm.word_multi(w) or pc.C.cnt(w) > P.ANCHOR_RC * D)


def _first_start(pc, D, p0, step):
    for j in range(1, WALK_MAX + 1):
        if _anchor_ok(pc, p0 + step * j * P.STRIDE, D):
            return j
    return None


def _firewall(pc, sides, bL):
    keep, rej = [], []
    for side, sg in sides:
        if len(sg) < P.SEG_MIN:
            if sg:
                rej.append(side + ":short")
            continue
        r = None
        if pc.fw.local(sg, pc.ref, bL):
            r = "local"
        elif pc.fw.numt(sg):
            r = "numt"
        elif pc.fw.sd(sg):
            r = "sd"
        elif P.satellite(sg):
            r = "sat"
        if r:
            rej.append(side + ":" + r)
        else:
            keep.append(side)
    return keep, rej


def _record(pl, lw, rw, rseg):
    """The one-site check's inputs at each side that exists."""
    if pl.get("lseg"):
        pl["chkL_T"] = lw["T"]
        pl["chkL_junction"] = lw["ctx"] + pl["lseg"][:INS_USE]
    if pl.get("two_sided"):
        pl["chkR_T"] = rw["T"]
        pl["chkR_junction"] = rw["ctx"] + P.rc(rseg)[:INS_USE]


def process_place(pc, cand, *, left_strides=1, right_after_bL=False, safe_right=False, walk_fn=P.walk, walkout=False,
                  fork=P.FORK):
    """Walk and decide one candidate. Return (place, left walk dict, right walk dict or None); the walk dicts are for the caller
    only. ``stage`` is where the place stopped: S2_through | S2_dead | S2_refN | S2_edge | S4_abstain | S4_reject | S4_pass."""
    C, ref, rm = pc.C, pc.ref, pc.rm
    if walkout and pc._paralog is None:
        raise ValueError("the walk-out needs the reference words (PlaceContext(words=RefWords(genome)))")
    D = cand["D"]
    pl = dict(chrom=pc.chrom, first=cand["first"], last=cand["last"], n_low=cand["n_low"], D=round(D, 2))
    lw, a = left_walk_at(C, ref, cand, left_strides, walk_fn)
    if lw["status"] != "break":
        pl["stage"] = "S2_" + lw["status"]
        return pl, lw, None
    pl["bL"] = lw["i"]
    why = []
    if not ((P.D_MIN if pc.d_min is None else pc.d_min) <= D <= P.D_MAX):
        why.append("depth")
    b = cand["last"] + P.STRIDE
    if right_after_bL:
        b1 = max(b, lw["i"] - P.PAIR_LO)
        if not safe_right:
            b = b1
        else:  # the moved anchor only if its word is clean and has a depth-like count; else the plain anchor
            w1 = ref[b1:b1 + P.K]
            if len(w1) == P.K and P.is_acgt(w1):
                c1 = C.cnt(w1)
                if P.LOWF * D <= c1 <= 2 * D:
                    b = b1
    for nm, p in (("anchorL", a), ("anchorR", b)):
        w = ref[p:p + P.K]
        if len(w) < P.K or not P.is_acgt(w) or rm.word_multi(w) or C.cnt(w) > P.ANCHOR_RC * D:
            why.append(nm)
    rf = rm.rep_frac(ref, lw["i"])
    pl["rep"] = None if rf is None else round(rf, 3)
    if rf is not None and rf >= P.REP_MAX:
        why.append("rep")
    jr, wrec, bL0 = 0, None, lw["i"]
    if why and walkout and set(why) <= {"anchorL", "anchorR"}:
        if left_strides != 1 or right_after_bL:
            raise ValueError("the walk-out belongs to the first pass (plain walk starts)")
        moved = _walk_out(pc, cand, D, a, b, why, lw)
        if moved is None:
            pl.update(stage="S4_abstain", why=",".join(why))
            return pl, lw, None
        lw, jr, rf2, wrec = moved
        pl["bL"] = lw["i"]
        pl["rep"] = None if rf2 is None else round(rf2, 3)
        pl["walked"] = 1
        pl["bL0"] = bL0
        pl["_abstain0"] = dict(why=",".join(why), bL=bL0, rep=None if rf is None else round(rf, 3))
    elif why:
        pl.update(stage="S4_abstain", why=",".join(why))
        return pl, lw, None
    # two-sided walk and extraction
    lseg, lstop = pc.extend(cand, "L", C, lw["ctx"], lw["alt"], lw["n_alt"], fork=fork)
    if jr:
        rw, _, _, _ = right_walk_at(C, ref, cand, b + jr * P.STRIDE, pc._paralog)
    else:
        rw, _, _, _ = right_walk_at(C, ref, cand, b, walk_fn)
    if wrec is not None:
        wrec.update(right_paralog=rw.get("n_par", 0), rstat=rw["status"])
    pl.update(fracL=round(lw["frac"], 3), lseg=lseg, lstop=lstop, rstat=rw["status"])
    rseg = ""
    if rw["status"] == "break":
        rs, rstop = pc.extend(cand, "R", C, rw["ctx"], rw["alt"], rw["n_alt"], fork=fork)
        rseg = P.rc(rs)
        pl.update(bR=rw["bR"], fracR=round(rw["frac"], 3), rstop=rstop)
    pl["rseg"] = rseg
    two = rw["status"] == "break" and P.PAIR_LO <= pl["bL"] - rw["bR"] <= P.PAIR_HI
    pl["two_sided"] = int(two)
    pl["delta"] = pl["bL"] - rw["bR"] if rw["status"] == "break" else None
    pl["tsd"] = max(0, pl["delta"]) if two else None
    pl["geno"] = P.genotype(min(lw["frac"], rw["frac"]) if rw["status"] == "break" else lw["frac"])
    pl["pos"] = pl["bL"]
    keep, rej = _firewall(pc, (("L", lseg), ("R", rseg)), pl["bL"])
    pl["seg_keep"] = ",".join(keep)
    pl["seg_reject"] = ",".join(rej)
    pl["stage"] = "S4_pass" if keep else "S4_reject"
    if pc.record_check and pl["stage"] == "S4_pass":
        _record(pl, lw, rw, rseg)
    return pl, lw, rw


def _walk_out(pc, cand, D, a, b, why, lw):
    """(left walk, right strides moved, repeat fraction, record) after moving the failing walk starts, or None."""
    st = pc.stats
    st["walkout_tried"] += 1
    jl = _first_start(pc, D, a, -1) if "anchorL" in why else 0
    jr = _first_start(pc, D, b, +1) if "anchorR" in why else 0
    if jl is None or jr is None:
        st["walkout_no_start"] += 1
        return None
    if jl:
        lw, _ = left_walk_at(pc.C, pc.ref, cand, 1 + jl, pc._paralog)
        if lw["status"] != "break":
            st["walkout_no_break"] += 1
            return None
    rf = pc.rm.rep_frac(pc.ref, lw["i"])
    if rf is not None and rf >= P.REP_MAX:
        st["walkout_rep"] += 1
        return None
    st["walkout_moved"] += 1
    rec = dict(chrom=pc.chrom, first=cand["first"], last=cand["last"], why=",".join(why), left_strides=jl, right_strides=jr,
               bL=lw["i"], left_paralog=lw.get("n_par", 0) if jl else 0)
    pc.walk_log.append(rec)
    return lw, jr, rf, rec


def right_only_place(pc, cand, walk_fn_plain, walk_fn_follow, fork=P.FORK):
    """Right-only fallback of the haplotype scan: the right walk alone from the plain anchor ``b = last + STRIDE`` (the plain walk;
    rerun with the alternative-following walk if it ends through or dead after counting n_both >= 1). A right break forms a
    one-sided place with bL = bR = the break. Returns (place or None, row or None, ran)."""
    C, ref, rm = pc.C, pc.ref, pc.rm
    b = cand["last"] + P.STRIDE
    rw, _, _, _ = right_walk_at(C, ref, cand, b, walk_fn_plain)
    follow = 0
    if rw.get("status") in ("through", "dead") and rw.get("n_both", 0) >= 1:
        follow = 1
        rw, _, _, _ = right_walk_at(C, ref, cand, b, walk_fn_follow)
    if rw.get("status") != "break":
        return None, None, 1 if follow else 0
    D = cand["D"]
    bR = rw["bR"]
    pl = dict(chrom=pc.chrom, first=cand["first"], last=cand["last"], n_low=cand["n_low"], D=round(D, 2))
    pl["bL"] = bR
    why = []
    if not ((P.D_MIN if pc.d_min is None else pc.d_min) <= D <= P.D_MAX):
        why.append("depth")
    w = ref[b:b + P.K]
    if len(w) < P.K or not P.is_acgt(w) or rm.word_multi(w) or C.cnt(w) > P.ANCHOR_RC * D:
        why.append("anchorR")
    rf = rm.rep_frac(ref, bR)
    pl["rep"] = None if rf is None else round(rf, 3)
    if rf is not None and rf >= P.REP_MAX:
        why.append("rep")
    if why:
        pl.update(stage="S4_abstain", why=",".join(why))
        return pl, None, 1
    rs, rstop = pc.extend(cand, "R", C, rw["ctx"], rw["alt"], rw["n_alt"], fork=fork)
    rseg = P.rc(rs)
    pl.update(fracL=round(rw["frac"], 3), lseg="", lstop="", rstat="break", bR=bR, fracR=round(rw["frac"], 3), rstop=rstop,
              rseg=rseg)
    pl["two_sided"] = 0
    pl["delta"] = None
    pl["tsd"] = None
    pl["geno"] = P.genotype(rw["frac"])
    pl["pos"] = bR
    keep, rej = _firewall(pc, (("R", rseg),), bR)
    pl["seg_keep"] = ",".join(keep)
    pl["seg_reject"] = ",".join(rej)
    pl["stage"] = "S4_pass" if keep else "S4_reject"
    if pc.record_check and pl["stage"] == "S4_pass":
        pl["chkR_T"] = rw["T"]
        pl["chkR_junction"] = rw["ctx"] + rs[:INS_USE]
    row = (pc.chrom, cand["first"], cand["last"], b, bR, pl["stage"], follow, rw.get("n_both"), rw.get("frac_walk", ""),
           rw.get("frac"), rw.get("T_walk", ""), rw.get("T"))
    return pl, row, 1
