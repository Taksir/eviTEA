"""Scans of one contig: the first pass and the haplotype scan (run detection, walk and decision).

First pass   ``_p4_lib`` run detection (31-mer counts at stride 10, local depth, runs of low words) and ``place.process_place`` for every
             run, with the options that belong to it: anchored extension, walk-out, poly-A crossing. With all of them off it
             is ``detector.scan_contig``.
Haplotype scan  (``het_scan`` safe-start | second-walk | full). In a sample with two haplotypes a walk that keeps the reference letter cannot see
             a junction on the other haplotype, and a junction that only half the reads show leaves fewer low words. The scan
             adds a fraction test to the run detection: a 41-mer whose count is below ``TAU`` x the sample's median 41-mer/31-mer ratio x its
             31-mer count flags the word (forward) or the word before it (backward), as a run of low 31-mer counts does.
             Per run:
               safe-start   the left walk starts up to four words before the run (the furthest start whose words exist, sit exactly
                   STRIDE apart, are valid and are not flagged), the right anchor after the left break if its word is clean and
                   has a depth-like count;
               second-walk  safe-start plus a second chance: a place that ends with no left break (through or dead) after a walk that passed at
                   least one substitution with reference-letter support is walked again with the context following the other
                   letter at such substitutions;
               full         second-walk plus an earlier start, 30 letters further left, following the other letter, for a run still without a
                   left break, and the right walk alone, as a one-sided place with bL = bR, for a run still without a left
                   break.
             ``Place.origin`` tells which of these made the place ('' for the plain walk).
"""

from __future__ import annotations

import numpy as np

from . import _p4_lib as P
from .extension import Extender
from .options import Options
from .place import PlaceContext, process_place, right_only_place
from .walks import walk_fn_for

__all__ = ["TAU", "X_LEN", "candidates_from_low", "count_x", "depth_floor", "low_masks", "pick_strides", "scan_contig_passes"]

TAU = 1 - P.REF_OK
X_LEN = P.K + 10  # 41-mer

ORIGINS = ("", "second_chance", "early_start", "right_only")
CACHE_MAX = 2_000_000


def candidates_from_low(pos, low, dep):
    """Return runs of flagged words <= RUN_GAP apart, as ``_p4_lib.candidates`` does after its mask. Keys: first, last, n_low, D, j0, j1."""
    idx = np.flatnonzero(low)
    out = []
    if len(idx) == 0:
        return out
    g0 = 0
    for t in range(1, len(idx) + 1):
        if t == len(idx) or pos[idx[t]] - pos[idx[t - 1]] > P.RUN_GAP:
            g = idx[g0:t]
            if len(g) >= P.MIN_LOW:
                out.append(dict(first=int(pos[g[0]]), last=int(pos[g[-1]]), n_low=int(len(g)),
                                D=float(np.median(dep[g])), j0=int(g[0]), j1=int(g[-1])))
            g0 = t
    return out


def low_masks(c, nf, dep, x):
    """Return the masks of flagged words from c (31-mer counts, -1 invalid), nf (N fraction), dep (local depth) and x (41-mer counts at the same starts, -1 invalid)."""
    with np.errstate(invalid="ignore"):
        ok = (c >= 0) & (nf <= P.NFRAC) & ~np.isnan(dep)
        old = ok & (c < P.LOWF * dep)
        sel = ok & (c >= P.D_MIN) & (x >= 0)
        r0 = float(np.median(x[sel] / c[sel])) if sel.any() else float("nan")
        xm = np.concatenate([[-1], x[:-1]]) if len(x) else x
        fwd = ok & (c > 0) & (x >= 0) & (x < TAU * r0 * c)
        bwd = ok & (c > 0) & (xm >= 0) & (xm < TAU * r0 * c)
        bwd[:1] = False
    return dict(old=old, fwd=fwd, bwd=bwd, low=old | fwd | bwd, r0=r0, ok=ok)


def count_x(C, ref, pos):
    """Return x, where x[j] is the count of the 41-mer at pos[j], or -1 if it does not fit or has a non-ACGT letter."""
    x = np.full(len(pos), -1, dtype=np.int64)
    for j in range(len(pos)):
        p = int(pos[j])
        if p + X_LEN <= len(ref):
            w = ref[p:p + X_LEN]
            if P.is_acgt(w):
                x[j] = C.cnt(w)
    return x


def pick_strides(j0, pos, ok, low):
    """Return the furthest safe left start (4, 3 or 2 strides before the run), else the plain 1."""
    for s in (4, 3, 2):
        if j0 - s < 0:
            continue
        w = slice(j0 - s, j0)
        if (np.diff(pos[j0 - s:j0 + 1]) == P.STRIDE).all() and ok[w].all() and not low[w].any():
            return s
    return 1


def depth_floor(depth, opts: Options):
    """Return the lower bound of the local-depth window of a contig: None (the fixed ``_p4_lib.D_MIN``) or, with the relative floor,
    the smaller of ``D_MIN`` and half the median local depth of the contig's valid words."""
    if not opts.relative_depth_floor:
        return None
    valid = depth[~np.isnan(depth)]
    if not len(valid):
        return None
    return min(float(P.D_MIN), 0.5 * float(np.median(valid)))


def _first_pass(counter, ref, chrom, rm, fw, opts: Options, words, cands, info, record_check, d_min=None):
    ext = Extender(anchored=opts.anchored_extension, cross=opts.polya_cross)
    pc = PlaceContext(counter, ref, rm, fw, chrom, ext, words if opts.walk_out else None, record_check, d_min)
    places = []
    for cand in cands:
        pl, _, _ = process_place(pc, cand, walkout=opts.walk_out)
        places.append(pl)
        if len(getattr(counter, "cache", ())) > CACHE_MAX:
            counter.clear()
    info["first_pass"] = dict(candidates=len(cands), walk=dict(pc.stats), extension=dict(ext.stats))
    return places, pc.walk_log, ext.cross_log


def _het_pass(counter, ref, chrom, rm, fw, opts: Options, cands, pos, m, info, record_check, d_min=None):
    mode = opts.het_scan
    anchored = opts.anchored_extension and opts.second_scan_anchored
    ext = Extender(anchored=anchored, cross=anchored and opts.polya_cross)
    pc = PlaceContext(counter, ref, rm, fw, chrom, ext, None, record_check, d_min)
    plain, follow = walk_fn_for(False), walk_fn_for(True)
    n = dict(second_chance_runs=0, second_chance_places=0, early_start_runs=0, early_start_places=0, right_only_runs=0,
             right_only_places=0)
    places = []
    for cand in cands:
        st = pick_strides(cand["j0"], pos, m["ok"], m["low"])
        if mode == "safe-start":
            pl, _, _ = process_place(pc, cand, left_strides=st, right_after_bL=True, safe_right=True)
            places.append(pl)
        else:
            pl, lw, rw = process_place(pc, cand, left_strides=st, right_after_bL=True, safe_right=True, walk_fn=plain)
            pl["origin"] = ""
            if opts.second_chance and pl["stage"] in ("S2_through", "S2_dead") and lw.get("n_both", 0) >= 1:
                n["second_chance_runs"] += 1
                pl2, lw2, rw2 = process_place(pc, cand, left_strides=st, right_after_bL=True, safe_right=True, walk_fn=follow)
                if lw2["status"] == "break":
                    pl, lw, rw = pl2, lw2, rw2
                    pl["origin"] = "second_chance"
                    n["second_chance_places"] += 1
            if mode == "full":
                if opts.early_start and pl["stage"] in ("S2_through", "S2_dead"):
                    a2 = cand["first"] - st * P.STRIDE - 30
                    if a2 >= 0 and P.is_acgt(ref[a2:a2 + P.K]):
                        n["early_start_runs"] += 1
                        pl3, lw3, rw3 = process_place(pc, cand, left_strides=st + 3, right_after_bL=True, safe_right=True,
                                                      walk_fn=follow)
                        if lw3["status"] == "break":
                            pl, lw, rw = pl3, lw3, rw3
                            pl["origin"] = "early_start"
                            n["early_start_places"] += 1
                if opts.right_fallback and pl["stage"].startswith("S2_"):
                    pl4, _row, _ran = right_only_place(pc, cand, plain, follow)
                    n["right_only_runs"] += 1
                    if pl4 is not None:
                        pl = pl4
                        pl["origin"] = "right_only"
                        n["right_only_places"] += 1
            places.append(pl)
        if len(getattr(counter, "cache", ())) > CACHE_MAX:
            counter.clear()
    info["het_pass"] = dict(candidates=len(cands), mode=mode, extension=dict(ext.stats), **n)
    return places


def scan_contig_passes(counter, sequence, chrom, rm, firewall, opts: Options, words=None):
    """Return (first-pass places, haplotype-scan places or None, info) of one contig. The 31-mer counts are read once; the
    haplotype scan reads the 41-mer counts besides. ``counter`` is a ``_p4_lib.Counter`` built with ``cache=False`` whose cache
    is switched on after the run detection, as in ``detector.scan_contig``."""
    ref = sequence.upper()
    info: dict = {}
    pos, c = P.scan_counts(counter, ref)
    depth, nfrac = P.local_depth(c)
    cands = P.candidates(pos, c, depth, nfrac)
    het = opts.het_scan != "off"
    if het:
        x = count_x(counter, ref, pos)
        m = low_masks(c, nfrac, depth, x)
        hcands = candidates_from_low(pos, m["low"], depth)
        info.update(r0=m["r0"], n_old=int(m["old"].sum()), n_fwd=int(m["fwd"].sum()), n_bwd=int(m["bwd"].sum()),
                    n_low=int(m["low"].sum()))
    counter.use = True
    record = opts.check_scope != "none"
    d_min = depth_floor(depth, opts)
    if d_min is not None:
        info["depth_floor"] = round(d_min, 2)
    first, walk_log, cross_log = _first_pass(counter, ref, chrom, rm, firewall, opts, words, cands, info, record, d_min)
    second = None
    if het:
        second = _het_pass(counter, ref, chrom, rm, firewall, opts, hcands, pos, m, info, record, d_min)
    return first, second, dict(info, walk_log=walk_log, cross_log=cross_log)
