"""Walk variants of the library-free layer, each a copy of ``_p4_lib.walk`` (previous version) with one added rule.

``walk_alt``   haplotype aware walk: at a substitution with reference-letter support the context can follow the other
               letter (``follow_alt``), and ``n_both`` counts such substitutions. With ``follow_alt=False`` every key of
               ``_p4_lib.walk`` is returned unchanged, plus ``n_both``.
``walk_repeat`` walk for a walk start moved out of a repeated word: a branch that would end the walk is passed when its
               word (the 30 letters before it + the branch letter) occurs in the reference, because another reference
               copy explains it.

``left_walk_at`` / ``right_walk_at`` are ``_p4_lib.left_walk`` / ``right_walk`` with the walk start given, so one place
function serves the first pass, the moved starts and the haplotype scan. With ``strides=1`` and the plain anchor they
equal the plain functions.
"""

from __future__ import annotations

import numpy as np

from . import _p4_lib as P

__all__ = ["RefWords", "left_walk_at", "right_walk_at", "walk_alt", "walk_fn_for", "walk_repeat"]


class RefWords:
    """Sorted canonical codes of every K-letter word of the reference (both strands)."""

    def __init__(self, genome):
        parts = []
        for s in genome.values():
            c, v = P.kmer_codes(s.upper(), P.K)
            parts.append(np.unique(c[v]))
        self.codes = np.unique(np.concatenate(parts)) if parts else np.zeros(0, dtype=np.uint64)

    def has(self, w: str) -> bool:
        c, v = P.kmer_codes(w.upper(), P.K)
        if len(c) != 1 or not v[0] or len(self.codes) == 0:
            return False
        j = int(np.searchsorted(self.codes, c[0]))
        return j < len(self.codes) and self.codes[j] == c[0]


def walk_alt(C, s, i, ctx, end, D, shortcut=True, follow_alt=False):
    """``_p4_lib.walk`` that also counts ``n_both`` = substitutions (3-letter rejoin test passed) with ``nref >= thr``. With
    ``follow_alt`` the context takes the alternative letter at such a substitution. At a break with ``follow_alt``, ``frac``
    and ``T`` are site values over both the walk's context and the reference context (``frac_walk``, ``T_walk`` keep the walk's
    own)."""
    thr = P.REF_OK * D
    n_snp = n_ind = n_both = 0
    guard = 4 * (end - i) + 1000
    while i < end:
        guard -= 1
        if guard < 0:
            return dict(status="dead", i=i, n_snp=n_snp, n_ind=n_ind, why="guard", n_both=n_both)
        rb = s[i]
        if not P.is_acgt(s[i:i + 4]) or i + 4 > len(s):
            return dict(status="refN", i=i, n_snp=n_snp, n_ind=n_ind, n_both=n_both)
        nref = C.cnt(ctx + rb)
        if shortcut and nref >= thr and C.cnt(ctx) - nref < thr:
            ctx = ctx[1:] + rb
            i += 1
            continue
        n = {b: (nref if b == rb else C.cnt(ctx + b)) for b in "ACGT"}
        T = sum(n.values())
        alt = max((b for b in "ACGT" if b != rb), key=lambda b: n[b])
        if n[alt] >= thr:
            if all(C.cnt((ctx + alt + s[i + 1:i + 1 + j])[-P.K:]) >= thr for j in (1, 2, 3)):
                n_snp += 1
                if nref >= thr:
                    n_both += 1
                    ctx = ctx[1:] + (alt if follow_alt else rb)
                else:
                    ctx = ctx[1:] + alt
                i += 1
                continue
            r = P.indel_rejoin(C, s, i, ctx, alt, D)
            if r is None:
                frac = nref / T if T else 0.0
                out = dict(status="break", i=i, ctx=ctx, alt=alt, n_alt=n[alt], nref=nref, T=T, frac=frac, n_snp=n_snp,
                           n_ind=n_ind, n_both=n_both)
                if follow_alt:
                    refctx = s[i - 30:i]
                    if ctx != refctx:
                        T_site = T + C.cnt(refctx)
                        nref_site = nref + C.cnt(refctx + s[i])
                        frac_site = nref_site / T_site if T_site else 0.0
                    else:
                        T_site, frac_site = T, frac
                    out.update(frac=frac_site, T=T_site, frac_walk=frac, T_walk=T)
                return out
            n_ind += 1
            if nref >= thr:
                ctx = ctx[1:] + rb
                i += 1
            else:
                ctx, i = r[0], r[1]
            continue
        if nref >= thr:
            ctx = ctx[1:] + rb
            i += 1
            continue
        return dict(status="dead", i=i, n_snp=n_snp, n_ind=n_ind, why="no_branch", n_both=n_both)
    return dict(status="through", i=i, n_snp=n_snp, n_ind=n_ind, n_both=n_both)


def walk_fn_for(follow_alt: bool):
    return lambda C, s, i, ctx, end, D, shortcut=True: walk_alt(C, s, i, ctx, end, D, shortcut, follow_alt)


def walk_repeat(words: RefWords):
    """The walk for moved starts: ``_p4_lib.walk`` plus the rule above. ``n_par`` counts the branches passed."""

    def walk(C, s, i, ctx, end, D, shortcut=True):
        thr = P.REF_OK * D
        n_snp = n_ind = n_par = 0
        guard = 4 * (end - i) + 1000
        while i < end:
            guard -= 1
            if guard < 0:
                return dict(status="dead", i=i, n_snp=n_snp, n_ind=n_ind, n_par=n_par, why="guard")
            rb = s[i]
            if not P.is_acgt(s[i:i + 4]) or i + 4 > len(s):
                return dict(status="refN", i=i, n_snp=n_snp, n_ind=n_ind, n_par=n_par)
            nref = C.cnt(ctx + rb)
            if shortcut and nref >= thr and C.cnt(ctx) - nref < thr:
                ctx = ctx[1:] + rb
                i += 1
                continue
            n = {b: (nref if b == rb else C.cnt(ctx + b)) for b in "ACGT"}
            T = sum(n.values())
            alt = max((b for b in "ACGT" if b != rb), key=lambda b: n[b])
            if n[alt] >= thr:
                if all(C.cnt((ctx + alt + s[i + 1:i + 1 + j])[-P.K:]) >= thr for j in (1, 2, 3)):
                    n_snp += 1
                    ctx = ctx[1:] + (rb if nref >= thr else alt)
                    i += 1
                    continue
                r = P.indel_rejoin(C, s, i, ctx, alt, D)
                if r is None:
                    if words.has(ctx + alt):
                        n_par += 1
                        if nref >= thr:
                            ctx = ctx[1:] + rb
                            i += 1
                            continue
                        return dict(status="dead", i=i, n_snp=n_snp, n_ind=n_ind, n_par=n_par, why="paralog")
                    return dict(status="break", i=i, ctx=ctx, alt=alt, n_alt=n[alt], nref=nref, T=T,
                                frac=nref / T if T else 0.0, n_snp=n_snp, n_ind=n_ind, n_par=n_par)
                n_ind += 1
                if nref >= thr:
                    ctx = ctx[1:] + rb
                    i += 1
                else:
                    ctx, i = r[0], r[1]
                continue
            if nref >= thr:
                ctx = ctx[1:] + rb
                i += 1
                continue
            return dict(status="dead", i=i, n_snp=n_snp, n_ind=n_ind, n_par=n_par, why="no_branch")
        return dict(status="through", i=i, n_snp=n_snp, n_ind=n_ind, n_par=n_par)

    return walk


def left_walk_at(C, ref, cand, strides=1, walk_fn=P.walk):
    """``_p4_lib.left_walk`` with the start ``strides`` words before the run (1 = plain)."""
    a = cand["first"] - strides * P.STRIDE
    if a < 0 or not P.is_acgt(ref[a:a + P.K]):
        return dict(status="edge", i=a), a
    end = min(len(ref) - 4, cand["last"] + P.K + P.INDEL_MAX)
    return walk_fn(C, ref, a + P.K, ref[a + 1:a + P.K], end, cand["D"]), a


def right_walk_at(C, ref, cand, b=None, walk_fn=P.walk):
    """``_p4_lib.right_walk`` with the anchor word start ``b`` given (default: the plain ``last + STRIDE``)."""
    if b is None:
        b = cand["last"] + P.STRIDE
    if b + P.K > len(ref) or not P.is_acgt(ref[b:b + P.K]):
        return dict(status="edge"), b, 0, 0
    lo = max(0, cand["first"] - 2 * P.INDEL_MAX - 2 * P.K - 10)
    hi = min(len(ref), b + P.K + 2 * P.INDEL_MAX + P.K + 10)
    R = P.rc(ref[lo:hi])
    r0 = hi - b
    end_r = min(len(R) - 4, hi - (cand["first"] - P.INDEL_MAX))
    w = walk_fn(C, R, r0, R[r0 - P.K + 1:r0], end_r, cand["D"])
    if w["status"] == "break":
        w["bR"] = (hi - 1 - w["i"]) + 1
    return w, b, lo, hi
