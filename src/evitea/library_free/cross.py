"""Rules that measure the genome and the sample instead of using values set on one species.

The merge's constants describe the genome they were set on: few repeats, few copies per family, flanks close to the reference.
Each rule here keeps the merge's test but measures what the test assumes, and acts only where the merge rules give up on a
place, so everywhere else a run writes the merge's places. All four act on the first pass (the haplotype scan keeps the merge
rules), and ``Options.merge()`` has them all off.

shared-word depth   A haplotype that differs from the reference at more than about one letter in 45 breaks most 31-letter
                    reference words, so the median word count around it reads near zero although the region is covered. Where the
                    median around a run is below the depth floor, the run is judged with the median count of the words in the same
                    window (+-DFLANK) that the sample carries (count above zero).
repeat cutoff       The repeat-content cutoff (REP_MAX: share of a place's +-REP_FLANK window made of 21-letter words with two or
                    more reference copies) was set on one reference. genome-share keeps the share of the genome the cutoff turns
                    off, not its value. The reference is cut into windows of 2 x REP_FLANK on contigs of at least 1 Mb. When a larger
                    share of its windows reach REP_MAX than in the anchor reference (``REPEAT_ANCHOR``), the cutoff is the repeat
                    content of the window at the anchor's share, counted from the most repetitive. It is never below REP_MAX, and on
                    the anchor reference it is REP_MAX.
poly-A tails        A segment that ends in a run of one letter stops the extension (poly-A retrotransposons). The first pass's
                    extension crosses such a run (``extension.Extender(cross=True)``), and its family rule counts a poly-A tail at a
                    junction toward the two TSDs (``priority.tail_mark``). Both act only on segments that the merge's extension
                    stops short and on families that the merge's family rule leaves out.
recurring sequence  A local rearrangement copies its own flank once. A sequence inserted at several places with target-site
                    duplications is a mobile element even when an old copy of it sits within 2 kb, which is common in genomes dense
                    in short repeats. A place whose segments the filter rejected as a local or segmental copy is kept when the place
                    itself has a TSD (or a poly-A tail with poly-A tails on) and a rejected segment groups with a member of a family
                    that passes the TSD rule. The one-site check then judges it like every other library-free call.

Each place a rule changed carries the rule's name in ``cross_rule`` (places table of the first pass).
"""

from __future__ import annotations

import numpy as np

from . import _p4_lib as P
from .place import _firewall, _record, process_place

__all__ = ["MIN_PROFILE_CONTIG", "REPEAT_ANCHOR", "CrossContext", "cross_place", "readmit_recurring", "repeat_cutoff",
           "repeat_profile", "shared_depth"]

MIN_PROFILE_CONTIG = 1_000_000
# Windows of the anchor reference (the D. melanogaster reference REP_MAX was set on: EviTEA's prepared dm6, reference.fa sha256
# eeb896ecb227ec5b2c779faa46f064cf781adddaae509c7f50b9d6503c194397, contigs of at least 1 Mb: chr2L, chr2R, chr3L, chr3R, chr4,
# chrX, chrY) and how many of them reach REP_MAX, measured with ``evitea repeat-profile``.
REPEAT_ANCHOR: dict | None = {"windows": 137_139, "at_cutoff": 21_303}
COPY_REJECTS = ("local", "sd")
SHORT_REJECTS = ("short", "sat")


# ---------------------------------------------------------------- repeat cutoff
def repeat_profile(genome, rm, min_len=None):
    """Repeat content of every window of 2 x REP_FLANK (centred at REP_FLANK, 3 x REP_FLANK, ...) of the contigs of at least
    ``min_len`` letters, in FASTA order; the value ``RefMulti.rep_frac`` gives at the window's centre. Windows with no valid word
    are left out."""
    k, f = P.REP_K, P.REP_FLANK
    min_len = MIN_PROFILE_CONTIG if min_len is None else min_len
    out = []
    for name, seq in genome.items():
        if len(seq) < min_len:
            continue
        codes, valid = P.kmer_codes(seq.upper(), k)
        multi = rm._isin(codes, k) & valid
        cv = np.concatenate([[0], np.cumsum(valid, dtype=np.int64)])
        cm = np.concatenate([[0], np.cumsum(multi, dtype=np.int64)])
        del codes, valid, multi
        xs = np.arange(f, len(seq) - f + 1, 2 * f, dtype=np.int64)
        lo, hi = xs - f, xs + f - k + 1
        nv = cv[hi] - cv[lo]
        nm = cm[hi] - cm[lo]
        keep = nv > 0
        out.append(nm[keep] / nv[keep])
    return np.concatenate(out) if out else np.zeros(0)


def repeat_cutoff(reps, anchor=None):
    """(cutoff, summary) from the windows' repeat content. ``anchor`` = dict(windows=, at_cutoff=) of the anchor reference
    (default ``REPEAT_ANCHOR``)."""
    anchor = REPEAT_ANCHOR if anchor is None else anchor
    if anchor is None:
        raise ValueError("--repeat-cutoff genome-share needs the anchor reference's window counts (cross.REPEAT_ANCHOR); "
                         "measure them with `evitea repeat-profile` or use --repeat-cutoff fixed")
    n = int(len(reps))
    n_at = int(np.count_nonzero(reps >= P.REP_MAX))
    a_n, a_at = int(anchor["windows"]), int(anchor["at_cutoff"])
    summary = dict(windows=n, at_rep_max=n_at, share_at_rep_max=round(n_at / n, 6) if n else None,
                   anchor_share=round(a_at / a_n, 6), rep_max=P.REP_MAX)
    if n == 0 or n_at * a_n <= a_at * n:
        cut = P.REP_MAX
    else:
        m = max(1, (a_at * n) // a_n)            # windows the cutoff turns off at the anchor's share (at least the top one)
        srt = np.sort(reps)[::-1]
        cut = max(P.REP_MAX, float(srt[m - 1]))
    summary["cutoff"] = cut
    summary["windows_at_cutoff"] = int(np.count_nonzero(reps >= cut))
    return cut, summary


# ---------------------------------------------------------------- shared-word depth
def shared_depth(pos, c, cand):
    """Median count of the sampled words within DFLANK of the run that the sample carries (count above zero), or None."""
    lo = int(np.searchsorted(pos, cand["first"] - P.DFLANK, side="left"))
    hi = int(np.searchsorted(pos, cand["last"] + P.DFLANK, side="right"))
    w = c[lo:hi]
    w = w[w > 0]
    return float(np.median(w)) if len(w) else None


# ---------------------------------------------------------------- first-pass places under the cross-species rules
class CrossContext:
    """What the rules need for one contig: the place context with the cross-species extender and cutoff, the word counts and
    the switches."""

    def __init__(self, pc, pos, c, opts, d_min=None):
        self.pc, self.pos, self.c = pc, pos, c
        self.shared_depth = opts.shared_depth
        self.polya = opts.polya_tails
        self.recurring = opts.recurring_sequence
        self.walk_out = opts.walk_out
        self.d_lo = P.D_MIN if d_min is None else d_min
        self.stats = dict(shared_depth=0, shared_depth_tried=0, repeat_cutoff=0, repeat_cutoff_tried=0, polya_tails=0,
                          polya_tails_tried=0)


def _changed(a, b):
    return any(a.get(k) != b.get(k) for k in ("stage", "why", "bL", "bR", "lseg", "rseg", "seg_keep", "seg_reject", "D"))


def _rejected_sides(pl, kinds):
    out = set()
    for r in str(pl.get("seg_reject") or "").split(","):
        if ":" in r and r.split(":")[1] in kinds:
            out.add(r.split(":")[0])
    return out


def _reextend(X, cand, pl, lw, rw):
    """The merge place with the sides its filter rejected as short or satellite extended again by the cross-species extender;
    None when no segment changes."""
    pc = X.pc
    sides = _rejected_sides(pl, SHORT_REJECTS)
    if not sides:
        return None
    new = dict(pl)
    changed = False
    if "L" in sides and pl.get("lseg"):
        s, stop = pc.extend(cand, "L", pc.C, lw["ctx"], lw["alt"], lw["n_alt"], fork=P.FORK)
        if s != pl["lseg"]:
            new["lseg"], new["lstop"], changed = s, stop, True
    if "R" in sides and rw is not None and rw.get("status") == "break" and pl.get("rseg"):
        rs, stop = pc.extend(cand, "R", pc.C, rw["ctx"], rw["alt"], rw["n_alt"], fork=P.FORK)
        if P.rc(rs) != pl["rseg"]:
            new["rseg"], new["rstop"], changed = P.rc(rs), stop, True
    if not changed:
        return None
    keep, rej = _firewall(pc, (("L", new["lseg"]), ("R", new["rseg"])), new["bL"])
    new["seg_keep"], new["seg_reject"] = ",".join(keep), ",".join(rej)
    new["stage"] = "S4_pass" if keep else "S4_reject"
    for k in ("chkL_T", "chkL_junction", "chkR_T", "chkR_junction"):
        new.pop(k, None)
    if pc.record_check and new["stage"] == "S4_pass":
        _record(new, lw, rw, new["rseg"])
    return new


def cross_place(X, cand, pl, lw, rw):
    """The first-pass place under the cross-species rules, given the merge place ``pl`` and its walks. Returns a new dict (the
    merge place is not changed)."""
    out, olw, orw = None, lw, rw
    D = cand["D"]
    if X.shared_depth and D < X.d_lo:
        X.stats["shared_depth_tried"] += 1
        d2 = shared_depth(X.pos, X.c, cand)
        if d2 is not None and d2 >= X.d_lo:
            p2, olw, orw = process_place(X.pc, dict(cand, D=d2), walkout=X.walk_out)
            if _changed(p2, pl):
                out = p2
                out["cross_rule"] = "shared_depth"
                X.stats["shared_depth"] += 1
            else:
                olw, orw = lw, rw
    elif pl.get("stage") == "S4_abstain" and X.pc.rep_max > P.REP_MAX and "depth" not in str(pl.get("why") or ""):
        why = str(pl.get("why") or "").split(",")
        rf = pl.get("rep") if "rep" in why else pl.get("_walkout_rep")
        # ``rep`` is rounded to 3 decimals in the place; the walk itself decides with the exact value
        if rf is not None and ("rep" in why or "_walkout_rep" in pl) and rf <= X.pc.rep_max + 0.001:
            X.stats["repeat_cutoff_tried"] += 1
            p2, olw, orw = process_place(X.pc, cand, walkout=X.walk_out)
            if _changed(p2, pl):
                out = p2
                out["cross_rule"] = "repeat_cutoff"
                X.stats["repeat_cutoff"] += 1
            else:
                olw, orw = lw, rw
    elif X.polya and pl.get("stage") in ("S4_pass", "S4_reject"):
        X.stats["polya_tails_tried"] += bool(_rejected_sides(pl, SHORT_REJECTS))
        p2 = _reextend(X, cand, pl, lw, rw)
        if p2 is not None:
            out = p2
            out["cross_rule"] = "polya_tails"
            X.stats["polya_tails"] += 1
    if out is None:
        out = dict(pl)
    if X.recurring and X.pc.record_check and out.get("stage") == "S4_reject" and _rejected_sides(out, COPY_REJECTS) \
            and "chkL_T" not in out and "chkR_T" not in out and olw is not None and olw.get("status") == "break":
        _record(out, olw, orw, out.get("rseg", ""))
    return out


# ---------------------------------------------------------------- recurring sequence (grouping)
def _has_site_mark(p, polya):
    from .priority import tail_mark

    t = p.get("tsd")
    return bool(p.get("two_sided")) and ((t is not None and P.TSD_LO <= t <= P.TSD_HI) or (polya and tail_mark(p)))


def readmit_recurring(places, gp, cl, fam, *, polya: bool, confirm_fn=P.confirm):
    """Rejected places (not duplicates, no mitochondrial side) with a TSD at the place and a side rejected as a local or segmental
    copy whose segment groups with a member of a family that passes the TSD rule (shared words, then ``confirm_fn``, as
    ``_p4_lib.group`` joins two places). Each is changed to S4_pass with those sides kept; its other sides stay rejected.
    Returns the changed places."""
    members = []
    for p, c in zip(gp, cl):
        if fam[c]["vtsd"]:
            for s in p["segs"]:
                m = P.mask_lc(s)
                members.append((m, P.kset(m)))
    if not members:
        return []
    index: dict = {}
    for mi, (_, ks) in enumerate(members):
        for w in ks:
            index.setdefault(w, []).append(mi)
    changed = []
    for p in places:
        if p.get("stage") != "S4_reject" or p.get("dup") or not _has_site_mark(p, polya):
            continue
        if "chkL_T" not in p and "chkR_T" not in p:
            continue
        rej = [r for r in str(p.get("seg_reject") or "").split(",") if ":" in r]
        if any(r.split(":")[1] == "numt" for r in rej):
            continue
        keep = []
        for side in ("L", "R"):
            seg = p.get("lseg" if side == "L" else "rseg") or ""
            if len(seg) < P.SEG_MIN or not ({f"{side}:{k}" for k in COPY_REJECTS} & set(rej)):
                continue
            m = P.mask_lc(seg)
            ks = P.kset(m)
            shared: dict = {}
            for w in ks:
                for mi in index.get(w, ()):
                    shared[mi] = shared.get(mi, 0) + 1
            tries = sorted(((n / min(len(ks), len(members[mi][1])), mi) for mi, n in shared.items()
                            if min(len(ks), len(members[mi][1])) and n >= P.G_SHARE * min(len(ks), len(members[mi][1]))),
                           key=lambda t: (-t[0], t[1]))
            if any(confirm_fn(m, members[mi][0]) for _, mi in tries):
                keep.append(side)
        if keep:
            p["stage"] = "S4_pass"
            p["seg_keep"] = ",".join(keep)
            p["seg_reject"] = ",".join(r for r in rej if r.split(":")[0] not in keep)
            p["cross_rule"] = "recurring_sequence"
            changed.append(p)
    return changed
