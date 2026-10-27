"""Grouping of places into insertions and families (clusters, family rule, family-identity filter), with the first-pass priority and the order of the second-chance parts.

Principle: a second-chance source (the anchored extension, the walk-out, the haplotype scan's second walk, early start and right
side alone) may add a place where the first pass has none, and never changes a first-pass duplicate choice, grouping or family
decision. Each source relaxes a test that protects precision, so on its own it can override a right decision with a wrong one.

dedupe        ``_p4_lib.dedupe`` (one place per insertion: places whose breaks chain within MERGE_D) with
              * first-pass priority: a place is ranked first by the stage it has when its anchored segments are ignored (a pass
                that rests on anchored segments only is a reject), and a walked-out place by the stage it had before the walk
                (S4_abstain), chained on the break its plain walk made;
              * the order of origins (``yield_order``): an early-start or right-side-alone place yields to the plain or
                second-chance walk of its chain; between those two the stage and the two-sided bonus decide.
group         ``_p4_lib.group`` (single linkage on shared words and a banded alignment) except that a pair holding a walked-out
              place may not join two clusters that both hold a place that is not walked out (first-pass priority).
cluster_members  the member segments the family-identity filter compares: only plain-extension segments when at least two members hold one, and
              none from a walked-out place (first-pass priority).
families      ``_p4_lib.families``; with ``polya_tsd`` a member whose junction holds a poly-A tail counts toward the two TSDs.
"""

from __future__ import annotations

from collections import defaultdict

from . import _p4_lib as P
from .extension import is_anchored

__all__ = ["cluster_members", "dedupe", "families", "stage_without_anchored", "group", "group_and_call", "revert_walk_out"]

RANK = {"S4_pass": 1, "S4_reject": 2, "S4_abstain": 3}
ORIGIN_RANK = {"": 0, None: 0, "second_chance": 0, "early_start": 1, "right_only": 1}
SEG_MIN = P.SEG_MIN


def _bonus(stage, p):
    return RANK[stage] - (0.5 if stage == "S4_pass" and p.get("two_sided") else 0)


def stage_without_anchored(p):
    """S4_pass only if a kept side has a plain-extension segment; a pass that rests on anchored sides only is S4_reject."""
    if p.get("stage") != "S4_pass":
        return p.get("stage")
    stop = {"L": p.get("lstop"), "R": p.get("rstop")}
    kept = str(p.get("seg_keep") or "").split(",")
    return "S4_pass" if any(s and not is_anchored(stop[s]) for s in kept) else "S4_reject"


def revert_walk_out(places):
    """Give every walked-out place back its plain abstention record (the detection kept it in ``_abstain0``)."""
    for p in places:
        if p.get("walked"):
            ab = p["_abstain0"]
            base = {k: p[k] for k in ("chrom", "first", "last", "n_low", "D")}
            p.clear()
            p.update(base, bL=ab["bL"], rep=ab["rep"], stage="S4_abstain", why=ab["why"])
    return places


def dedupe(places, *, priority: bool, yield_order: bool = False):
    """One place per insertion; sets ``dup`` (0/1) and ``_rep`` on every place with a left break."""
    chain = (lambda p: p["bL0"] if p.get("walked") else p["bL"]) if priority else (lambda p: p["bL"])
    idx = [i for i, p in enumerate(places) if p.get("stage") in RANK and p.get("bL") is not None]
    idx.sort(key=lambda i: (places[i]["chrom"], chain(places[i]), i))
    groups, cur = [], []
    for i in idx:
        if cur and (places[i]["chrom"] != places[cur[-1]]["chrom"] or chain(places[i]) - chain(places[cur[-1]]) > P.MERGE_D):
            groups.append(cur)
            cur = []
        cur.append(i)
    if cur:
        groups.append(cur)

    def key(i):
        p = places[i]
        fz = "S4_abstain" if (priority and p.get("walked")) else (stage_without_anchored(p) if priority else p["stage"])
        if yield_order:
            return (RANK[fz], ORIGIN_RANK[p.get("origin")], _bonus(fz, p), _bonus(p["stage"], p), chain(p), i)
        if priority:
            return (_bonus(fz, p), _bonus(p["stage"], p), chain(p), i)
        return (_bonus(p["stage"], p), p["bL"], i)

    for g in groups:
        best = min(g, key=key)
        for i in g:
            places[i]["dup"] = int(i != best)
            places[i]["_rep"] = places[best]
    return places


def group(places, *, priority: bool, confirm_fn=P.confirm):
    """Single linkage; with ``priority`` a pair holding a walked-out place may not join two clusters that both hold a place that
    is not walked out. Without any walked-out place it is ``_p4_lib.group``."""
    if not (priority and any(p.get("walked") for p in places)):
        return P.group(places, confirm_fn)
    segs = []
    for pi, p in enumerate(places):
        for s in p["segs"]:
            m = P.mask_lc(s)
            segs.append((pi, m, P.kset(m)))
    idx = {}
    for si, (_, _, ks) in enumerate(segs):
        for w in ks:
            idx.setdefault(w, []).append(si)
    shared = {}
    for si, (pi, _, ks) in enumerate(segs):
        for w in ks:
            for sj in idx[w]:
                if sj > si and segs[sj][0] != pi:
                    shared[(si, sj)] = shared.get((si, sj), 0) + 1
    pairs = []
    for (si, sj), n in shared.items():
        mn = min(len(segs[si][2]), len(segs[sj][2]))
        if mn and n >= P.G_SHARE * mn:
            pairs.append((-n / mn, segs[si][0], segs[sj][0], si, sj))
    pairs.sort()
    uf = P.UF(len(places))
    plain = [not p.get("walked") for p in places]
    has_plain = list(plain)
    nal = 0
    for _, pa, pb, si, sj in pairs:
        ra, rb = uf.f(pa), uf.f(pb)
        if ra == rb:
            continue
        if has_plain[ra] and has_plain[rb] and not (plain[pa] and plain[pb]):
            continue
        nal += 1
        if confirm_fn(segs[si][1], segs[sj][1]):
            uf.u(pa, pb)
            has_plain[uf.f(pa)] = has_plain[ra] or has_plain[rb]
    return [uf.f(i) for i in range(len(places))], nal


def tail_mark(p):
    """A poly-A tail at a junction: a run of LC_RUN T in the first K letters of the left segment, or of LC_RUN A in the last K
    letters of the right segment (segments are in forward orientation; the right one ends at the junction)."""
    return ("T" * P.LC_RUN) in (p.get("lseg") or "")[:P.K] or ("A" * P.LC_RUN) in (p.get("rseg") or "")[-P.K:]


def families(places, cl, *, polya_tsd: bool):
    if not polya_tsd:
        return P.families(places, cl)
    out = {}
    for p, c in zip(places, cl):
        f = out.setdefault(c, dict(n=0, two=0, tsd=0))
        f["n"] += 1
        f["two"] += p.get("two_sided", 0)
        t = p.get("tsd")
        f["tsd"] += int((t is not None and P.TSD_LO <= t <= P.TSD_HI) or tail_mark(p))
    for f in out.values():
        f["pair"] = int(f["n"] >= 2 and f["two"] >= 2)
        f["vtsd"] = int(f["pair"] and f["tsd"] >= 2)
    return out


def _group_once(places, opts, polya_tsd):
    dedupe(places, priority=opts.first_pass_priority, yield_order=opts.het_yield)
    gp = [p for p in places if p["stage"] == "S4_pass" and not p.get("dup")]
    for p in gp:
        keep = p["seg_keep"].split(",")
        p["segs"] = [s for s, side in ((p.get("lseg", ""), "L"), (p.get("rseg", ""), "R")) if side in keep]
    cl, _ = group(gp, priority=opts.first_pass_priority)
    return gp, cl, families(gp, cl, polya_tsd=polya_tsd)


def group_and_call(places, opts, *, first: bool = True):
    """Group every S4_pass place (``detector.group_and_call`` with the rules above) and add cluster, track, fam_pair, fam_vtsd.
    ``first``: the places are the first pass's, where the cross-species rules of ``opts`` apply (poly-A tails in the family rule,
    the recurring sequence: ``cross.py``); the haplotype scan's are grouped with the merge rules."""
    if not opts.walk_out:
        revert_walk_out(places)
    polya_tsd = opts.polya_tsd or (first and opts.polya_tails)
    gp, cl, fam = _group_once(places, opts, polya_tsd)
    if first and opts.recurring_sequence:
        from .cross import readmit_recurring

        if readmit_recurring(places, gp, cl, fam, polya=opts.polya_tails):
            gp, cl, fam = _group_once(places, opts, polya_tsd)
    for p, c in zip(gp, cl):
        p["cluster"] = c
        p["track"] = P.classify_place(p, fam[c])
        p["fam_pair"] = fam[c]["pair"]
        p["fam_vtsd"] = fam[c]["vtsd"]
    return gp


def cluster_members(places, segments, *, priority: bool):
    """cluster -> member segment lists (segments of >= 60 letters) of the grouped places, in table order. With ``priority`` a
    cluster with at least two members holding a plain-extension segment keeps only those, and a walked-out place gives none.
    ``places`` are rows of text fields."""
    allm = defaultdict(list)
    fv = defaultdict(list)
    for p in places:
        if not (p["stage"] == "S4_pass" and p["dup"] != "1" and p["cluster"] != ""):
            continue
        d = segments.get((p["chrom"], p["bL"], p["cluster"]), {})
        if len(d.get("lseg", "")) != int(p["lseg_len"]) or len(d.get("rseg", "")) != int(p["rseg_len"]):
            raise ValueError(f"segments do not match the places table at {p['chrom']}:{p['bL']}")
        allm[p["cluster"]].append([d[s] for s in ("lseg", "rseg") if len(d.get(s, "")) >= SEG_MIN])
        if priority:
            if p.get("walked") == "1":
                fv[p["cluster"]].append([])
                continue
            stop = {"lseg": p.get("lstop"), "rseg": p.get("rstop")}
            fv[p["cluster"]].append([d[s] for s in ("lseg", "rseg") if len(d.get(s, "")) >= SEG_MIN and not is_anchored(stop[s])])
    if not priority:
        return dict(allm)
    return {c: (fv[c] if sum(1 for m in fv[c] if m) >= 2 else ms) for c, ms in allm.items()}
