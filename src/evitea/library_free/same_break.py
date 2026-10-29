"""One place per break: places that each won their duplicate choice and share a left break are one insertion.

With first-pass priority, ``priority.dedupe`` chains a walked-out place on the break its plain walk made, not on the break of
the walk-out. A walk that starts further out can end on the break of a place in another chain. Both places then win their
duplicate choice and both are grouped: the insertion counts twice in its family, gives two calls at one break, and the two
places carry one name in the segments table (``chrom:bL:side:cl<cluster>``).

Of each set of chosen places with one break, the place that the duplicate choice ranks first stays (``priority.dedupe`` on the
set, chained on the shared break). The other places of the set chose for their duplicate groups an insertion that is the kept
place's, so every place of those groups is set aside before the grouping, with stage ``S4_same_break`` and ``dup`` 1. The
other duplicate groups do not change, because a whole group leaves the chain. Without a walked-out place every chain
coordinate is the break, and nothing is set aside.

The step reads only the place records, so saved scan pieces stay valid.
"""

from __future__ import annotations

from collections import defaultdict

from . import priority as PR
from .options import Options

__all__ = ["STAGE", "set_aside_same_break"]

STAGE = "S4_same_break"


def set_aside_same_break(places: list[dict], opts: Options) -> int:
    """Set aside the duplicate groups whose chosen place shares its break with a chosen place ranked before it; return how
    many places."""
    if not (opts.walk_out and opts.first_pass_priority) or not any(p.get("walked") for p in places):
        return 0
    PR.dedupe(places, priority=True, yield_order=opts.het_yield)
    members = defaultdict(list)
    chosen = defaultdict(list)
    for p in places:
        if p.get("stage") in PR.RANK and p.get("bL") is not None:
            members[id(p["_rep"])].append(p)
            if not p["dup"]:
                chosen[(p["chrom"], p["bL"])].append(p)
    n = 0
    for same in chosen.values():
        if len(same) < 2:
            continue
        one = [dict(p, bL0=p["bL"]) for p in same]
        PR.dedupe(one, priority=True, yield_order=opts.het_yield)
        kept = same[next(i for i, q in enumerate(one) if not q["dup"])]
        for p, q in zip(same, one):
            if q["dup"]:
                for r in members[id(p)]:
                    r.update(stage=STAGE, dup=1, _rep=kept)
                    n += 1
    return n
