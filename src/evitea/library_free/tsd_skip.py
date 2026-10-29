"""TSD skip of the family-identity filter: a diverged family that still carries insertion evidence is kept.

The family-identity filter (``--rule-c``) drops a family whose members share little inserted DNA (median identity at or below
RULE_C_MAX): old, diverged resident copies. Some real families are that diverged too (a mixed family of mammalian SINEs, MITEs
with RIP mutations). A resident copy has no target-site duplication of its own that the sample's reads can show; a family of
insertions has one at most places. The skip asks the family's trusted members:

* trusted members are the family's calls that a library call confirms (NEAR_EVITEA) or that pass the one-site check at the
  reference-path fraction of an insertion (not above REFERENCE_PATH_CEILING). A call that the filter would have dropped is judged
  only on these, since the others are mostly check failures of old copies.
* the family is kept (not dropped by the filter) when it has at least ``minimum`` trusted members and at least ``share`` of them
  have a TSD of 2 to 20 bp.

A kept family goes on to the other statuses (NEAR_EVITEA, DROPPED_CHECK, ...) as if the filter were off. The first pass and the
haplotype scan are judged the same way, each on its own families; the scan keeps its trusted re-grouping as well.

The bar (``--rule-c-tsd-share``, default 0.3) sits in the gap measured on the training runs: families that the filter should drop
carry TSDs at 0.11 to 0.16 of their trusted members, real families at 0.39 or more. It was chosen after those runs; held-out
runs test it. ``--rule-c-tsd-skip off`` gives the filter alone. The skip reads the family calls only (outside
``runtime.CODE_FILES``), so saved scan pieces are reused with it on or off. It needs the one-site check
(``--check-scope`` other than none): the trusted members come from it.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

__all__ = ["Skip", "has_tsd", "judge", "lines", "summary"]

TSD_MIN, TSD_MAX = 2, 20     # target-site duplication of a family member (bp)
MIN_TRUSTED = 5              # trusted members a family needs to be judged
SHARE = 0.3                  # share of the trusted members with a TSD at which the family is kept


@dataclass(frozen=True)
class Skip:
    """``--rule-c-tsd-skip`` / ``--rule-c-tsd-share`` / ``--rule-c-tsd-min``."""

    on: bool = True
    share: float = SHARE
    minimum: int = MIN_TRUSTED

    def __post_init__(self):
        if not (0.0 < self.share <= 1.0):
            raise ValueError(f"--rule-c-tsd-share must be above 0 and at most 1, got {self.share!r}")
        if not isinstance(self.minimum, int) or isinstance(self.minimum, bool) or self.minimum < 1:
            raise ValueError(f"--rule-c-tsd-min must be an integer of 1 or more, got {self.minimum!r}")

    def as_dict(self) -> dict:
        return dict(on=self.on, share=self.share, minimum=self.minimum, tsd=[TSD_MIN, TSD_MAX])


def has_tsd(tsd: str | None) -> bool:
    """Is the TSD length (text of the places table) 2 to 20 bp?"""

    try:
        return TSD_MIN <= int(float(tsd)) <= TSD_MAX
    except (TypeError, ValueError):
        return False


def judge(scan: str, families: Mapping[str, tuple[float, int, Iterable[str]]], skip: Skip) -> tuple[frozenset, list[dict]]:
    """(clusters kept, report). ``families``: cluster -> (family similarity c, calls in the family, TSD texts of its trusted members)
    for the families the filter would drop. The report has one entry per family, most trusted members first."""

    kept, report = set(), []
    for cluster, (c, calls, tsds) in families.items():
        tsds = list(tsds)
        trusted = len(tsds)
        with_tsd = sum(has_tsd(t) for t in tsds)
        share = with_tsd / trusted if trusted else None
        keep = trusted >= skip.minimum and share >= skip.share
        if keep:
            kept.add(cluster)
        report.append(dict(scan=scan, cluster=cluster, c=None if c != c else round(c, 4), calls=calls, trusted=trusted, tsd=with_tsd,
                           tsd_share=None if share is None else round(share, 4), kept=keep))
    report.sort(key=lambda r: (-r["trusted"], r["cluster"]))
    return frozenset(kept), report


def summary(skip: Skip, reports: Iterable[list[dict]], *, note: str = "") -> dict:
    """The summary block: settings, counts, and every family the filter would drop with its TSD share (families without a trusted
    member have no share; they are counted, not listed)."""

    families = [r for rep in reports for r in rep]
    listed = [r for r in families if r["trusted"]]
    out = dict(skip.as_dict(), judged=len(families), kept=sum(r["kept"] for r in families),
               without_trusted_member=len(families) - len(listed), families=listed)
    if note:
        out["note"] = note
    return out


def lines(report: Iterable[dict], minimum: int) -> list[str]:
    """Progress lines: the families that can decide (at least ``minimum`` trusted members)."""

    return [f"library-free: TSD skip: {r['scan']} family {r['cluster']}: {r['trusted']} trusted of {r['calls']} calls, "
            f"TSD share {r['tsd_share']:.2f}, {'kept' if r['kept'] else 'dropped by the family-identity filter'}"
            for r in report if r["trusted"] >= minimum]
