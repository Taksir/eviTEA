"""Confidence tiers (``evitea call-tiers``): one label per call, output only.

A tier says how much evidence stands behind a call and whether the call sits
where this method is known to err.  Nothing here adds, removes or moves a
call: the call files of a run are only read, and ``call_tiers.tsv`` repeats
the first six columns of the combined call table (``event_id``, ``chr``,
``start``, ``end``, ``family``, ``source``) in the same order.

Tiers (1 = most confident; the first rule that applies decides):

* library-free call (``source`` ``LIBRARY_FREE``): tier 5 ``RESIDENT`` when a
  reference TE copy of any family lies within ``NEAR_BP`` (500) bases of the
  call (its family has no name to compare), else tier 3 ``LIBRARY_FREE``;
* tier 5 ``RESIDENT``: the call interval overlaps a reference copy of the
  call's own family;
* tier 4 ``LOW_SUPPORT``: total support below ``MIN_SUPPORT`` (4) fragments;
* tier 2 ``REPEAT_RESCUE``: the repeat-aware rescue stage;
* tier 1 ``DIRECT``: every other call (exact junction on one or both sides,
  read pairs on both sides, anchor-pair extension) with support of 4 or more
  that is not inside a copy of its own family.

``STRICT_TIERS`` (tier 1) is the documented strict mode: keep only those calls.

The evidence class of a call is the stage that made it, read from the native
``evidence_type`` and reported beside the tier.  Reference copies come from
the prepared repeat annotation (BED6, FlyBase ``FBte`` ids) and the TE family
map; both interval ends are treated as inclusive when a call and a copy are
compared, and an ``FBte`` id missing from the map keeps its own name.
"""

from __future__ import annotations

import bisect
import csv
import io
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import astuple, dataclass, fields
from pathlib import Path

from .errors import InputContractError
from .shadow_layers import ShadowCallset

__all__ = [
    "EVIDENCE_CLASS",
    "MIN_SUPPORT",
    "NEAR_BP",
    "ResidentIndex",
    "STRICT_TIERS",
    "TIER_COLUMNS",
    "TIER_NAMES",
    "TierRow",
    "assign_tier",
    "call_tiers_tsv_bytes",
    "evidence_class",
    "load_resident_index",
    "tier_rows",
]

MIN_SUPPORT = 4
NEAR_BP = 500
TIER_NAMES = {1: "DIRECT", 2: "REPEAT_RESCUE", 3: "LIBRARY_FREE", 4: "LOW_SUPPORT", 5: "RESIDENT"}
STRICT_TIERS = (1,)
SOURCES = ("EVITEA", "LIBRARY_FREE")

EVIDENCE_CLASS = {
    "two_sided": "EXACT_JUNCTION_TWO_SIDED",
    "one_sided_left": "EXACT_JUNCTION_ONE_SIDED",
    "one_sided_right": "EXACT_JUNCTION_ONE_SIDED",
    "two_sided_mate_cluster_rescue": "READ_PAIRS_TWO_SIDED",
    "v3_anchor_pe_bilateral": "ANCHOR_PAIR_EXTENSION",
    "v3_anchor_pe_asymmetric_shadow": "ANCHOR_PAIR_EXTENSION",
    "v3_full_memory_anchor_pe_shadow": "ANCHOR_PAIR_EXTENSION",
    "v3_internal21_strict_a_anchor_pe_shadow": "ANCHOR_PAIR_EXTENSION",
    "msbwt_r9e_repeat_aware_crossvalidated_te_divergence": "REPEAT_AWARE_RESCUE",
}
LIBRARY_FREE_CLASS = "LIBRARY_FREE"


def evidence_class(evidence_type: str) -> str:
    """The stage class of a native ``evidence_type``; an unknown type is an error (a new stage needs a tier decision)."""

    try:
        return EVIDENCE_CLASS[evidence_type]
    except KeyError:
        raise InputContractError(f"no confidence tier is defined for evidence_type {evidence_type!r}") from None


def assign_tier(source: str, evidence_type: str | None, total_support: int | None, *, inside_own_family: bool,
                near_any: bool) -> int:
    """The tier of one call (module docstring)."""

    if source == "LIBRARY_FREE":
        return 5 if near_any else 3
    if source != "EVITEA":
        raise InputContractError(f"unknown call source {source!r}")
    klass = evidence_class(evidence_type)
    if isinstance(total_support, bool) or not isinstance(total_support, int) or total_support < 0:
        raise InputContractError(f"total support must be a non-negative integer, got {total_support!r}")
    if inside_own_family:
        return 5
    if total_support < MIN_SUPPORT:
        return 4
    return 2 if klass == "REPEAT_AWARE_RESCUE" else 1


class ResidentIndex:
    """Reference TE copies per contig: is a call inside a copy of its own family, is any copy within ``NEAR_BP``."""

    def __init__(self, records: Iterable[tuple[str, int, int, str]], near: int = NEAR_BP):
        self.near = near
        by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        for chrom, start, end, family in records:
            by_chrom[chrom].append((int(start), int(end), family))
        self.by_chrom = {chrom: sorted(v) for chrom, v in by_chrom.items()}
        self.starts = {chrom: [s for s, _, _ in v] for chrom, v in self.by_chrom.items()}
        self.max_len = {chrom: max(e - s for s, e, _ in v) for chrom, v in self.by_chrom.items()}

    def context(self, chrom: str, start: int, end: int, family: str | None) -> tuple[bool, bool]:
        """(inside a copy of ``family`` — overlap, inclusive ends —, any copy within the near distance)."""

        copies = self.by_chrom.get(chrom)
        if not copies:
            return (False, False)
        starts = self.starts[chrom]
        lo = bisect.bisect_left(starts, start - self.near - self.max_len[chrom])
        hi = bisect.bisect_right(starts, end + self.near)
        inside = near = False
        for a, b, f in copies[lo:hi]:
            if a <= end and start <= b:
                d = 0
            else:
                d = a - end if a > end else start - b
            if d <= self.near:
                near = True
                if d == 0 and family is not None and f == family:
                    inside = True
        return (inside, near)


def _family_map(path: Path) -> dict[str, str]:
    """FlyBase-style table (FBte_ID, FlyBase_Family) or two columns record<TAB>family; labels kept as written."""

    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if rows and "FBte_ID" in rows[0] and "FlyBase_Family" in rows[0]:
        i, j = rows[0].index("FBte_ID"), rows[0].index("FlyBase_Family")
        return {r[i]: r[j] for r in rows[1:] if len(r) > max(i, j)}
    return {r[0]: r[1] for r in rows if len(r) >= 2}


def load_resident_index(repeat_annotation: Path, family_map: Path) -> ResidentIndex:
    """Reference copies from the prepared BED6 repeat annotation, named through the family map."""

    names = _family_map(family_map)
    records = []
    with Path(repeat_annotation).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, record = line.rstrip("\r\n").split("\t")[:4]
            records.append((chrom, int(start), int(end), names.get(record, record)))
    return ResidentIndex(records)


@dataclass(frozen=True, slots=True)
class TierRow:
    event_id: str
    chrom: str
    start: str
    end: str
    family: str
    source: str
    evidence_class: str
    total_support: int | None
    resident_own_family: int | None
    resident_near_any: int
    tier: int
    tier_name: str
    strict: int


_COLUMN_NAMES = {"chrom": "chr"}
TIER_COLUMNS = tuple(_COLUMN_NAMES.get(f.name, f.name) for f in fields(TierRow))


def _support(value: str | None, event_id: str) -> int:
    try:
        number = int(str(value).strip())
    except ValueError:
        raise InputContractError(f"call {event_id}: total_support {value!r} is not an integer") from None
    if number < 0:
        raise InputContractError(f"call {event_id}: total_support {number} is negative")
    return number


def _added_rows(added: bytes | None) -> list[dict[str, str]]:
    if not added:
        return []
    return list(csv.DictReader(io.StringIO(added.decode("utf-8-sig"), newline=""), delimiter="\t"))


def tier_rows(native: bytes, added: bytes | None, resident: ResidentIndex) -> tuple[TierRow, ...]:
    """One row per call of the combined call table: the native calls, then the added library-free calls."""

    callset = ShadowCallset.from_bytes(native)
    missing = {"evidence_type", "total_support"} - set(callset.fields)
    if missing:
        raise InputContractError(f"native call table lacks {sorted(missing)}")
    out = []
    for row in callset.rows:
        klass = evidence_class(row["evidence_type"])
        support = _support(row["total_support"], row["event_id"])
        inside, near = resident.context(row["chr"], int(row["start"]), int(row["end"]), row["family"])
        tier = assign_tier("EVITEA", row["evidence_type"], support, inside_own_family=inside, near_any=near)
        out.append(TierRow(row["event_id"], row["chr"], row["start"], row["end"], row["family"], "EVITEA", klass,
                           support, int(inside), int(near), tier, TIER_NAMES[tier], int(tier in STRICT_TIERS)))
    for row in _added_rows(added):
        _inside, near = resident.context(row["chr"], int(row["start"]), int(row["end"]), None)
        tier = assign_tier("LIBRARY_FREE", None, None, inside_own_family=False, near_any=near)
        out.append(TierRow(row["event_id"], row["chr"], row["start"], row["end"], row["family"], "LIBRARY_FREE",
                           LIBRARY_FREE_CLASS, None, None, int(near), tier, TIER_NAMES[tier],
                           int(tier in STRICT_TIERS)))
    return tuple(out)


def call_tiers_tsv_bytes(rows: Iterable[TierRow]) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(TIER_COLUMNS)
    for row in rows:
        writer.writerow(["" if value is None else value for value in astuple(row)])
    return handle.getvalue().encode("utf-8")
