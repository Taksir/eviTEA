"""Apply breakpoint refinement to the wide calls of a native callset."""

from __future__ import annotations

import csv
import io
from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from dataclasses import astuple, dataclass, fields

from ..callset import _families
from ..errors import InputContractError
from ..families import canonical_family
from ..shadow_layers import ShadowCallset
from . import MIN_WIDTH
from . import walker as W

__all__ = [
    "BREAKPOINT_COLUMNS",
    "MIN_WIDTH",
    "BreakpointRow",
    "RefinementResult",
    "SiteWalk",
    "breakpoint_row",
    "breakpoints_tsv_bytes",
    "refine_native",
    "refined_coordinates",
    "refinement_summary",
    "walk_native",
]


@dataclass(frozen=True, slots=True)
class BreakpointRow:
    """One wide call: its native range [s0, e0), the walk outcome, and the reported range [s1, e1)."""

    event_id: str
    evidence_type: str
    chrom: str
    s0: int
    e0: int
    outcome: str
    reason: str
    tsd: int | None
    bp_lo: int | None
    bp_hi: int | None
    decision: str
    s1: int
    e1: int
    call_families: str
    site_family: str
    family_agree: int | None


BREAKPOINT_COLUMNS = tuple(f.name for f in fields(BreakpointRow))


@dataclass(frozen=True, slots=True)
class SiteWalk:
    """The walk of one native call's range [s0, e0): site outcome and both side records (``walker.site``)."""

    event_id: str
    evidence_type: str
    chrom: str
    s0: int
    e0: int
    outcome: str
    reason: str
    tsd: int | None
    L: dict
    R: dict
    call_families: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RefinementResult:
    """Breakpoint rows plus the optional extras: one-sided pins and per-call annotations (``None`` = switch off)."""

    rows: tuple[BreakpointRow, ...]
    pins: dict[str, tuple[str, int]] | None = None
    annotations: tuple | None = None

    @classmethod
    def of(cls, value) -> "RefinementResult":
        return value if isinstance(value, cls) else cls(tuple(value))


def _site_families(outcome: str, L: dict, R: dict) -> list[str]:
    if outcome not in W.TE_JUNCTION:
        return []
    return sorted({s.get("fam") for s in (L, R) if W.is_te(s) and s.get("fam")})


def walk_native(
    native: bytes,
    *,
    reference: Mapping[str, str],
    counter: W.KmerCounter,
    classifier: W.SegmentClassifier,
    rm_index: Mapping,
    min_width: int = MIN_WIDTH,
    include_narrow: bool = False,
    progress: Callable[[str], None] | None = None,
) -> tuple[SiteWalk, ...]:
    """Walk every native call with end - start > min_width (every call with ``include_narrow``); native row order."""

    rows = ShadowCallset.from_bytes(native).rows
    out: list[SiteWalk] = []
    fw_cache: dict[str, tuple[str, str]] = {}
    todo = [row for row in rows if include_narrow or int(row["end"]) - int(row["start"]) > min_width]
    for n, row in enumerate(todo, 1):
        chrom, s0, e0 = row["chr"], int(row["start"]), int(row["end"])
        if chrom not in reference:
            raise InputContractError(f"call {row['event_id']!r} lies on {chrom!r}, absent from the reference")
        if chrom not in fw_cache:
            fw_cache.clear()
            fw = reference[chrom].upper()
            fw_cache[chrom] = (fw, W.rc(fw))
        fw, rv = fw_cache[chrom]
        outcome, reason, tsd, L, R = W.site(counter, fw, rv, s0, e0, classifier, chrom, rm_index)
        out.append(SiteWalk(row["event_id"], row.get("evidence_type", ""), chrom, s0, e0, outcome, reason, tsd, L, R,
                            tuple(_families(row))))
        if len(counter.cache) > 2_000_000:
            counter.cache.clear()
        if progress is not None and (n % 50 == 0 or n == len(todo)):
            progress(f"breakpoint refinement: {n}/{len(todo)} calls walked")
    return tuple(out)


def breakpoint_row(walk: SiteWalk) -> BreakpointRow:
    """The breakpoint-table row of one walked call (refinement rule of ``walker.refine_decision``)."""

    interval = W.breakpoint(walk.outcome, walk.L, walk.R)
    lo, hi = interval if interval else (None, None)
    decision, s1, e1 = W.refine_decision(walk.s0, walk.e0, walk.outcome, lo, hi)
    site_fams = _site_families(walk.outcome, walk.L, walk.R)
    agree = (int(bool({canonical_family(f) for f in site_fams} & {canonical_family(f) for f in walk.call_families}))
             if site_fams else None)
    return BreakpointRow(walk.event_id, walk.evidence_type, walk.chrom, walk.s0, walk.e0, walk.outcome, walk.reason,
                         walk.tsd, lo, hi, decision, s1, e1, ";".join(walk.call_families), ";".join(site_fams), agree)


def refine_native(
    native: bytes,
    *,
    reference: Mapping[str, str],
    counter: W.KmerCounter,
    classifier: W.SegmentClassifier,
    rm_index: Mapping,
    min_width: int = MIN_WIDTH,
    progress: Callable[[str], None] | None = None,
) -> tuple[BreakpointRow, ...]:
    """Walk every native call with end - start > min_width; native row order; no call is added or removed."""

    walks = walk_native(native, reference=reference, counter=counter, classifier=classifier, rm_index=rm_index,
                        min_width=min_width, progress=progress)
    return tuple(breakpoint_row(walk) for walk in walks)


def refined_coordinates(rows: Iterable[BreakpointRow]) -> dict[str, tuple[int, int]]:
    """Map event_id to the refined [start, end) for REFINED rows only."""

    return {row.event_id: (row.s1, row.e1) for row in rows if row.decision == "REFINED"}


def breakpoints_tsv_bytes(
    rows: Iterable[BreakpointRow], pins: Mapping[str, tuple[str, int]] | None = None
) -> bytes:
    """The breakpoint table; with ``pins`` (pinning on) two columns, pin_side and pin_pos, are appended."""

    from .pinning import PIN_COLUMNS

    handle = io.StringIO(newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(BREAKPOINT_COLUMNS if pins is None else (*BREAKPOINT_COLUMNS, *PIN_COLUMNS))
    for row in rows:
        values = ["" if value is None else value for value in astuple(row)]
        if pins is not None:
            values += list(pins.get(row.event_id, ("", "")))
        writer.writerow(values)
    return handle.getvalue().encode("utf-8")


def refinement_summary(rows: Iterable[BreakpointRow]) -> dict:
    rows = tuple(rows)
    return {
        "min_width": MIN_WIDTH,
        "wide_calls": len(rows),
        "decisions": dict(sorted(Counter(row.decision for row in rows).items())),
        "outcomes": dict(sorted(Counter(row.outcome for row in rows).items())),
    }
