"""Output-only call annotations: a zygosity label and the TSD letters (``call_annotations.tsv``).

Both come from the breakpoint walk of the call's range (narrow calls are
walked too when the annotations are on).  Nothing here adds, removes or moves
a call; the call files are byte-identical with the annotations on or off.

* **Zygosity** uses the sides whose walk ended in a host-anchored TE junction
  (``walker.is_te``).  At that junction the walker records the share of reads
  that still follow the reference (``frac``).  The smallest share over those
  sides is reported; ``HOMOZYGOUS`` when it is at most 0.15 (the walker's own
  ``SUPP`` threshold for a supported divergence), else ``HETEROZYGOUS`` (a TE
  side has at most 0.85 by construction).  A side is only informative at a
  flank depth of at least 10 reads (``walker.D_MIN``; ``walker.D_MAX`` is
  unbounded unless ``--depth-cap`` sets a bound), which is the depth floor; with
  no TE side the label is empty.
* **TSD**: for a two-sided TE junction the target-site duplication is the
  reference between the right walk's first non-reference base and the left
  walk's (``fw[R.pos + 1 : L.pos]``, length ``L.pos - R.pos - 1``).  Empty and
  flagged ``NO_TSD`` (length 0), ``OVERLAP`` (negative), ``ONE_SIDED`` (one TE
  side: no second junction to bound the letters, length still reported when the
  walker measured it) or empty for calls without a TE junction.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Iterable
from dataclasses import astuple, dataclass, fields

from . import MIN_WIDTH
from . import walker as W

__all__ = [
    "ANNOTATION_COLUMNS",
    "AnnotationRow",
    "HOM_MAX_REF_FRACTION",
    "annotate_walk",
    "call_annotations_tsv_bytes",
    "tsd_field",
    "zygosity",
]

HOM_MAX_REF_FRACTION = W.SUPP


@dataclass(frozen=True, slots=True)
class AnnotationRow:
    event_id: str
    chrom: str
    start: int
    end: int
    walked: str
    outcome: str
    zygosity: str
    zygosity_ref_fraction: float | None
    zygosity_depth: float | None
    zygosity_sides: str
    tsd_length: int | None
    tsd: str
    tsd_flag: str


ANNOTATION_COLUMNS = tuple(f.name for f in fields(AnnotationRow))


def zygosity(L: dict, R: dict) -> tuple[str, float | None, float | None, str]:
    """(label, smallest reference-path fraction, smallest flank depth, sides used) over the TE-junction sides."""

    sides = [(name, side) for name, side in (("L", L), ("R", R)) if side.get("status") == "ok" and W.is_te(side)]
    if not sides:
        return ("", None, None, "")
    frac = min(side["frac"] for _, side in sides)
    depth = min(side["D"] for _, side in sides)
    label = "HOMOZYGOUS" if frac <= HOM_MAX_REF_FRACTION else "HETEROZYGOUS"
    return (label, frac, depth, "".join(name for name, _ in sides))


def tsd_field(outcome: str, tsd: int | None, L: dict, R: dict, fw: str) -> tuple[int | None, str, str]:
    """(length, letters, flag) of the target-site duplication."""

    if outcome == "TE_JUNCTION_TWO_SIDED":
        if tsd > 0:
            return (tsd, fw[R["pos"] + 1:L["pos"]].upper(), "OK")
        return (tsd, "", "NO_TSD" if tsd == 0 else "OVERLAP")
    if outcome == "TE_JUNCTION_ONE_SIDED_TSD":
        return (tsd, "", "ONE_SIDED")
    return (None, "", "")


def annotate_walk(walk, fw: str, *, min_width: int = MIN_WIDTH) -> AnnotationRow:
    """The annotation row of one walked call; ``fw`` is the call's contig (any case)."""

    label, frac, depth, used = zygosity(walk.L, walk.R)
    length, letters, flag = tsd_field(walk.outcome, walk.tsd, walk.L, walk.R, fw)
    walked = "WIDE" if walk.e0 - walk.s0 > min_width else "NARROW"
    return AnnotationRow(walk.event_id, walk.chrom, walk.s0, walk.e0, walked, walk.outcome, label, frac, depth, used,
                         length, letters, flag)


def call_annotations_tsv_bytes(rows: Iterable[AnnotationRow]) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(ANNOTATION_COLUMNS)
    for row in rows:
        writer.writerow(["" if value is None else value for value in astuple(row)])
    return handle.getvalue().encode("utf-8")
