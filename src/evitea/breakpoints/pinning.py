"""One-sided breakpoint pinning of calls the breakpoint switch leaves wide.

A wide call whose walk found no TE junction (decision ``KEEP_NO_TE_JUNCTION``)
gets a reported position when exactly one side's walk ended in a TE junction
label (``TE_GENOMIC`` or ``TE_CONSENSUS``, not ``UNANCHORED``) with a recorded
position; the position is that side's ``pos`` as the walker records it.  The
rule reads the side labels only.

Positions only: the call, its range in ``calls.tsv`` and ``calls.bed`` and the
native callset never change; the pin goes to ``breakpoints.tsv`` columns
``pin_side`` and ``pin_pos``.
"""

from __future__ import annotations

from collections.abc import Sequence

from ..errors import InputContractError

__all__ = ["PIN_COLUMNS", "TE_LABELS", "one_sided_position", "pin_rows"]

PIN_COLUMNS = ("pin_side", "pin_pos")
TE_LABELS = ("TE_GENOMIC", "TE_CONSENSUS")


def _is_te_label(label) -> bool:
    text = "" if label is None else str(label)
    return text.split(":")[-1] in TE_LABELS and not text.startswith("UNANCHORED")


def one_sided_position(L: dict, R: dict) -> tuple[str, int] | None:
    """('L' | 'R', position) when exactly one side carries a TE junction label with a position, else None."""

    l_te, r_te = _is_te_label(L.get("lab")), _is_te_label(R.get("lab"))
    if l_te == r_te:
        return None
    side, record = ("L", L) if l_te else ("R", R)
    pos = record.get("pos")
    return None if pos is None or pos == "" else (side, int(float(pos)))


def pin_rows(rows: Sequence, walks: Sequence) -> dict[str, tuple[str, int]]:
    """Map event_id to (side, position) for rows left wide without a TE junction; rows and walks in the same order."""

    if [row.event_id for row in rows] != [walk.event_id for walk in walks]:
        raise InputContractError("breakpoint rows and walks must name the same calls in the same order")
    out: dict[str, tuple[str, int]] = {}
    for row, walk in zip(rows, walks):
        if row.decision != "KEEP_NO_TE_JUNCTION":
            continue
        pin = one_sided_position(walk.L, walk.R)
        if pin is not None:
            out[row.event_id] = pin
    return out
