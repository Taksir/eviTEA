"""Canonical outputs of the default engine, derived from the native callset.

The scientific result of an EviTEA v1 run is the **native callset**: the
frozen Route-B v3 callset schema (15 columns), written byte for byte as
``calls.native.tsv``.  It is the file whose SHA-256 reproduces the frozen
development callsets.  ``calls.tsv``, ``calls.bed`` and ``unresolved.tsv`` are
deterministic *views* of it and of the integration layers' fail-closed
events; they add no scientific decision.

Native coordinates are kept as written by the producers:

* ``start == end``: a breakpoint (zero-based interbase point); split-junction
  and all v3 extension events.
* ``start < end``: a breakpoint range (mate-rescue and repeat-aware rescue
  events), reported unchanged.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import io
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

from .contracts import stable_json_dumps
from .errors import InputContractError
from .shadow_layers import ShadowCallset, UnresolvedShadowEvent

__all__ = [
    "CALL_COLUMNS",
    "CALL_ORIGINS",
    "BED_COLUMNS",
    "UNRESOLVED_COLUMNS",
    "CallsetState",
    "calls_bed_bytes",
    "calls_tsv_bytes",
    "unresolved_tsv_bytes",
]

# Native ``evidence_type`` -> public call origin.  Unknown evidence types fail
# closed: every native row must come from a frozen producer.
CALL_ORIGINS: Mapping[str, str] = MappingProxyType(
    {
        "two_sided": "INCUMBENT_SPLIT_JUNCTION",
        "one_sided_left": "INCUMBENT_SPLIT_JUNCTION",
        "one_sided_right": "INCUMBENT_SPLIT_JUNCTION",
        "two_sided_mate_cluster_rescue": "INCUMBENT_MATE_RESCUE",
        "msbwt_r9e_repeat_aware_crossvalidated_te_divergence": "INCUMBENT_REPEAT_AWARE",
        "v3_anchor_pe_bilateral": "EXTENSION_ANCHOR_PE_BILATERAL",
        "v3_anchor_pe_asymmetric_shadow": "EXTENSION_ANCHOR_PE_ASYMMETRIC",
        "v3_full_memory_anchor_pe_shadow": "EXTENSION_FULL_MEMORY",
        "v3_internal21_strict_a_anchor_pe_shadow": "EXTENSION_INTERNAL21",
    }
)

CALL_COLUMNS = (
    "sample_id",
    "call_id",
    "chrom",
    "family",
    "families",
    "coordinate_kind",
    "start",
    "end",
    "left_coord",
    "right_coord",
    "call_origin",
    "evidence_type",
    "left_support",
    "right_support",
    "total_support",
    "support_pair_ids",
)
BED_COLUMNS = (
    "chrom",
    "start0",
    "end0",
    "call_id",
    "score",
    "strand",
    "family",
    "coordinate_kind",
    "call_origin",
)
UNRESOLVED_COLUMNS = (
    "sample_id",
    "unresolved_id",
    "stage",
    "reason",
    "chrom",
    "lo",
    "hi",
    "family",
    "participant_ids",
)
POINT = "INTERBASE_POINT"
RANGE = "RANGE"


@dataclass(frozen=True, slots=True)
class CallsetState:
    """What a call engine hands to output finalization."""

    sample_id: str
    native: bytes
    unresolved: tuple[UnresolvedShadowEvent, ...] = ()
    layer_rows: Mapping[str, int] = field(default_factory=dict)
    layer_sha256: Mapping[str, str] = field(default_factory=dict)
    counts: Mapping[str, Any] = field(default_factory=dict)
    profile: str = "production"

    def __post_init__(self) -> None:
        if not self.sample_id:
            raise InputContractError("sample_id must be non-empty text")
        ShadowCallset.from_bytes(self.native)  # schema check

    @property
    def callset(self) -> ShadowCallset:
        return ShadowCallset.from_bytes(self.native)

    @property
    def native_sha256(self) -> str:
        return hashlib.sha256(self.native).hexdigest()

    def summary(self) -> dict[str, Any]:
        rows = self.callset.rows
        return {
            "sample_id": self.sample_id,
            "profile": self.profile,
            "final_call_count": len(rows),
            "native_callset_sha256": self.native_sha256,
            "calls_by_origin": dict(sorted(Counter(_origin(row) for row in rows).items())),
            "layer_rows": dict(self.layer_rows),
            "layer_sha256": dict(self.layer_sha256),
            "unresolved_count": len(self.unresolved),
            "unresolved_reasons": dict(sorted(Counter(row.reason for row in self.unresolved).items())),
            "engine_counts": dict(self.counts),
        }


def _origin(row: Mapping[str, str]) -> str:
    evidence = row.get("evidence_type", "")
    try:
        return CALL_ORIGINS[evidence]
    except KeyError as exc:
        raise InputContractError(
            f"native call {row.get('event_id')!r} has unknown evidence_type {evidence!r}"
        ) from exc


def _int(row: Mapping[str, str], key: str) -> int:
    try:
        return int(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise InputContractError(
            f"native call {row.get('event_id')!r} has non-integer {key}"
        ) from exc


def _families(row: Mapping[str, str]) -> list[str]:
    text = row.get("families", "")
    try:
        value = ast.literal_eval(text) if text else [row.get("family", "")]
    except (SyntaxError, ValueError) as exc:
        raise InputContractError(f"native call {row.get('event_id')!r} has malformed families") from exc
    if isinstance(value, str):
        value = [value]
    return sorted({str(item) for item in value})


def _pair_ids(row: Mapping[str, str]) -> list[int]:
    text = row.get("fragment_ids", "") or ""
    try:
        return sorted({int(item) for item in text.split(",") if item.strip()})
    except ValueError as exc:
        raise InputContractError(f"native call {row.get('event_id')!r} has malformed fragment_ids") from exc


def _opt_int(row: Mapping[str, str], key: str) -> int | None:
    text = row.get(key, "")
    if text in ("", None):
        return None
    value = _int(row, key)
    return None if value < 0 else value


def _call_view(sample_id: str, row: Mapping[str, str]) -> dict[str, Any]:
    start, end = _int(row, "start"), _int(row, "end")
    if end < start:
        raise InputContractError(f"native call {row.get('event_id')!r} has end < start")
    return {
        "sample_id": sample_id,
        "call_id": row["event_id"],
        "chrom": row["chr"],
        "family": row["family"],
        "families": _families(row),
        "coordinate_kind": POINT if start == end else RANGE,
        "start": start,
        "end": end,
        "left_coord": _opt_int(row, "left_coord"),
        "right_coord": _opt_int(row, "right_coord"),
        "call_origin": _origin(row),
        "evidence_type": row.get("evidence_type", ""),
        "left_support": _opt_int(row, "left_support"),
        "right_support": _opt_int(row, "right_support"),
        "total_support": _opt_int(row, "total_support"),
        "support_pair_ids": _pair_ids(row),
    }


def _text(value: Any) -> str:
    if value is None:
        return "."
    if isinstance(value, (list, tuple)):
        return stable_json_dumps(list(value))
    text = str(value)
    if any(character in text for character in "\t\r\n"):
        raise InputContractError("TSV scalars cannot contain tabs or newlines")
    return text


def _tsv(columns: tuple[str, ...], rows: Iterable[Mapping[str, Any]]) -> bytes:
    handle = io.StringIO(newline="")
    handle.write("\t".join(columns) + "\n")
    for row in rows:
        handle.write("\t".join(_text(row.get(column)) for column in columns) + "\n")
    return handle.getvalue().encode("utf-8")


def _sorted_views(
    state: CallsetState, refined: Mapping[str, tuple[int, int]] | None = None
) -> list[dict[str, Any]]:
    views = [_call_view(state.sample_id, row) for row in state.callset.rows]
    ids = [view["call_id"] for view in views]
    if len(set(ids)) != len(ids):
        raise InputContractError("native callset has duplicate event ids")
    if refined:
        unknown = sorted(set(refined) - set(ids))
        if unknown:
            raise InputContractError(f"refined breakpoints name unknown calls: {unknown[:5]}")
        for view in views:
            if view["call_id"] in refined:
                start, end = (int(x) for x in refined[view["call_id"]])
                if not view["start"] <= start < end <= view["end"]:
                    raise InputContractError(f"refined call {view['call_id']} leaves its native range")
                view.update(start=start, end=end, coordinate_kind=RANGE)
    return sorted(views, key=lambda v: (v["chrom"], v["start"], v["end"], v["family"], v["call_id"]))


def calls_tsv_bytes(
    state: CallsetState, refined: Mapping[str, tuple[int, int]] | None = None
) -> bytes:
    """``refined`` (call id -> [start, end)) replaces the native range of refined calls."""

    return _tsv(CALL_COLUMNS, _sorted_views(state, refined))


def calls_bed_bytes(
    state: CallsetState,
    contig_lengths: Mapping[str, int],
    refined: Mapping[str, tuple[int, int]] | None = None,
) -> bytes:
    rows = []
    for view in _sorted_views(state, refined):
        try:
            length = int(contig_lengths[view["chrom"]])
        except KeyError as exc:
            raise InputContractError(f"missing contig length for {view['chrom']}") from exc
        start, end = view["start"], view["end"]
        if view["coordinate_kind"] == POINT:
            if start < length:
                start0, end0 = start, start + 1
            elif start == length and length > 0:
                start0, end0 = start - 1, start
            else:
                raise InputContractError(f"call {view['call_id']} lies outside {view['chrom']}")
        else:
            start0, end0 = max(0, start), min(end, length)
            if start0 >= end0:
                raise InputContractError(f"range call {view['call_id']} lies outside {view['chrom']}")
        rows.append({
            "chrom": view["chrom"], "start0": start0, "end0": end0, "call_id": view["call_id"],
            "score": 0, "strand": ".", "family": view["family"],
            "coordinate_kind": view["coordinate_kind"], "call_origin": view["call_origin"],
        })
    return _tsv(BED_COLUMNS, rows)


def unresolved_tsv_bytes(state: CallsetState) -> bytes:
    rows = [
        {
            "sample_id": state.sample_id,
            "unresolved_id": f"{row.layer}:{row.event_key}",
            "stage": row.layer,
            "reason": row.reason,
            "chrom": row.chrom,
            "lo": row.lo,
            "hi": row.hi,
            "family": row.family,
            "participant_ids": sorted(row.participants),
        }
        for row in state.unresolved
    ]
    rows.sort(key=lambda r: (r["chrom"], r["lo"], r["hi"], r["stage"], r["unresolved_id"]))
    return _tsv(UNRESOLVED_COLUMNS, rows)


def native_rows(data: bytes) -> list[dict[str, str]]:
    """Parse native bytes (utility for tests and adapters)."""

    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""), delimiter="\t"))
