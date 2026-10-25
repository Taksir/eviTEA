"""Extension layers appended to the library-layer callset.

Four ordered extension layers are composed on top of the byte-exact
library-layer callset::

    library layer -> bilateral paired-end anchor support
        -> asymmetric paired-end anchor support
        -> full-memory paired-end anchor support
        -> internal 21-mer exact junction evidence

Every layer keeps the complete previous callset as an exact byte prefix.
It appends only events that the callset does not already represent, in the
schema of the library-layer callset.  The layers differ in candidate admission,
event identifiers, family alias handling, distance geometry, and collision
policy.  These differences are fixed, not options:

==============================  ===========  ===========  ===============  ===============
policy                          bilateral    asymmetric   full-memory      internal 21-mer
==============================  ===========  ===========  ===============  ===============
admitted paired-end states      BILATERAL    one side     any but NEITHER  (locus record)
family aliases                  no           no           yes              yes
distance geometry               point-range  point-range  point-range      range-range
nearby different-family call    ignored      ignored      ignored          unresolved
two new events within 500 bp    not checked  not checked  same family      any family
==============================  ===========  ===========  ===============  ===============

If two new events collide, neither is appended.  The pair is returned as an
``UnresolvedShadowEvent``.

This module reads no report files.  It takes only the typed inputs below.
"""

from __future__ import annotations

import ast
import csv
import io
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

from .constants import EVENT_EQUIV_DISTANCE, FROZEN_FAMILY_ALIASES, MIN_FOLDS, MIN_PAIRS, N_FOLDS
from .errors import InputContractError

REQUIRED_SHADOW_FIELDS = frozenset({"event_id", "chr", "start", "end", "family", "families"})
_SIDES = ("LEFT", "RIGHT")
_B_STATES = ("BILATERAL", "LEFT_ONLY", "RIGHT_ONLY", "NEITHER")


# ---------------------------------------------------------------------------
# The append-only callset
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ShadowCallset:
    """Byte-exact callset: every layer appends to ``data`` without rewriting it."""

    data: bytes
    fields: tuple[str, ...]
    rows: tuple[Mapping[str, str], ...]

    @classmethod
    def from_bytes(cls, data: bytes) -> "ShadowCallset":
        if not data:
            raise InputContractError("callset is empty")
        text = data.decode("utf-8-sig")
        reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        missing = REQUIRED_SHADOW_FIELDS - set(fields)
        if missing:
            raise InputContractError(f"callset schema lacks {sorted(missing)}")
        return cls(data, fields, tuple(dict(row) for row in reader))

    @property
    def row_count(self) -> int:
        return len(self.rows)

    def append(self, rows: Sequence[Mapping[str, str]]) -> "ShadowCallset":
        """Append rows in the base schema, keeping the base bytes as an exact prefix."""

        head = self.data if self.data.endswith((b"\n", b"\r")) else self.data + b"\n"
        handle = io.StringIO(newline="")
        writer = csv.DictWriter(
            handle, fieldnames=list(self.fields), delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in self.fields})
        data = head + handle.getvalue().encode("utf-8")
        if data[: len(self.data)] != self.data:  # pragma: no cover - structural guard
            raise InputContractError("appended callset lost its exact base prefix")
        return ShadowCallset.from_bytes(data)


# ---------------------------------------------------------------------------
# Typed layer inputs
# ---------------------------------------------------------------------------


def _ids(values: Iterable[int]) -> tuple[int, ...]:
    out = tuple(sorted({int(value) for value in values}))
    if any(value < 0 for value in out):
        raise InputContractError("pair identifiers must be non-negative")
    return out


def _side(value: str) -> str:
    if value not in _SIDES:
        raise InputContractError(f"anchor host side must be LEFT or RIGHT: {value!r}")
    return value


def _passes(pair_ids: Sequence[int]) -> bool:
    return len(pair_ids) >= MIN_PAIRS and len({p % N_FOLDS for p in pair_ids}) >= MIN_FOLDS


def paired_end_state(left_pair_ids: Sequence[int], right_pair_ids: Sequence[int]) -> str:
    left, right = _passes(left_pair_ids), _passes(right_pair_ids)
    if left and right:
        return "BILATERAL"
    if left:
        return "LEFT_ONLY"
    if right:
        return "RIGHT_ONLY"
    return "NEITHER"


@dataclass(frozen=True, slots=True)
class AnchorPECandidate:
    """One replicated exact-junction anchor of a corrected cluster, scored for paired-end anchor support."""

    cluster_id: int
    family: str
    chrom: str
    breakpoint0: int
    anchor_host_side: str
    a_pair_ids: tuple[int, ...]
    left_b_pair_ids: tuple[int, ...]
    right_b_pair_ids: tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "a_pair_ids", _ids(self.a_pair_ids))
        object.__setattr__(self, "left_b_pair_ids", _ids(self.left_b_pair_ids))
        object.__setattr__(self, "right_b_pair_ids", _ids(self.right_b_pair_ids))
        _side(self.anchor_host_side)
        if self.breakpoint0 < 0 or self.cluster_id < 0:
            raise InputContractError("paired-end anchor candidate coordinates must be non-negative")
        if not _passes(self.a_pair_ids):
            raise InputContractError("exact junction evidence requires >=2 pairs across >=2 folds")

    @property
    def b_state(self) -> str:
        return paired_end_state(self.left_b_pair_ids, self.right_b_pair_ids)

    @property
    def distinct_b_pair_ids(self) -> tuple[int, ...]:
        return _ids((*self.left_b_pair_ids, *self.right_b_pair_ids))

    @property
    def all_pair_ids(self) -> tuple[int, ...]:
        return _ids((*self.a_pair_ids, *self.left_b_pair_ids, *self.right_b_pair_ids))


@dataclass(frozen=True, slots=True)
class FullMemoryEvent:
    """Hostile-valid proposal-memory exact-junction event with its held-out paired-end anchor support."""

    candidate_id: str
    family: str
    chrom: str
    breakpoint0: int
    anchor_host_side: str
    a_pair_ids: tuple[int, ...]
    left_b_pair_ids: tuple[int, ...]
    right_b_pair_ids: tuple[int, ...]
    b_state: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "a_pair_ids", _ids(self.a_pair_ids))
        object.__setattr__(self, "left_b_pair_ids", _ids(self.left_b_pair_ids))
        object.__setattr__(self, "right_b_pair_ids", _ids(self.right_b_pair_ids))
        _side(self.anchor_host_side)
        if self.b_state not in _B_STATES:
            raise InputContractError(f"unknown paired-end state {self.b_state!r}")
        if not self.candidate_id:
            raise InputContractError("full-memory candidate_id must be non-empty")

    @property
    def b_pair_mentions(self) -> int:
        return len(self.left_b_pair_ids) + len(self.right_b_pair_ids)


@dataclass(frozen=True, slots=True)
class Internal21Instance:
    """One confirmed exact key of an internal 21-mer locus, in one rotation."""

    rotation: str
    breakpoint0: int
    anchor_host_side: str
    a_pair_ids: tuple[int, ...]
    left_b_pair_ids: tuple[int, ...]
    right_b_pair_ids: tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "a_pair_ids", _ids(self.a_pair_ids))
        object.__setattr__(self, "left_b_pair_ids", _ids(self.left_b_pair_ids))
        object.__setattr__(self, "right_b_pair_ids", _ids(self.right_b_pair_ids))
        _side(self.anchor_host_side)
        if not self.rotation:
            raise InputContractError("internal 21-mer instance rotation must be non-empty")


@dataclass(frozen=True, slots=True)
class Internal21Locus:
    """A ``NOVEL_UNCONFOUNDED`` internal 21-mer locus and its instances."""

    family: str
    chrom: str
    locus_lo: int
    locus_hi: int
    instances: tuple[Internal21Instance, ...]

    def __post_init__(self) -> None:
        if self.locus_lo > self.locus_hi:
            raise InputContractError("internal 21-mer locus bounds are inverted")
        object.__setattr__(self, "instances", tuple(self.instances))

    def representative(self) -> tuple[Internal21Instance, int]:
        """Return the representative exact key (chosen without truth data) and its rotation count.

        The key seen in the most confirmed rotations wins, then the lower
        breakpoint, then LEFT before RIGHT.  The instance returned has the
        lexicographically earliest rotation of that key.
        """

        members = [
            row for row in self.instances
            if self.locus_lo <= row.breakpoint0 <= self.locus_hi
        ]
        if not members:
            raise InputContractError(
                f"internal 21-mer locus {self.family} {self.chrom}:{self.locus_lo}-{self.locus_hi} has no instance"
            )
        by_key: dict[tuple[int, str], list[Internal21Instance]] = defaultdict(list)
        for row in members:
            by_key[(row.breakpoint0, row.anchor_host_side)].append(row)
        key = min(
            by_key,
            key=lambda k: (-len({r.rotation for r in by_key[k]}), k[0], 0 if k[1] == "LEFT" else 1),
        )
        chosen = min(by_key[key], key=lambda r: r.rotation)
        return chosen, len({r.rotation for r in by_key[key]})


# ---------------------------------------------------------------------------
# Shared geometry and family matching
# ---------------------------------------------------------------------------

_ALIASES = dict(FROZEN_FAMILY_ALIASES)


def _alias(name: str, use_aliases: bool) -> str:
    name = str(name or "").strip()
    return _ALIASES.get(name, name) if use_aliases else name


def call_families(row: Mapping[str, str], *, use_aliases: bool) -> frozenset[str]:
    """Return the family set of an existing call, from ``families`` or else ``family``."""

    raw = str(row.get("families", "") or "").strip()
    out: set[str] = set()
    if raw:
        try:
            parsed = ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            parsed = None
        if isinstance(parsed, (list, tuple, set)):
            out = {_alias(v, use_aliases) for v in parsed if str(v).strip()}
        if not out:
            out = {_alias(v, use_aliases) for v in re.split(r"[;,]", raw) if v.strip()}
    if not out and str(row.get("family", "") or "").strip():
        out = {_alias(row["family"], use_aliases)}
    return frozenset(out)


def _same_family(family: str, row: Mapping[str, str], use_aliases: bool) -> bool:
    fam = _alias(family, use_aliases)
    return fam in call_families(row, use_aliases=use_aliases) or fam == _alias(
        row.get("family", ""), use_aliases
    )


def range_distance(lo_a: int, hi_a: int, lo_b: int, hi_b: int) -> int:
    """Return the distance between two closed integer ranges (0 if they overlap)."""

    if hi_a < lo_b:
        return lo_b - hi_a
    if hi_b < lo_a:
        return lo_a - hi_b
    return 0


def _row_range(row: Mapping[str, str]) -> tuple[int, int] | None:
    try:
        return int(row["start"]), int(row["end"])
    except (KeyError, TypeError, ValueError):
        return None


def nearest_call(
    calls: Iterable[Mapping[str, str]],
    *,
    chrom: str,
    family: str,
    lo: int,
    hi: int,
    same_family: bool,
    use_aliases: bool,
    tie_break_event_id: bool,
) -> tuple[int | None, Mapping[str, str] | None]:
    """Return (distance, row) of the nearest same-family (or different-family) call on ``chrom`` to ``[lo, hi]``."""

    best: int | None = None
    best_row: Mapping[str, str] | None = None
    for row in calls:
        if row.get("chr") != chrom:
            continue
        if _same_family(family, row, use_aliases) is not same_family:
            continue
        span = _row_range(row)
        if span is None:
            continue
        distance = range_distance(lo, hi, span[0], span[1])
        if best is None or distance < best or (
            tie_break_event_id
            and distance == best
            and str(row.get("event_id", "")) < str(best_row.get("event_id", ""))  # type: ignore[union-attr]
        ):
            best, best_row = distance, row
    return best, best_row


# ---------------------------------------------------------------------------
# Row construction in the library-layer schema
# ---------------------------------------------------------------------------


def _call_row(
    *,
    event_id: str,
    family: str,
    chrom: str,
    breakpoint0: int,
    anchor_host_side: str,
    left_support: int,
    right_support: int,
    total_support: int,
    fragment_ids: Iterable[int],
    evidence_type: str,
    source: str,
    locus_id: str,
) -> dict[str, str]:
    return {
        "event_id": event_id,
        "chr": chrom,
        "start": str(breakpoint0),
        "end": str(breakpoint0),
        "family": family,
        "families": repr([family]),
        "left_coord": str(breakpoint0) if anchor_host_side == "LEFT" else "-1",
        "right_coord": str(breakpoint0) if anchor_host_side == "RIGHT" else "-1",
        "delta": "0",
        "left_support": str(left_support),
        "right_support": str(right_support),
        "total_support": str(total_support),
        "evidence_type": evidence_type,
        "score": str(total_support),
        "fragment_ids": ",".join(str(value) for value in _ids(fragment_ids)),
        "source": source,
        "route": source,
        "method": source,
        "call_source": source,
        "locus_id": locus_id,
        "support_pairs": str(total_support),
    }


# ---------------------------------------------------------------------------
# Layer results
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ShadowLedgerRow:
    layer: str
    event_key: str
    status: str
    event_id: str = ""
    family: str = ""
    chrom: str = ""
    breakpoint0: int | None = None
    nearest_same_family_distance: int | None = None
    nearest_same_family_event_id: str = ""
    nearest_different_family_distance: int | None = None
    nearest_different_family_event_id: str = ""
    detail: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class UnresolvedShadowEvent:
    layer: str
    event_key: str
    reason: str
    chrom: str
    family: str
    lo: int
    hi: int
    participants: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ShadowLayerResult:
    layer: str
    callset: ShadowCallset
    appended_event_ids: tuple[str, ...]
    ledger: tuple[ShadowLedgerRow, ...]
    unresolved: tuple[UnresolvedShadowEvent, ...]


ADD = "ADDED"
ABSORBED = "ABSORBED_EXISTING_SAME_FAMILY"
UNRESOLVED_DIFF = "UNRESOLVED_NEARBY_DIFFERENT_FAMILY_CALL"
UNRESOLVED_SAME_NEW = "UNRESOLVED_NEW_SAME_FAMILY_MULTIPLICITY"
UNRESOLVED_CROSS_NEW = "UNRESOLVED_NEW_CROSS_FAMILY_COMPETITION"


def _new_to_new(
    layer: str,
    pending: list[tuple[str, str, str, int, dict[str, str]]],
    *,
    same_family_only: bool,
    use_aliases: bool,
) -> tuple[list[tuple[str, str, str, int, dict[str, str]]], list[UnresolvedShadowEvent]]:
    """Exclude new events within 500 bp of another new event, and return them as unresolved."""

    conflicts: dict[int, set[int]] = defaultdict(set)
    for i, (_, fam_a, chrom_a, bp_a, _) in enumerate(pending):
        for j in range(i + 1, len(pending)):
            _, fam_b, chrom_b, bp_b, _ = pending[j]
            if chrom_a != chrom_b or abs(bp_a - bp_b) > EVENT_EQUIV_DISTANCE:
                continue
            same = _alias(fam_a, use_aliases) == _alias(fam_b, use_aliases)
            if same_family_only and not same:
                continue
            conflicts[i].add(j)
            conflicts[j].add(i)
    kept: list[tuple[str, str, str, int, dict[str, str]]] = []
    unresolved: list[UnresolvedShadowEvent] = []
    for i, item in enumerate(pending):
        if not conflicts[i]:
            kept.append(item)
            continue
        key, fam, chrom, bp, _ = item
        peers = sorted(conflicts[i])
        same = all(
            _alias(pending[j][1], use_aliases) == _alias(fam, use_aliases) for j in peers
        )
        unresolved.append(
            UnresolvedShadowEvent(
                layer=layer,
                event_key=key,
                reason=UNRESOLVED_SAME_NEW if same else UNRESOLVED_CROSS_NEW,
                chrom=chrom,
                family=fam,
                lo=bp,
                hi=bp,
                participants=tuple(pending[j][0] for j in peers),
            )
        )
    return kept, unresolved


# ---------------------------------------------------------------------------
# Layers
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ShadowLabels:
    """Sample labels embedded in event identifiers, one per layer.

    ``for_sample`` builds the default labels for a sample.
    """

    bilateral_prefix: str
    asymmetric: str
    full_memory: str
    internal21: str

    @classmethod
    def for_sample(cls, sample_id: str) -> "ShadowLabels":
        tag = re.sub(r"[^A-Za-z0-9]+", "", sample_id).upper() or "SAMPLE"
        return cls(f"{tag}_", tag, tag, tag)


def _anchor_pe_layer(
    layer: str,
    base: ShadowCallset,
    candidates: Iterable[AnchorPECandidate],
    admitted_states: frozenset[str],
    event_id: callable,
    evidence_type: str,
    source: str,
    locus_id: callable,
) -> ShadowLayerResult:
    admitted = [row for row in candidates if row.b_state in admitted_states]
    ordered = sorted(
        admitted, key=lambda r: (r.chrom, r.breakpoint0, r.family, r.cluster_id)
    )
    keys = {(r.chrom, r.breakpoint0, r.family) for r in ordered}
    if len(keys) != len(ordered):
        raise InputContractError(f"{layer}: duplicate (chrom, breakpoint, family) key")
    existing_ids = {row.get("event_id", "") for row in base.rows}
    ledger: list[ShadowLedgerRow] = []
    new_rows: list[dict[str, str]] = []
    for index, row in enumerate(ordered, start=1):
        eid = event_id(row, index)
        if eid in existing_ids:
            raise InputContractError(f"{layer}: event_id collision {eid}")
        distance, nearest = nearest_call(
            base.rows, chrom=row.chrom, family=row.family, lo=row.breakpoint0,
            hi=row.breakpoint0, same_family=True, use_aliases=False,
            tie_break_event_id=False,
        )
        represented = distance is not None and distance <= EVENT_EQUIV_DISTANCE
        ledger.append(
            ShadowLedgerRow(
                layer=layer,
                event_key=f"{row.family}:{row.chrom}:{row.breakpoint0}:{row.anchor_host_side}",
                status=ABSORBED if represented else ADD,
                event_id=eid,
                family=row.family,
                chrom=row.chrom,
                breakpoint0=row.breakpoint0,
                nearest_same_family_distance=distance,
                nearest_same_family_event_id="" if nearest is None else nearest.get("event_id", ""),
                detail={"cluster_id": row.cluster_id, "b_state": row.b_state},
            )
        )
        if represented:
            continue
        distinct_b = row.distinct_b_pair_ids
        new_rows.append(
            _call_row(
                event_id=eid, family=row.family, chrom=row.chrom,
                breakpoint0=row.breakpoint0, anchor_host_side=row.anchor_host_side,
                left_support=len(row.left_b_pair_ids), right_support=len(row.right_b_pair_ids),
                total_support=len(distinct_b), fragment_ids=row.all_pair_ids,
                evidence_type=evidence_type, source=source, locus_id=locus_id(row),
            )
        )
    callset = base.append(new_rows)
    return ShadowLayerResult(
        layer, callset, tuple(r["event_id"] for r in new_rows), tuple(ledger), ()
    )


def append_bilateral_layer(
    base: ShadowCallset, candidates: Iterable[AnchorPECandidate], labels: ShadowLabels
) -> ShadowLayerResult:
    """Append events with exact junction evidence and paired-end anchor support on both sides."""

    prefix = labels.bilateral_prefix
    return _anchor_pe_layer(
        "BILATERAL", base, candidates, frozenset({"BILATERAL"}),
        lambda r, i: f"V3ANCHORPE_{prefix}{r.cluster_id:04d}_{i:03d}",
        "v3_anchor_pe_bilateral", "V3_ANCHOR_PE",
        lambda r: f"V3AP_{prefix}CLUSTER_{r.cluster_id}",
    )


def append_asymmetric_layer(
    base: ShadowCallset, candidates: Iterable[AnchorPECandidate], labels: ShadowLabels
) -> ShadowLayerResult:
    """Append events with exact junction evidence and paired-end anchor support on exactly one side."""

    tag = labels.asymmetric
    return _anchor_pe_layer(
        "ASYMMETRIC", base, candidates, frozenset({"LEFT_ONLY", "RIGHT_ONLY"}),
        lambda r, i: f"V3ANCHORPE_ASYM_{tag}_{r.cluster_id:04d}_{i:03d}",
        "v3_anchor_pe_asymmetric_shadow", "V3_ANCHOR_PE_ASYMMETRIC_SHADOW",
        lambda r: f"V3AP_ASYM_{tag}_CLUSTER_{r.cluster_id}",
    )


def append_full_memory_layer(
    base: ShadowCallset,
    events: Iterable[FullMemoryEvent],
    labels: ShadowLabels,
    *,
    sample_id: str,
) -> ShadowLayerResult:
    """Append full-memory events (exact junction evidence, paired-end anchor support on at least one side)."""

    layer = "FULL_MEMORY"
    valid = [row for row in events if row.b_state != "NEITHER"]
    ordered = sorted(
        valid,
        key=lambda r: (r.chrom, _alias(r.family, True), r.breakpoint0, r.candidate_id),
    )
    ledger: list[ShadowLedgerRow] = []
    pending: list[tuple[str, str, str, int, dict[str, str]]] = []
    for row in ordered:
        family = _alias(row.family, True)
        eid = f"V3FULLMEM_AP_{labels.full_memory}_{row.candidate_id.split('_')[-1]}"
        distance, nearest = nearest_call(
            base.rows, chrom=row.chrom, family=family, lo=row.breakpoint0,
            hi=row.breakpoint0, same_family=True, use_aliases=True,
            tie_break_event_id=False,
        )
        represented = distance is not None and distance <= EVENT_EQUIV_DISTANCE
        ledger.append(
            ShadowLedgerRow(
                layer=layer,
                event_key=row.candidate_id,
                status=ABSORBED if represented else ADD,
                event_id=eid, family=family, chrom=row.chrom, breakpoint0=row.breakpoint0,
                nearest_same_family_distance=distance,
                nearest_same_family_event_id="" if nearest is None else nearest.get("event_id", ""),
                detail={"sample_id": sample_id, "b_state": row.b_state},
            )
        )
        if represented:
            continue
        pending.append(
            (
                eid, family, row.chrom, row.breakpoint0,
                _call_row(
                    event_id=eid, family=family, chrom=row.chrom,
                    breakpoint0=row.breakpoint0, anchor_host_side=row.anchor_host_side,
                    left_support=len(row.left_b_pair_ids),
                    right_support=len(row.right_b_pair_ids),
                    total_support=row.b_pair_mentions,
                    fragment_ids=(*row.a_pair_ids, *row.left_b_pair_ids, *row.right_b_pair_ids),
                    evidence_type="v3_full_memory_anchor_pe_shadow",
                    source="V3_FULL_MEMORY_ANCHOR_PE_SHADOW",
                    locus_id=f"V3FMAP_{sample_id}_{row.candidate_id}",
                ),
            )
        )
    kept, unresolved = _new_to_new(layer, pending, same_family_only=True, use_aliases=True)
    callset = base.append([item[4] for item in kept])
    return ShadowLayerResult(
        layer, callset, tuple(item[0] for item in kept), tuple(ledger), tuple(unresolved)
    )


def append_internal21_layer(
    base: ShadowCallset, loci: Iterable[Internal21Locus], labels: ShadowLabels
) -> ShadowLayerResult:
    """Append internal 21-mer loci that make exact junction evidence observable."""

    layer = "INTERNAL21"
    ordered = sorted(
        loci, key=lambda r: (r.chrom, _alias(r.family, True), r.locus_lo, r.locus_hi)
    )
    ledger: list[ShadowLedgerRow] = []
    unresolved: list[UnresolvedShadowEvent] = []
    pending: list[tuple[str, str, str, int, dict[str, str]]] = []
    for locus in ordered:
        family = _alias(locus.family, True)
        key = f"{family}:{locus.chrom}:{locus.locus_lo}-{locus.locus_hi}"
        common = dict(chrom=locus.chrom, family=family, lo=locus.locus_lo, hi=locus.locus_hi,
                      use_aliases=True, tie_break_event_id=True)
        same_d, same_r = nearest_call(base.rows, same_family=True, **common)
        diff_d, diff_r = nearest_call(base.rows, same_family=False, **common)
        nearest = dict(
            nearest_same_family_distance=same_d,
            nearest_same_family_event_id="" if same_r is None else same_r.get("event_id", ""),
            nearest_different_family_distance=diff_d,
            nearest_different_family_event_id="" if diff_r is None else diff_r.get("event_id", ""),
        )
        if same_d is not None and same_d <= EVENT_EQUIV_DISTANCE:
            ledger.append(ShadowLedgerRow(layer, key, ABSORBED, family=family, chrom=locus.chrom, **nearest))
            continue
        if diff_d is not None and diff_d <= EVENT_EQUIV_DISTANCE:
            ledger.append(ShadowLedgerRow(layer, key, UNRESOLVED_DIFF, family=family, chrom=locus.chrom, **nearest))
            unresolved.append(
                UnresolvedShadowEvent(
                    layer, key, UNRESOLVED_DIFF, locus.chrom, family, locus.locus_lo,
                    locus.locus_hi, (nearest["nearest_different_family_event_id"],),
                )
            )
            continue
        chosen, recurrence = locus.representative()
        ledger.append(
            ShadowLedgerRow(
                layer, key, ADD, family=family, chrom=locus.chrom,
                breakpoint0=chosen.breakpoint0,
                detail={"rotation": chosen.rotation, "exact_key_rotation_count": recurrence},
                **nearest,
            )
        )
        a, left, right = chosen.a_pair_ids, chosen.left_b_pair_ids, chosen.right_b_pair_ids
        total = len(set(a) | set(left) | set(right))
        pending.append(
            (
                key, family, locus.chrom, chosen.breakpoint0,
                dict(
                    family=family, chrom=locus.chrom, breakpoint0=chosen.breakpoint0,
                    anchor_host_side=chosen.anchor_host_side, left_support=len(left),
                    right_support=len(right), total_support=total,
                    fragment_ids=(*a, *left, *right),
                    evidence_type="v3_internal21_strict_a_anchor_pe_shadow",
                    source="V3_INTERNAL21_STRICT_A_ANCHORPE_SHADOW",
                ),
            )
        )
    kept, new_unresolved = _new_to_new(layer, pending, same_family_only=False, use_aliases=True)
    unresolved.extend(new_unresolved)
    rows: list[dict[str, str]] = []
    event_ids: list[str] = []
    for index, (_, _, _, _, spec) in enumerate(kept, start=1):
        eid = f"V3INTA_AP_{labels.internal21}_{index:06d}"
        rows.append(_call_row(event_id=eid, locus_id=eid, **spec))
        event_ids.append(eid)
    callset = base.append(rows)
    return ShadowLayerResult(layer, callset, tuple(event_ids), tuple(ledger), tuple(unresolved))


@dataclass(frozen=True, slots=True)
class ExtensionLayerInputs:
    bilateral_candidates: tuple[AnchorPECandidate, ...]
    full_memory_events: tuple[FullMemoryEvent, ...]
    internal21_loci: tuple[Internal21Locus, ...]


def compose_frozen_shadow(
    incumbent: bytes,
    inputs: ExtensionLayerInputs,
    labels: ShadowLabels,
    *,
    sample_id: str,
) -> tuple[ShadowLayerResult, ...]:
    """Apply the four extension layers in order and return each layer's result.

    The bilateral and asymmetric layers share one candidate set.  Each admits
    candidates by paired-end state.
    """

    base = ShadowCallset.from_bytes(incumbent)
    bilateral = append_bilateral_layer(base, inputs.bilateral_candidates, labels)
    asymmetric = append_asymmetric_layer(bilateral.callset, inputs.bilateral_candidates, labels)
    full_memory = append_full_memory_layer(
        asymmetric.callset, inputs.full_memory_events, labels, sample_id=sample_id
    )
    internal21 = append_internal21_layer(full_memory.callset, inputs.internal21_loci, labels)
    return bilateral, asymmetric, full_memory, internal21
