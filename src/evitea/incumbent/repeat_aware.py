"""Credit recurrent edits to retain candidate-family groups, and build repeat-aware events per locus."""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

from ..constants import repeat_aware_edit_rate_passes
from ..contracts import BreakpointSide, ClosedBreakpointInterval0, RepeatAwareEvent
from ..errors import InputContractError
from ..families import canonical_family

MIN_OTHER_PHYSICAL_GROUPS = 2
_CIGAR_TOKEN = re.compile(r"([1-9][0-9]*)([MIDNSHP=X])")
_VARIANT_KINDS = frozenset({"SUB", "INS", "DEL"})


def _nonempty(value: str, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise InputContractError(f"{field} must be non-empty text")
    return value


def _nonnegative(value: int, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise InputContractError(f"{field} must be a non-negative integer")
    return value


def _positive(value: int, field: str) -> int:
    _nonnegative(value, field)
    if value == 0:
        raise InputContractError(f"{field} must be positive")
    return value


def _canonical_pair_ids(values: Iterable[int], field: str) -> tuple[int, ...]:
    rows = tuple(values)
    for value in rows:
        _nonnegative(value, field)
    return tuple(sorted(set(rows)))


@dataclass(frozen=True, slots=True, order=True)
class RepeatAwareVariant:
    """One candidate-family alignment edit in TE-consensus coordinates."""

    kind: str
    ref_pos0: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        if self.kind not in _VARIANT_KINDS:
            raise InputContractError(
                "repeat-aware variant kind must be SUB, INS, or DEL"
            )
        _nonnegative(self.ref_pos0, "ref_pos0")
        if not isinstance(self.ref, str) or not isinstance(self.alt, str):
            raise InputContractError("repeat-aware variant alleles must be text")
        if self.kind == "SUB" and (len(self.ref) != 1 or len(self.alt) != 1):
            raise InputContractError("SUB requires one reference and alternate base")
        if self.kind == "INS" and (self.ref or not self.alt):
            raise InputContractError("INS requires an empty ref and non-empty alt")
        if self.kind == "DEL" and (not self.ref or self.alt):
            raise InputContractError("DEL requires a non-empty ref and empty alt")

    @property
    def cost(self) -> int:
        if self.kind == "SUB":
            return 1
        return len(self.alt) if self.kind == "INS" else len(self.ref)

    @property
    def exact_key(self) -> tuple[str, int, str, str]:
        return (self.kind, self.ref_pos0, self.ref, self.alt)

    @property
    def site_key(self) -> tuple[str, int, int]:
        span = len(self.ref) if self.kind == "DEL" else 1
        return (self.kind, self.ref_pos0, span)

    @property
    def label(self) -> str:
        return f"{self.kind}@{self.ref_pos0}:{self.ref}>{self.alt}"


def _cigar_tokens(cigar: str) -> tuple[tuple[int, str], ...]:
    if not isinstance(cigar, str) or not cigar or cigar == "*":
        raise InputContractError(
            "repeat-aware alignment CIGAR must be mapped SAM text"
        )
    tokens: list[tuple[int, str]] = []
    cursor = 0
    for match in _CIGAR_TOKEN.finditer(cigar):
        if match.start() != cursor:
            raise InputContractError(f"invalid repeat-aware CIGAR: {cigar}")
        tokens.append((int(match.group(1)), match.group(2)))
        cursor = match.end()
    if cursor != len(cigar) or not tokens:
        raise InputContractError(f"invalid repeat-aware CIGAR: {cigar}")
    return tuple(tokens)


def alignment_variant_ops(
    query_sequence: str,
    reference_sequence: str,
    reference_start0: int,
    cigar: str,
) -> tuple[RepeatAwareVariant, ...]:
    """Reconstruct substitution/indel operations from one TE alignment."""

    query = _nonempty(query_sequence, "query_sequence").upper()
    reference = _nonempty(reference_sequence, "reference_sequence").upper()
    _nonnegative(reference_start0, "reference_start0")
    query0 = 0
    reference0 = reference_start0
    operations: list[RepeatAwareVariant] = []
    for length, operation in _cigar_tokens(cigar):
        if operation in {"M", "=", "X"}:
            if (
                query0 + length > len(query)
                or reference0 + length > len(reference)
            ):
                raise InputContractError(
                    "repeat-aware alignment extends beyond a sequence"
                )
            for offset in range(length):
                ref_base = reference[reference0 + offset]
                query_base = query[query0 + offset]
                if ref_base != query_base:
                    operations.append(
                        RepeatAwareVariant(
                            "SUB", reference0 + offset, ref_base, query_base
                        )
                    )
            query0 += length
            reference0 += length
        elif operation == "I":
            if query0 + length > len(query):
                raise InputContractError(
                    "repeat-aware insertion extends beyond query"
                )
            operations.append(
                RepeatAwareVariant(
                    "INS", reference0, "", query[query0 : query0 + length]
                )
            )
            query0 += length
        elif operation == "D":
            if reference0 + length > len(reference):
                raise InputContractError(
                    "repeat-aware deletion extends beyond reference"
                )
            operations.append(
                RepeatAwareVariant(
                    "DEL",
                    reference0,
                    reference[reference0 : reference0 + length],
                    "",
                )
            )
            reference0 += length
        elif operation == "N":
            reference0 += length
            if reference0 > len(reference):
                raise InputContractError(
                    "repeat-aware skip extends beyond reference"
                )
        elif operation == "S":
            query0 += length
            if query0 > len(query):
                raise InputContractError(
                    "repeat-aware soft clip extends beyond query"
                )
        elif operation in {"H", "P"}:
            continue
    if query0 != len(query):
        raise InputContractError(
            "repeat-aware CIGAR does not consume the complete query"
        )
    return tuple(operations)


def equivalent_indel_positions(
    reference_sequence: str,
    operation: RepeatAwareVariant,
) -> tuple[int, ...]:
    """Return all repeat-shift-equivalent TE-consensus positions."""

    reference = _nonempty(reference_sequence, "reference_sequence").upper()
    if operation.kind not in {"INS", "DEL"}:
        return (operation.ref_pos0,)
    positions = {operation.ref_pos0}
    if operation.kind == "INS":
        if operation.ref_pos0 > len(reference):
            raise InputContractError("insertion position lies beyond reference")
        position = operation.ref_pos0
        sequence = operation.alt
        while position > 0 and reference[position - 1] == sequence[-1]:
            position -= 1
            sequence = sequence[-1] + sequence[:-1]
            positions.add(position)
        position = operation.ref_pos0
        sequence = operation.alt
        while position < len(reference) and reference[position] == sequence[0]:
            sequence = sequence[1:] + sequence[0]
            position += 1
            positions.add(position)
    else:
        length = len(operation.ref)
        if operation.ref_pos0 + length > len(reference):
            raise InputContractError("deletion lies beyond reference")
        position = operation.ref_pos0
        deleted = operation.ref
        while position > 0 and reference[position - 1] == deleted[-1]:
            position -= 1
            deleted = reference[position] + deleted[:-1]
            positions.add(position)
        position = operation.ref_pos0
        deleted = operation.ref
        while (
            position + length < len(reference)
            and reference[position + length] == deleted[0]
        ):
            deleted = deleted[1:] + reference[position + length]
            position += 1
            positions.add(position)
    return tuple(sorted(positions))


def indel_haplotype_admissibility(
    reference_sequence: str,
    operation: RepeatAwareVariant,
    recurrent_operations: Sequence[RepeatAwareVariant],
) -> tuple[bool, str, tuple[int, ...], tuple[str, ...]]:
    """Apply the repeat-shift and linked-substitution safeguard."""

    if operation.kind not in {"INS", "DEL"}:
        return True, "NON_INDEL", (operation.ref_pos0,), ()
    positions = equivalent_indel_positions(reference_sequence, operation)
    if len(positions) == 1:
        return True, "PASS_UNIQUE_INDEL_PLACEMENT", positions, ()
    lo0 = min(positions)
    hi0 = max(positions) + max(0, len(operation.ref) - 1)
    anchor_lo0 = max(0, lo0 - 1)
    anchor_hi0 = min(len(reference_sequence) - 1, hi0 + 1)
    anchors = tuple(
        recurrent.label
        for recurrent in recurrent_operations
        if recurrent.kind == "SUB"
        and anchor_lo0 <= recurrent.ref_pos0 <= anchor_hi0
    )
    if anchors:
        return (
            True,
            "PASS_REPEAT_EQUIVALENT_INDEL_LINKED_SUB_ANCHOR",
            positions,
            anchors,
        )
    return (
        False,
        "REJECT_REPEAT_EQUIVALENT_INDEL_UNANCHORED",
        positions,
        (),
    )


@dataclass(frozen=True, slots=True)
class RepeatAwareGroupObservation:
    group_id: str
    locus_id: str
    side: BreakpointSide
    candidate_family: str
    representative_pair_id: int
    te_facing_length: int
    family_edit_costs: Mapping[str, int]
    te_record_id: str | None = None
    alignment_nm: int | None = None
    operations: tuple[RepeatAwareVariant, ...] = ()

    def __post_init__(self) -> None:
        _nonempty(self.group_id, "group_id")
        _nonempty(self.locus_id, "locus_id")
        if not isinstance(self.side, BreakpointSide):
            raise InputContractError(
                "repeat-aware group side must be LEFT or RIGHT"
            )
        object.__setattr__(
            self, "candidate_family", canonical_family(self.candidate_family)
        )
        _nonnegative(self.representative_pair_id, "representative_pair_id")
        _positive(self.te_facing_length, "te_facing_length")
        costs: dict[str, int] = {}
        for family, cost in self.family_edit_costs.items():
            canonical = canonical_family(family)
            _nonnegative(cost, "family edit cost")
            previous = costs.get(canonical)
            if previous is not None and previous != cost:
                raise InputContractError("aliased family costs disagree")
            costs[canonical] = cost
        if self.candidate_family not in costs:
            raise InputContractError("candidate family is absent from edit costs")
        object.__setattr__(
            self, "family_edit_costs", MappingProxyType(dict(sorted(costs.items())))
        )
        if (self.te_record_id is None) != (self.alignment_nm is None):
            raise InputContractError(
                "profile record and edit count must appear together"
            )
        if self.te_record_id is None:
            if self.operations:
                raise InputContractError(
                    "variant operations require an alignment profile"
                )
        else:
            _nonempty(self.te_record_id, "te_record_id")
            _nonnegative(self.alignment_nm, "alignment_nm")
        if not isinstance(self.operations, tuple) or not all(
            isinstance(operation, RepeatAwareVariant)
            for operation in self.operations
        ):
            raise InputContractError(
                "operations must be RepeatAwareVariant values"
            )

    @property
    def unique_best_family(self) -> str | None:
        best_cost = min(self.family_edit_costs.values())
        winners = tuple(
            family
            for family, cost in self.family_edit_costs.items()
            if cost == best_cost
        )
        return winners[0] if len(winners) == 1 else None

    @property
    def candidate_nm(self) -> int:
        return self.family_edit_costs[self.candidate_family]

    @property
    def baseline_pass(self) -> bool:
        return (
            self.unique_best_family == self.candidate_family
            and repeat_aware_edit_rate_passes(
                self.candidate_nm, self.te_facing_length
            )
        )

    @property
    def rate_only_failure(self) -> bool:
        return (
            self.unique_best_family == self.candidate_family
            and not repeat_aware_edit_rate_passes(
                self.candidate_nm, self.te_facing_length
            )
        )

    @property
    def profile_complete(self) -> bool:
        return self.te_record_id is not None and self.alignment_nm is not None

    @property
    def valid_nm_reconstruction(self) -> bool:
        return (
            self.profile_complete
            and self.alignment_nm == self.candidate_nm
            and sum(operation.cost for operation in self.operations)
            == self.alignment_nm
            and len({operation.exact_key for operation in self.operations})
            == len(self.operations)
        )


@dataclass(frozen=True, slots=True)
class RepeatAwareLocusObservation:
    sample_id: str
    locus_id: str
    legacy_event_id: str
    family: str
    chrom: str
    interval0: ClosedBreakpointInterval0
    existing_left_pair_ids: tuple[int, ...] = ()
    existing_right_pair_ids: tuple[int, ...] = ()
    left_coord0: int | None = None
    right_coord0: int | None = None

    def __post_init__(self) -> None:
        _nonempty(self.sample_id, "sample_id")
        _nonempty(self.locus_id, "locus_id")
        _nonempty(self.legacy_event_id, "legacy_event_id")
        object.__setattr__(self, "family", canonical_family(self.family))
        _nonempty(self.chrom, "chrom")
        if not isinstance(self.interval0, ClosedBreakpointInterval0):
            raise InputContractError("interval0 must be a closed interbase interval")
        object.__setattr__(
            self,
            "existing_left_pair_ids",
            _canonical_pair_ids(self.existing_left_pair_ids, "left pair_id"),
        )
        object.__setattr__(
            self,
            "existing_right_pair_ids",
            _canonical_pair_ids(self.existing_right_pair_ids, "right pair_id"),
        )
        for value, field in (
            (self.left_coord0, "left_coord0"),
            (self.right_coord0, "right_coord0"),
        ):
            if value is not None:
                _nonnegative(value, field)


@dataclass(frozen=True, slots=True)
class RepeatAwareSampleContext:
    sample_id: str
    loci: tuple[RepeatAwareLocusObservation, ...]

    def __post_init__(self) -> None:
        _nonempty(self.sample_id, "sample_id")
        if any(locus.sample_id != self.sample_id for locus in self.loci):
            raise InputContractError("repeat-aware loci mix sample identity")


@dataclass(frozen=True, slots=True)
class RepeatAwarePreparedAssets:
    te_sequences: Mapping[str, str]

    def __post_init__(self) -> None:
        sequences: dict[str, str] = {}
        for record_id, sequence in self.te_sequences.items():
            _nonempty(record_id, "TE record ID")
            sequences[record_id] = _nonempty(sequence, "TE sequence").upper()
        object.__setattr__(
            self, "te_sequences", MappingProxyType(dict(sorted(sequences.items())))
        )


@dataclass(frozen=True, slots=True)
class RepeatAwareObservations:
    groups: tuple[RepeatAwareGroupObservation, ...]
    aggregation_group_ids: frozenset[str] | None = None

    def __post_init__(self) -> None:
        if self.aggregation_group_ids is None:
            object.__setattr__(
                self,
                "aggregation_group_ids",
                frozenset(group.group_id for group in self.groups),
            )
            return
        values = frozenset(self.aggregation_group_ids)
        if not all(isinstance(value, str) and value for value in values):
            raise InputContractError(
                "aggregation_group_ids must contain non-empty group IDs"
            )
        object.__setattr__(self, "aggregation_group_ids", values)


@dataclass(frozen=True, slots=True)
class RepeatAwareGroupDecision:
    group_id: str
    baseline_pass: bool
    eligible_rate_only_failure: bool
    crossvalidated_pass: bool
    status: str
    adjusted_nm: int
    recurrent_variant_labels: tuple[str, ...] = ()
    credited_variant_labels: tuple[str, ...] = ()
    training_pair_ids: tuple[int, ...] = ()

    @property
    def final_pass(self) -> bool:
        return self.baseline_pass or self.crossvalidated_pass


@dataclass(frozen=True, slots=True)
class RepeatAwareLocusDecision:
    locus_id: str
    left_pair_ids: tuple[int, ...]
    right_pair_ids: tuple[int, ...]
    passes_bilateral_support: bool


@dataclass(frozen=True, slots=True)
class RepeatAwareResult:
    events: tuple[RepeatAwareEvent, ...]
    group_decisions: tuple[RepeatAwareGroupDecision, ...]
    locus_decisions: tuple[RepeatAwareLocusDecision, ...]


def _operations_match_reference(
    reference: str,
    operations: Sequence[RepeatAwareVariant],
) -> bool:
    for operation in operations:
        if operation.kind == "INS":
            if operation.ref_pos0 > len(reference):
                return False
        elif operation.kind == "SUB":
            if (
                operation.ref_pos0 >= len(reference)
                or reference[operation.ref_pos0] != operation.ref
            ):
                return False
        elif (
            operation.ref_pos0 + len(operation.ref) > len(reference)
            or reference[
                operation.ref_pos0 : operation.ref_pos0 + len(operation.ref)
            ]
            != operation.ref
        ):
            return False
    return True


def _profile_usable(
    group: RepeatAwareGroupObservation,
    assets: RepeatAwarePreparedAssets,
) -> bool:
    if (
        not group.valid_nm_reconstruction
        or group.te_record_id not in assets.te_sequences
    ):
        return False
    return _operations_match_reference(
        assets.te_sequences[group.te_record_id], group.operations
    )


def _admissible_variants(
    pool: Iterable[RepeatAwareGroupObservation],
) -> dict[tuple[str, int, str, str], set[int]]:
    site_alleles: dict[
        tuple[str, int, int],
        dict[tuple[str, int, str, str], set[int]],
    ] = defaultdict(lambda: defaultdict(set))
    for training in pool:
        for operation in training.operations:
            site_alleles[operation.site_key][operation.exact_key].add(
                training.representative_pair_id
            )
    admissible: dict[tuple[str, int, str, str], set[int]] = {}
    for alleles in site_alleles.values():
        ranked = sorted(alleles.items(), key=lambda item: (-len(item[1]), item[0]))
        if not ranked or len(ranked[0][1]) < MIN_OTHER_PHYSICAL_GROUPS:
            continue
        if len(ranked) > 1 and len(ranked[1][1]) == len(ranked[0][1]):
            continue
        admissible[ranked[0][0]] = ranked[0][1]
    return admissible


def _crossvalidate_group(
    group: RepeatAwareGroupObservation,
    groups: Sequence[RepeatAwareGroupObservation],
    assets: RepeatAwarePreparedAssets,
) -> RepeatAwareGroupDecision:
    if group.baseline_pass:
        return RepeatAwareGroupDecision(
            group.group_id,
            True,
            False,
            False,
            "PASS_BASELINE_UNIQUE_BEST_RATE_3_PER_35",
            group.candidate_nm,
        )
    if not group.rate_only_failure:
        return RepeatAwareGroupDecision(
            group.group_id,
            False,
            False,
            False,
            "NOT_ELIGIBLE_UNIQUE_BEST_RATE_ONLY_FAILURE",
            group.candidate_nm,
        )
    if not _profile_usable(group, assets):
        return RepeatAwareGroupDecision(
            group.group_id,
            False,
            True,
            False,
            "NO_VALID_ALIGNMENT_PROFILE",
            group.candidate_nm,
        )
    pool = tuple(
        candidate
        for candidate in groups
        if candidate.group_id != group.group_id
        and candidate.representative_pair_id != group.representative_pair_id
        and candidate.locus_id == group.locus_id
        and candidate.candidate_family == group.candidate_family
        and candidate.te_record_id == group.te_record_id
        and _profile_usable(candidate, assets)
    )
    admissible = _admissible_variants(pool)
    heldout = {operation.exact_key: operation for operation in group.operations}
    recurrent_keys = tuple(sorted(key for key in heldout if key in admissible))
    recurrent_operations = tuple(heldout[key] for key in recurrent_keys)
    credited_keys: list[tuple[str, int, str, str]] = []
    for key in recurrent_keys:
        operation = heldout[key]
        admissible_indel, _status, _positions, _anchors = (
            indel_haplotype_admissibility(
                assets.te_sequences[group.te_record_id],
                operation,
                recurrent_operations,
            )
        )
        if admissible_indel:
            credited_keys.append(key)
    credit = sum(heldout[key].cost for key in credited_keys)
    adjusted_nm = max(0, group.alignment_nm - credit)
    rate_pass = repeat_aware_edit_rate_passes(
        adjusted_nm, group.te_facing_length
    )
    crossvalidated_pass = bool(credited_keys) and rate_pass
    if crossvalidated_pass:
        status = "PASS_REPEAT_AWARE_COHERENT_DIVERGENCE"
    elif not recurrent_keys:
        status = "NO_RECURRENT_VARIANT_IN_HELDOUT"
    elif not credited_keys:
        status = "ALL_RECURRENT_CREDIT_REMOVED_BY_INDEL_ANCHOR_GATE"
    else:
        status = "ADJUSTED_RATE_STILL_FAILS"
    training_pair_ids = tuple(
        sorted(
            {
                pair_id
                for key in credited_keys
                for pair_id in admissible[key]
            }
        )
    )
    return RepeatAwareGroupDecision(
        group.group_id,
        False,
        True,
        crossvalidated_pass,
        status,
        adjusted_nm,
        tuple(heldout[key].label for key in recurrent_keys),
        tuple(heldout[key].label for key in credited_keys),
        training_pair_ids,
    )


def run_repeat_aware_incumbent(
    sample_context: RepeatAwareSampleContext,
    prepared_assets: RepeatAwarePreparedAssets,
    shared_observations: RepeatAwareObservations,
) -> RepeatAwareResult:
    """Run the repeat-aware core and preserve event ranges."""

    if not isinstance(sample_context, RepeatAwareSampleContext):
        raise InputContractError(
            "sample_context must be RepeatAwareSampleContext"
        )
    if not isinstance(prepared_assets, RepeatAwarePreparedAssets):
        raise InputContractError(
            "prepared_assets must be RepeatAwarePreparedAssets"
        )
    if not isinstance(shared_observations, RepeatAwareObservations):
        raise InputContractError(
            "shared_observations must be RepeatAwareObservations"
        )
    loci = tuple(sample_context.loci)
    groups = tuple(shared_observations.groups)
    if len({locus.locus_id for locus in loci}) != len(loci):
        raise InputContractError("duplicate repeat-aware locus_id")
    if len({locus.legacy_event_id for locus in loci}) != len(loci):
        raise InputContractError("duplicate legacy_event_id")
    if len({group.group_id for group in groups}) != len(groups):
        raise InputContractError("duplicate repeat-aware group_id")
    group_ids = {group.group_id for group in groups}
    unknown_aggregation_ids = (
        shared_observations.aggregation_group_ids - group_ids
    )
    if unknown_aggregation_ids:
        raise InputContractError(
            "aggregation_group_ids refer to unknown repeat-aware groups"
        )
    loci_by_id = {locus.locus_id: locus for locus in loci}
    for group in groups:
        try:
            locus = loci_by_id[group.locus_id]
        except KeyError as exc:
            raise InputContractError(
                "repeat-aware group refers to an unknown locus"
            ) from exc
        if group.candidate_family != locus.family:
            raise InputContractError(
                "repeat-aware group and locus families disagree"
            )

    decisions = tuple(
        _crossvalidate_group(group, groups, prepared_assets)
        for group in sorted(groups, key=lambda item: item.group_id)
    )
    decision_by_group = {decision.group_id: decision for decision in decisions}
    groups_by_locus: dict[str, list[RepeatAwareGroupObservation]] = defaultdict(list)
    for group in groups:
        groups_by_locus[group.locus_id].append(group)

    events: list[RepeatAwareEvent] = []
    locus_decisions: list[RepeatAwareLocusDecision] = []
    for locus in sorted(loci, key=lambda item: item.locus_id):
        left = set(locus.existing_left_pair_ids)
        right = set(locus.existing_right_pair_ids)
        for group in groups_by_locus[locus.locus_id]:
            if not decision_by_group[group.group_id].final_pass:
                continue
            if group.group_id not in shared_observations.aggregation_group_ids:
                continue
            destination = left if group.side is BreakpointSide.LEFT else right
            destination.add(group.representative_pair_id)
        left_ids = tuple(sorted(left))
        right_ids = tuple(sorted(right))
        passes = len(left_ids) >= 2 and len(right_ids) >= 2
        locus_decisions.append(
            RepeatAwareLocusDecision(locus.locus_id, left_ids, right_ids, passes)
        )
        if not passes:
            continue
        all_pair_ids = tuple(sorted(left | right))
        events.append(
            RepeatAwareEvent(
                sample_id=sample_context.sample_id,
                legacy_event_id=locus.legacy_event_id,
                family=locus.family,
                chrom=locus.chrom,
                interval0=locus.interval0,
                left_coord0=locus.left_coord0,
                right_coord0=locus.right_coord0,
                left_support_pairs=len(left_ids),
                right_support_pairs=len(right_ids),
                total_support_pairs=len(all_pair_ids),
                fragment_pair_ids=all_pair_ids,
                evidence_type="msbwt_repeat_aware_crossvalidated_te_divergence",
            )
        )
    return RepeatAwareResult(tuple(events), decisions, tuple(locus_decisions))
