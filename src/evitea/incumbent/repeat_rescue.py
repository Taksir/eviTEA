"""Retain a group's call after credit for recurrent TE-copy edits, and reaggregate loci.

Groups rejected by confirmation *only* because the unique-best candidate
family exceeded ``35 * NM <= 3 * length`` are retained when their edits are
explained by exact TE-copy variants recurring in >= 2 *other* physical groups
at the same locus, family and TE record (leave-one-group-out).  Repeat-shift
equivalent indels need a linked recurrent substitution anchor.  Passing groups
are reaggregated with the 2-LEFT + 2-RIGHT physical-pair rule.

Reaggregation pairs each side's linked group hashes (in linking order) with
that side's representative pair IDs (ascending) when admitting passing groups.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..constants import repeat_aware_edit_rate_passes
from .events import IncumbentCall, IncumbentCallClass
from .junction_groups import JunctionLocus
from .repeat_aware import RepeatAwareVariant, alignment_variant_ops, indel_haplotype_admissibility
from .te_certification import GroupCertification, TEAlignment
from .te_library import normalize_family

MIN_MODEL_OTHER_GROUPS = 2
ELIGIBLE_STATUS = "REJECT_UNIQUE_BEST_BUT_RATE_ABOVE_3_PER_35"


@dataclass(frozen=True, slots=True)
class DivergenceProfile:
    group_hash: str
    locus_id: str
    family: str
    record_id: str
    nm: int
    operations: tuple[RepeatAwareVariant, ...]
    valid: bool
    representative_pair_id: int


def choose_candidate_alignment(
    alignments: Sequence[TEAlignment], candidate_family: str, record_families: Mapping[str, str]
) -> TEAlignment | None:
    family = normalize_family(candidate_family)
    rows = [row for row in alignments if normalize_family(record_families.get(row.record_id, "")) == family]
    if not rows:
        return None
    return min(rows, key=lambda row: (row.nm, row.record_id, int(row.reverse), row.pos1, row.cigar, row.sequence))


@dataclass(frozen=True, slots=True)
class DivergenceDecision:
    group_hash: str
    crossvalidated_pass: bool
    status: str
    adjusted_nm: int | None


def crossvalidate(
    certifications: Sequence[GroupCertification],
    te_records: Mapping[str, str],
    record_families: Mapping[str, str],
) -> dict[str, DivergenceDecision]:
    profiles: dict[str, DivergenceProfile] = {}
    for row in certifications:
        group = row.group
        family = normalize_family(group.candidate_family)
        alignment = choose_candidate_alignment(row.alignments, family, record_families)
        if alignment is None or alignment.record_id not in te_records:
            continue
        try:
            operations = alignment_variant_ops(
                alignment.sequence, te_records[alignment.record_id], alignment.pos1 - 1, alignment.cigar
            )
        except Exception:
            continue
        profiles[group.group_hash] = DivergenceProfile(
            group.group_hash, group.locus_id, family, alignment.record_id, alignment.nm, operations,
            sum(op.cost for op in operations) == alignment.nm, group.representative_pair_id,
        )
    pools: dict[tuple[str, str, str], list[DivergenceProfile]] = defaultdict(list)
    for profile in profiles.values():
        if profile.valid and profile.representative_pair_id >= 0:
            pools[(profile.locus_id, profile.family, profile.record_id)].append(profile)

    decisions: dict[str, DivergenceDecision] = {}
    for row in certifications:
        group_hash = row.group.group_hash
        eligible = not row.passes and row.status == ELIGIBLE_STATUS
        if not eligible:
            decisions[group_hash] = DivergenceDecision(group_hash, False, "NOT_ELIGIBLE", row.best_nm)
            continue
        profile = profiles.get(group_hash)
        if profile is None or not profile.valid or profile.representative_pair_id < 0:
            decisions[group_hash] = DivergenceDecision(group_hash, False, "NO_VALID_ALIGNMENT_PROFILE", row.best_nm)
            continue
        pool = [
            other for other in pools[(profile.locus_id, profile.family, profile.record_id)]
            if other.group_hash != group_hash and other.representative_pair_id != profile.representative_pair_id
        ]
        site_alleles: dict[tuple, dict[tuple, set[int]]] = defaultdict(lambda: defaultdict(set))
        for other in pool:
            for op in other.operations:
                site_alleles[op.site_key][op.exact_key].add(other.representative_pair_id)
        admissible: dict[tuple, set[int]] = {}
        for alleles in site_alleles.values():
            ranked = sorted(alleles.items(), key=lambda item: (-len(item[1]), item[0]))
            if not ranked or len(ranked[0][1]) < MIN_MODEL_OTHER_GROUPS:
                continue
            if len(ranked) > 1 and len(ranked[1][1]) == len(ranked[0][1]):
                continue
            admissible[ranked[0][0]] = ranked[0][1]
        held_out = {op.exact_key: op for op in profile.operations}
        recurrent = [key for key in held_out if key in admissible]
        recurrent_ops = [held_out[key] for key in recurrent]
        reference = te_records[profile.record_id]
        credited = [key for key in recurrent
                    if indel_haplotype_admissibility(reference, held_out[key], recurrent_ops)[0]]
        adjusted = max(0, profile.nm - sum(held_out[key].cost for key in credited))
        passes = bool(credited) and repeat_aware_edit_rate_passes(adjusted, row.segment_length)
        if passes:
            status = "PASS_REPEAT_AWARE_COHERENT_DIVERGENCE"
        elif not recurrent:
            status = "NO_RECURRENT_VARIANT_IN_HELDOUT"
        elif not credited:
            status = "ALL_RECURRENT_CREDIT_REMOVED_BY_INDEL_ANCHOR_GATE"
        else:
            status = "ADJUSTED_RATE_STILL_FAILS"
        decisions[group_hash] = DivergenceDecision(group_hash, passes, status, adjusted)
    return decisions


@dataclass(frozen=True, slots=True)
class ReaggregatedLocus:
    locus_id: str
    left_ids: tuple[int, ...]
    right_ids: tuple[int, ...]
    passes: bool


def reaggregate(ledger: Sequence[JunctionLocus], group_pass: Mapping[str, bool]) -> tuple[ReaggregatedLocus, ...]:
    out: list[ReaggregatedLocus] = []
    for row in ledger:
        left = set(row.existing_left_ids) | set(row.strict_left_ids)
        right = set(row.existing_right_ids) | set(row.strict_right_ids)
        for hashes, representatives, destination in (
            (row.left_group_hashes, row.left_representatives, left),
            (row.right_group_hashes, row.right_representatives, right),
        ):
            if len(hashes) != len(representatives):
                raise ValueError(f"parallel group/representative lists differ for {row.locus_id}")
            for group_hash, pair_id in zip(hashes, representatives):
                if group_pass.get(group_hash, False):
                    destination.add(pair_id)
        out.append(ReaggregatedLocus(row.locus_id, tuple(sorted(left)), tuple(sorted(right)),
                                     len(left) >= 2 and len(right) >= 2))
    return tuple(out)


def rescue_calls(
    sample_id: str, ledger: Sequence[JunctionLocus], aggregated: Sequence[ReaggregatedLocus]
) -> tuple[IncumbentCall, ...]:
    by_locus = {row.locus_id: row for row in aggregated}
    calls: list[IncumbentCall] = []
    for row in ledger:
        state = by_locus[row.locus_id]
        if not state.passes:
            continue
        start = int(math.floor(float(row.candidate_start)))
        end = int(math.ceil(float(row.candidate_end))) + 1
        if end <= start:
            end = start + 1
        union = tuple(sorted(set(state.left_ids) | set(state.right_ids)))
        calls.append(IncumbentCall(
            sample_id, f"V2R9E_REPEATDIV_{len(calls) + 1:06d}", IncumbentCallClass.REPEAT_AWARE_RESCUE,
            row.chrom, start, end, row.family, None, None, None, len(state.left_ids), len(state.right_ids),
            "msbwt_r9e_repeat_aware_crossvalidated_te_divergence", float(len(union)), union,
        ))
    return tuple(calls)
