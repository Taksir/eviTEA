"""Physical-pair family decoding without row-level support inflation."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

from ..contracts import (
    FamilyResolutionState,
    ReadMate,
    SafeProbeRecordStatus,
    TEEligibilityState,
    TEFamilyEvidence,
    TEProbe,
    TEProbeHit,
)
from ..errors import InputContractError


def decode_pair_family(
    sample_id: str,
    pair_id: int,
    probe_hits: Iterable[TEProbeHit],
    codebook: Mapping[str, TEProbe],
    *,
    probe_memberships: Mapping[str, tuple[TEProbe, ...]] | None = None,
) -> TEFamilyEvidence:
    """Collapse observations to one physical-pair state.

    Duplicate occurrence rows are retained upstream but do not create new
    probe, block, record, family, mate, or pair support here.
    """

    if not sample_id:
        raise InputContractError("sample_id must be non-empty text")
    if isinstance(pair_id, bool) or not isinstance(pair_id, int) or pair_id < 0:
        raise InputContractError("pair_id must be a non-negative integer")
    rows = tuple(probe_hits)
    probes: list[TEProbe] = []
    mates: set[ReadMate] = set()
    for hit in rows:
        if hit.sample_id != sample_id or hit.pair_id != pair_id:
            raise InputContractError("probe hits mix sample or physical-pair identity")
        try:
            probe = codebook[hit.probe_id]
        except KeyError as exc:
            raise InputContractError(f"unknown safe probe: {hit.probe_id}") from exc
        probes.extend(
            (probe,)
            if probe_memberships is None
            else probe_memberships.get(hit.probe_id, (probe,))
        )
        mates.add(hit.mate_no)

    if not probes:
        return TEFamilyEvidence(
            sample_id=sample_id,
            pair_id=pair_id,
            family_state=FamilyResolutionState.NO_SAFE_EVIDENCE,
            eligibility_state=TEEligibilityState.NO_SAFE_EVIDENCE,
        )

    families = tuple(sorted({probe.family for probe in probes}))
    probe_ids = tuple(sorted({probe.probe_id for probe in probes}))
    block_ids = tuple(sorted({probe.independent_block_id for probe in probes}))
    record_ids = tuple(sorted({probe.te_record_id for probe in probes}))
    common = dict(
        sample_id=sample_id,
        pair_id=pair_id,
        families_observed=families,
        probe_ids=probe_ids,
        independent_block_ids=block_ids,
        contributing_te_record_ids=record_ids,
    )
    if len(families) != 1:
        return TEFamilyEvidence(
            **common,
            family_state=FamilyResolutionState.FAMILY_CONFLICT,
            eligibility_state=TEEligibilityState.FAMILY_CONFLICT,
        )
    family = families[0]
    if len(mates) != 1:
        return TEFamilyEvidence(
            **common,
            family_state=FamilyResolutionState.RESOLVED_SINGLETON,
            resolved_family=family,
            eligibility_state=TEEligibilityState.BOTH_MATES_TE_POSITIVE,
        )
    mate = next(iter(mates))

    blocks_by_multi_record: dict[str, set[str]] = {}
    for probe in probes:
        if probe.record_safe_status is SafeProbeRecordStatus.MULTI_INDEPENDENT_BLOCK:
            blocks_by_multi_record.setdefault(probe.te_record_id, set()).add(
                probe.independent_block_id
            )
    eligible = any(len(blocks) >= 2 for blocks in blocks_by_multi_record.values())
    if eligible:
        eligibility = TEEligibilityState.ELIGIBLE
    elif all(
        probe.record_safe_status is SafeProbeRecordStatus.UNSUPPORTED_TERMINAL_SAFE
        for probe in probes
    ):
        eligibility = TEEligibilityState.UNSUPPORTED_SAFE_TE_RECORD
    else:
        eligibility = TEEligibilityState.INSUFFICIENT_INDEPENDENT_BLOCKS
    return TEFamilyEvidence(
        **common,
        family_state=FamilyResolutionState.RESOLVED_SINGLETON,
        resolved_family=family,
        te_positive_mate=mate,
        eligibility_state=eligibility,
    )
