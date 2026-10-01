"""Synthetic end-to-end safe-probe observations through physical-pair decoding."""

import json
from collections import defaultdict
from pathlib import Path

from evitea.contracts import (
    BaseSpan0,
    FamilyResolutionState,
    ProbeClass,
    ReadMate,
    SafeProbeRecordStatus,
    TEEligibilityState,
    TEEventOrientation,
    TEProbe,
    TEProbeHit,
)
from evitea.evidence.fragments import decode_pair_family

FIXTURE_DIR = Path(__file__).parents[1] / "fixtures" / "phase4" / "F-TE-001"


def _load():
    input_doc = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    expected = json.loads((FIXTURE_DIR / "expected.json").read_text(encoding="utf-8"))
    codebook = {
        row["probe_id"]: TEProbe(
            probe_id=row["probe_id"],
        probe_class=ProbeClass.TERMINAL_21MER,
            sequence=row["sequence"],
            family=row["family"],
            te_record_id=row["te_record_id"],
            source_span0=BaseSpan0(*row["source_span0"]),
            probe_source_strand=row["probe_source_strand"],
            independent_block_id=row["independent_block_id"],
            record_safe_status=SafeProbeRecordStatus[row["record_safe_status"]],
        )
        for row in input_doc["codebook"]
    }
    hits_by_pair = defaultdict(list)
    for row in input_doc["hits"]:
        hits_by_pair[row["pair_id"]].append(
            TEProbeHit(
                sample_id="synthetic-sample",
                pair_id=row["pair_id"],
                mate_no=ReadMate(row["mate_no"]),
                probe_id=row["probe_id"],
                query_span0=BaseSpan0(*row["query_span0"]),
                hit_orientation=TEEventOrientation[row["hit_orientation"]],
            )
        )
    return input_doc, expected, codebook, hits_by_pair


def test_pairs_keep_conflict_both_mates_eligibility_and_duplicate_multiplicity():
    fixture, expected, codebook, hits_by_pair = _load()
    assert fixture["pair_ids"] == list(range(7))

    decoded = {
        pair_id: decode_pair_family(
            "synthetic-sample", pair_id, hits_by_pair[pair_id], codebook
        )
        for pair_id in fixture["pair_ids"]
    }

    for pair_id, evidence in decoded.items():
        target = expected[str(pair_id)]
        if "family_state" in target:
            assert evidence.family_state.value == target["family_state"]
        assert evidence.eligibility_state.value == target["eligibility_state"]
        if "resolved_family" in target:
            assert evidence.resolved_family == target["resolved_family"]
        if "te_positive_mate" in target:
            actual_mate = (
                None
                if evidence.te_positive_mate is None
                else int(evidence.te_positive_mate)
            )
            assert actual_mate == target["te_positive_mate"]
        if "families_observed" in target:
            assert list(evidence.families_observed) == target["families_observed"]

    assert decoded[0].family_state is FamilyResolutionState.RESOLVED_SINGLETON
    assert decoded[0].eligibility_state is TEEligibilityState.ELIGIBLE
    assert len(hits_by_pair[0]) == 4
    assert decoded[0].probe_ids == ("a0", "a1")
    assert decoded[0].independent_block_ids == ("TE_A:block0", "TE_A:block1")
    assert decoded[3].family_state is FamilyResolutionState.FAMILY_CONFLICT
    assert decoded[2].eligibility_state is TEEligibilityState.BOTH_MATES_TE_POSITIVE

