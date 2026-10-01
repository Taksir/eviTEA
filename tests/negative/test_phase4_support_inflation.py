"""Negative guard: repeated rows for one probe cannot manufacture a second block."""

import json
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


def test_four_occurrence_rows_for_one_block_remain_ineligible():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    row = next(item for item in fixture["codebook"] if item["probe_id"] == "a0")
    codebook = {
        "a0": TEProbe(
            probe_id="a0",
        probe_class=ProbeClass.TERMINAL_21MER,
            sequence=row["sequence"],
            family=row["family"],
            te_record_id=row["te_record_id"],
            source_span0=BaseSpan0(*row["source_span0"]),
            probe_source_strand=row["probe_source_strand"],
            independent_block_id=row["independent_block_id"],
            record_safe_status=SafeProbeRecordStatus[row["record_safe_status"]],
        )
    }
    repeated = [
        TEProbeHit(
            sample_id="synthetic-sample",
            pair_id=6,
            mate_no=ReadMate.READ1,
            probe_id="a0",
            query_span0=BaseSpan0(0, 21),
            hit_orientation=TEEventOrientation.FORWARD,
        )
        for _ in range(4)
    ]
    assert len(repeated) == 4

    evidence = decode_pair_family(
        "synthetic-sample", 6, repeated, codebook
    )

    assert evidence.family_state is FamilyResolutionState.RESOLVED_SINGLETON
    assert evidence.eligibility_state is TEEligibilityState.INSUFFICIENT_INDEPENDENT_BLOCKS
    assert evidence.probe_ids == ("a0",)
    assert evidence.independent_block_ids == ("TE_A:block0",)

