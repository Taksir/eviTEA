"""Truth-blind synthetic selection of safe, independent 21-mer probes."""

import json
from pathlib import Path

from evitea.contracts import BaseSpan0, ProbeClass, SafeProbeRecordStatus
from evitea.evidence.codebook import (
    SafeProbeCandidate,
    build_safe_probes,
    classify_te_record_probe_status,
    validate_probe_independence,
)

FIXTURE_DIR = Path(__file__).parents[1] / "fixtures" / "phase4" / "F-TE-001"


def _candidates():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    return tuple(
        SafeProbeCandidate(
            candidate_id=row["candidate_id"],
            sequence=row["sequence"],
            family=row["family"],
            te_record_id=row["te_record_id"],
            source_span0=BaseSpan0(*row["source_span0"]),
            probe_source_strand=row["probe_source_strand"],
            independent_block_id=row["independent_block_id"],
            compatible_families=tuple(row["compatible_families"]),
            reference_occurrences=row["reference_occurrences"],
        )
        for row in fixture["candidates"]
    )


def test_codebook_keeps_only_family_unique_reference_absent_independent_probes():
    probes = build_safe_probes(_candidates(), ProbeClass.TERMINAL_21MER)
    by_id = {probe.probe_id: probe for probe in probes}

    assert set(by_id) == {"a0", "a1", "c0"}
    assert by_id["a0"].record_safe_status is SafeProbeRecordStatus.MULTI_INDEPENDENT_BLOCK
    assert by_id["a1"].record_safe_status is SafeProbeRecordStatus.MULTI_INDEPENDENT_BLOCK
    assert by_id["c0"].record_safe_status is SafeProbeRecordStatus.SINGLE_INDEPENDENT_BLOCK
    validate_probe_independence(probes)


def test_empty_safe_probe_set_has_explicit_unsupported_status():
    assert (
        classify_te_record_probe_status(())
        is SafeProbeRecordStatus.UNSUPPORTED_TERMINAL_SAFE
    )

