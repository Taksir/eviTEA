"""Synthetic Bowtie2 host-state classification boundaries."""

import json
from pathlib import Path

from evitea.adjudication.host import (
    adjudicate_host_alignments,
    classify_host_mapping,
)
from evitea.backends.protocols import RawAlignment
from evitea.contracts import CensorState, HostMappingState, ReadMate

FIXTURE_DIR = Path(__file__).parents[1] / "fixtures" / "phase4" / "F-HOST-001"


def _alignment(row):
    return RawAlignment(
        query_name=row["query_name"],
        pair_id=row["pair_id"],
        mate=ReadMate(row["mate"]),
        flag=row["flag"],
        reference_name=row["reference_name"],
        position0=row["position0"],
        mapq=row["mapq"],
        cigar=row["cigar"],
        sequence=row["sequence"],
    )


def test_host_states_count_reported_placements_and_keep_mapq_boundary():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    expected = json.loads((FIXTURE_DIR / "expected.json").read_text(encoding="utf-8"))

    for case, rows in fixture["alignments"].items():
        state = classify_host_mapping(tuple(_alignment(row) for row in rows))
        assert state.value == expected[case]

    threshold = tuple(
        _alignment(row) for row in fixture["alignments"]["unique_at_threshold"]
    )
    low = tuple(
        _alignment(row) for row in fixture["alignments"]["one_low_mapq"]
    )
    multiple = tuple(
        _alignment(row) for row in fixture["alignments"]["two_reported_placements"]
    )
    unmapped = tuple(
        _alignment(row) for row in fixture["alignments"]["unmapped"]
    )
    assert classify_host_mapping(threshold) is HostMappingState.UNIQUE_HIGH
    assert classify_host_mapping(low) is HostMappingState.NONUNIQUE_OR_LOW
    assert classify_host_mapping(multiple) is HostMappingState.NONUNIQUE_OR_LOW
    assert classify_host_mapping(unmapped) is HostMappingState.UNMAPPED


def test_adjudication_preserves_resource_censor_and_each_mapped_placement():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    rows = tuple(
        _alignment(row)
        for row in fixture["alignments"]["two_reported_placements"]
    )

    result = adjudicate_host_alignments(
        "synthetic-sample",
        rows,
        censor_state=CensorState.HOST_RESOURCE_CENSORED,
    )

    assert result.state is HostMappingState.NONUNIQUE_OR_LOW
    assert result.censor_state is CensorState.HOST_RESOURCE_CENSORED
    assert len(result.placements) == 2
    assert [placement.pair_id for placement in result.placements] == [2, 2]
    assert [placement.mate_no for placement in result.placements] == [ReadMate.READ1] * 2
    assert [placement.placement_rank for placement in result.placements] == [1, 2]

