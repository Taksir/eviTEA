"""Synthetic exact-projection parity for forward/reverse occurrences."""

import json
from pathlib import Path

import pytest

from evitea.backends.protocols import RawReadHit
from evitea.contracts import (
    BaseSpan0,
    ProbeClass,
    ReadMate,
    SafeProbeRecordStatus,
    TEEventOrientation,
    TEProbe,
)
from evitea.errors import BackendCompatibilityError
from evitea.evidence.codebook import reverse_complement
from evitea.evidence.projection import SampleReadIndex, project_probes

FIXTURE_DIR = Path(__file__).parents[1] / "fixtures" / "phase4" / "F-TE-001"


class _FixtureBackend:
    def __init__(self, hits_by_query, count_overrides=None):
        self.hits_by_query = hits_by_query
        self.count_overrides = count_overrides or {}
        self.located = []
        self.counted = []

    def locate_reads(self, sequence):
        self.located.append(sequence)
        return tuple(self.hits_by_query.get(sequence, ()))

    def count_exact(self, sequence):
        self.counted.append(sequence)
        return self.count_overrides.get(sequence, len(self.hits_by_query.get(sequence, ())))

    def build(self, reads, output_dir):
        raise AssertionError("projection fixture must not build an index")


def _fixture():
    document = json.loads((FIXTURE_DIR / "projection.json").read_text(encoding="utf-8"))
    sequence = document["probe_sequence"]
    probe = TEProbe(
        probe_id="projection-probe",
        probe_class=ProbeClass.TERMINAL_21MER,
        sequence=sequence,
        family="FamilyA",
        te_record_id="TE_A",
        source_span0=BaseSpan0(0, 21),
        probe_source_strand="+",
        independent_block_id="TE_A:block0",
        record_safe_status=SafeProbeRecordStatus.MULTI_INDEPENDENT_BLOCK,
    )
    sample_sequences = {
        (row["pair_id"], ReadMate(row["mate"])): row["sequence"]
        for row in document["sample_sequences"]
    }
    sample_index = SampleReadIndex("synthetic-sample", sample_sequences)

    hits_by_query = {}
    for orientation, query in (
        (TEEventOrientation.FORWARD, sequence),
        (TEEventOrientation.REVERSE, reverse_complement(sequence)),
    ):
        expected = [
            row for row in document["expected_hits"]
            if row["hit_orientation"] == orientation.name
        ]
        hits_by_query[query] = tuple(
            RawReadHit(
                sequence=query,
                backend_row=ordinal0,
                backend_read_id=100 + ordinal0,
                pair_id=row["pair_id"],
                mate=ReadMate(row["mate"]),
                occurrence_ordinal0=ordinal0,
            )
            for ordinal0, row in enumerate(expected)
        )
    return document, probe, sample_index, hits_by_query


def test_projection_matches_exact_forward_reverse_occurrence_fixture():
    document, probe, sample_index, hits_by_query = _fixture()
    backend = _FixtureBackend(hits_by_query)

    hits = tuple(project_probes(sample_index, [probe], backend))

    observed = [
        {
            "pair_id": hit.pair_id,
            "mate": int(hit.mate_no),
            "query_span0": [hit.query_span0.start0, hit.query_span0.end0],
            "hit_orientation": hit.hit_orientation.name,
        }
        for hit in hits
    ]
    assert observed == document["expected_hits"]
    assert backend.located == backend.counted
    assert len(hits) == sum(len(rows) for rows in hits_by_query.values()) == 3


def test_projection_rejects_count_locate_multiplicity_disagreement():
    _document, probe, sample_index, hits_by_query = _fixture()
    forward_query = probe.sequence
    backend = _FixtureBackend(
        hits_by_query,
        count_overrides={forward_query: len(hits_by_query[forward_query]) + 1},
    )

    with pytest.raises(BackendCompatibilityError):
        tuple(project_probes(sample_index, [probe], backend))


def test_projection_rejects_backend_count_that_disagrees_with_read_occurrences():
    _document, probe, sample_index, hits_by_query = _fixture()
    forward_query = probe.sequence
    hits_by_query[forward_query] = hits_by_query[forward_query][:1]
    backend = _FixtureBackend(hits_by_query)

    with pytest.raises(BackendCompatibilityError):
        tuple(project_probes(sample_index, [probe], backend))

