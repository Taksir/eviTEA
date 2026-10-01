"""Candidate-acquisition stages on small deterministic inputs."""

from __future__ import annotations

import random

import pytest

from evitea.errors import BackendCompatibilityError
from evitea.incumbent.anchor_census import build_anchor_catalogue, census_anchor_intervals
from evitea.incumbent.candidates import CandidateLocus
from evitea.incumbent.parallel import map_locus_batches
from evitea.incumbent.sequence import revcomp
from evitea.incumbent.sparse_panel import (
    CHANNEL_FALLBACK,
    CHANNEL_TERMINAL,
    route_sparse_panel,
    select_sparse_panel,
)

from incumbent_fakes import paired_sample


def _records():
    rng = random.Random(11)
    return {
        "FBte0000001": "".join(rng.choice("ACGT") for _ in range(400)),
        "FBte0000002": "".join(rng.choice("ACGT") for _ in range(60)) + "N" + "".join(rng.choice("ACGT") for _ in range(40)),
    }


def test_catalogue_layout_numbering_and_orientation():
    records = _records()
    catalogue = build_anchor_catalogue(records, {"FBte0000001": "roo"})
    keys = [(q.k, q.query_seq) for q in catalogue.queries]
    assert keys == sorted(keys)
    assert catalogue.queries[0].query_id == "Q0000001"
    long_sources = [s for s in catalogue.sources if s.te_id == "FBte0000001" and s.k == 21]
    regions = {s.region for s in long_sources}
    assert regions == {"5p_terminal", "3p_terminal", "internal"}
    assert {s.start for s in long_sources if s.region == "internal"} == {120, 170, 220, 259}
    assert 379 in {s.start for s in long_sources if s.region == "3p_terminal"}
    for source in long_sources[:4]:
        forward = records[source.te_id][source.start : source.start + 21]
        assert source.query_seq == (forward if source.orientation == "F" else revcomp(forward))
    short = [s for s in catalogue.sources if s.te_id == "FBte0000002"]
    assert {s.region for s in short} >= {"short_internal_fallback"}
    assert any(not s.queryable and s.query_id == "" for s in short)
    assert {s.family for s in short} == {"FBte0000002"}


def test_sparse_panel_prefers_narrowest_and_adds_fallback(monkeypatch):
    records = _records()
    catalogue = build_anchor_catalogue(records, {"FBte0000001": "roo", "FBte0000002": "jockey"})
    for number, query in enumerate(catalogue.queries):
        width = 3 if query.source_families == ("roo",) else 100
        query.interval_lo, query.interval_hi = 10 * number, 10 * number + width
    panel = select_sparse_panel(catalogue)
    assert panel and all(entry.width <= 64 or entry.channel_mask & CHANNEL_FALLBACK for entry in panel)
    fallback = [entry for entry in panel if entry.channel_mask & CHANNEL_FALLBACK]
    assert len(fallback) == 4 and all(entry.source_families == {"jockey"} for entry in fallback)
    assert any(entry.channel_mask & CHANNEL_TERMINAL for entry in panel if "roo" in entry.source_families)
    import evitea.incumbent.sparse_panel as sparse_panel

    monkeypatch.setattr(sparse_panel, "PANEL_ROW_HARD_CAP", sum(entry.width for entry in panel) - 1)
    with pytest.raises(BackendCompatibilityError):
        select_sparse_panel(catalogue)


def test_panel_routing_unions_mates_and_families():
    records = _records()
    seq = records["FBte0000001"]
    pairs = [(seq[0:40], "T" * 40), ("A" * 40, revcomp(seq[0:40])), ("C" * 40, "G" * 40)]
    index = paired_sample(pairs)
    catalogue = build_anchor_catalogue({"FBte0000001": seq[:200]}, {"FBte0000001": "roo"})
    census_anchor_intervals(catalogue, index)
    panel = select_sparse_panel(catalogue)
    rows = {row.pair_id: row for row in route_sparse_panel(panel, index)}
    assert set(rows) == {0, 1}
    assert rows[0].mate1_hit and not rows[0].mate2_hit and rows[0].families == ("roo",)
    assert rows[1].mate2_hit and not rows[1].mate1_hit and rows[1].channel_mask & CHANNEL_TERMINAL


def test_locus_batches_are_order_preserving_across_processes():
    loci = [CandidateLocus(f"L{i:07d}", "chr2L", "fam", 100.0 + i, 200.0 + i, (i,)) for i in range(7)]
    state = {"reads": {}, "reference": {"chr2L": "A" * 2000}, "family_kmers": {}, "te_index": None}
    serial = map_locus_batches("split_junction", loci, state, workers=1, batch_size=2)
    pooled = map_locus_batches("split_junction", loci, state, workers=2, batch_size=2)
    assert serial == pooled == ([], {"NO_CANDIDATE_FAMILY_KMER_SET": 7})
