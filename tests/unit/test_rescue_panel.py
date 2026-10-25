"""The widest panel (1,024 rows) is always the one kept, because the pair-aware reselection is built for it."""

from __future__ import annotations

import random

import pytest

from evitea.errors import BackendCompatibilityError
from evitea.incumbent import coverage_panel
from evitea.incumbent.anchor_census import build_anchor_catalogue, census_anchor_intervals
from evitea.incumbent.pair_aware_panel import FROZEN_RESCUE_CAP, select_pair_aware_panel
from evitea.incumbent.sequence import revcomp
from evitea.incumbent.sparse_panel import select_sparse_panel

from incumbent_fakes import paired_sample


def _rare_word_inputs():
    """Two TE records whose k-mers each occur a few times in the reads: every k-mer is rare."""

    rng = random.Random(9)
    records = {f"TE{i}": "".join(rng.choice("ACGT") for _ in range(400)) for i in (1, 2)}
    families = {"TE1": "famA", "TE2": "famB"}
    pairs = []
    for sequence in records.values():
        for start in range(0, 280, 40):
            host = "".join(rng.choice("ACGT") for _ in range(60))
            pairs.append((sequence[start:start + 60], revcomp(host)))
    index = paired_sample(pairs)
    catalogue = build_anchor_catalogue(records, families)
    census_anchor_intervals(catalogue, index)
    panel = select_sparse_panel(catalogue)
    short = coverage_panel.short_exact_census(catalogue.records, catalogue.family_by_te(), index)
    h1 = coverage_panel.hard_gap_h1_loci(
        catalogue.records, catalogue.family_by_te(), coverage_panel.baseline_anchor_loci(catalogue, panel), short
    )
    return catalogue, panel, short, h1, index


def test_rare_words_keep_the_1024_row_panel(monkeypatch: pytest.MonkeyPatch) -> None:
    catalogue, panel, short, h1, index = _rare_word_inputs()
    scores = {}
    original = coverage_panel._optimize_panel

    def every_cap(*args, **kwargs):
        for cap in coverage_panel.PANEL_CAPS:  # the three caps tie here, so the 256-row panel would win on score alone
            scores[cap] = original(cap, *args[1:], **kwargs)[1]
        return original(*args, **kwargs)

    monkeypatch.setattr(coverage_panel, "_optimize_panel", every_cap)
    rescue = coverage_panel.optimize_rescue_panel(catalogue, panel, short, h1, index)
    assert scores[256] == scores[512] == scores[1024]
    assert max(coverage_panel.PANEL_CAPS, key=lambda cap: scores[cap]) == 256
    assert rescue.cap == FROZEN_RESCUE_CAP == coverage_panel.RESCUE_PANEL_CAP == 1024
    select_pair_aware_panel(catalogue, panel, short, rescue, index)  # accepted


def test_a_1024_row_panel_over_the_row_limit_stops_with_a_clear_error(monkeypatch: pytest.MonkeyPatch) -> None:
    catalogue, panel, short, h1, index = _rare_word_inputs()
    monkeypatch.setattr(coverage_panel, "MATERIALIZE_ROW_HARD_CAP", 0)
    with pytest.raises(BackendCompatibilityError, match="1024-row coverage panel needs .* materialization limit"):
        coverage_panel.optimize_rescue_panel(catalogue, panel, short, h1, index)
