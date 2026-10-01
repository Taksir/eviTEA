"""INTERNAL21 acquisition classification and fixed-proposal observability."""

from __future__ import annotations

from evitea.full_memory.strict_a import build_b_index
from evitea.internal21.acquisition import classify_internal_pairs, pair_probe_sequences
from evitea.internal21.observability import novel_unconfounded_loci, observe_increment

P1, P2, P3 = "A" * 20 + "C", "C" * 20 + "G", "G" * 20 + "T"


def test_classification_states_and_pair_probes():
    family = {P1: "F", P2: "F", P3: "G"}
    hits = {1: {1: {P1, P2}, 2: set()}, 2: {1: {P1}, 2: {P2}}, 3: {1: set(), 2: {P1, P3}}}
    meta, states = classify_internal_pairs(hits, family)
    assert set(meta) == {1} and meta[1]["te_mate"] == 1 and meta[1]["host_mate"] == 2
    assert states == {"ONE_MATE_FAMILY_RESOLVED": 1, "BOTH_MATES_TE_POSITIVE": 1, "FAMILY_CONFLICT": 1}
    assert len(pair_probe_sequences(meta[1])) == 4


def witness(pid, bp=1000, side="RIGHT", family="F"):
    return {"pair_id": pid, "fold": pid % 5, "family": family, "chrom": "chr2L", "breakpoint0": bp, "anchor_host_side": side}


def test_increment_requires_internal_a_and_is_classified():
    # held-out folds (0, 1): A from pairs 5 (fold 0) and 11 (fold 1); B two pairs per side in folds 0/1.
    b_rows = [{"pair_id": p, "fold": p % 5, "family": "F", "chrom": "chr2L", "start": s, "end": s + 100, "strand": st}
              for p, s, st in ((20, 500, "+"), (26, 600, "+"), (30, 1200, "-"), (36, 1300, "-"))]
    index, _ = build_b_index(b_rows)
    proposals = [{"sample": "S", "rotation": "H01", "candidate_id": "S_H01_000001", "family": "F", "chrom": "chr2L",
                  "candidate_lo": 900, "candidate_hi": 1100, "discovery_pair_ids": "2;3"}]
    old = [witness(5)]
    internal = [witness(11)]
    events, loci = observe_increment("S", proposals, old, internal, index, prior_events=[])
    assert [(e["rotation"], e["breakpoint0"], e["internal_A_pair_ids"]) for e in events] == [("H01", 1000, [11])]
    assert [x.classification for x in loci] == ["NOVEL_UNCONFOUNDED"]
    (locus,) = novel_unconfounded_loci(loci)
    assert locus.representative()[0].rotation == "H01"
    prior = [{"sample": "S", "family": "F", "chrom": "chr2L", "breakpoint0": 1400, "validated_anchorpe": 1}]
    _, loci = observe_increment("S", proposals, old, internal, index, prior_events=prior)
    assert loci[0].classification == "PRIOR_SAME_FAMILY_REPRESENTATION"
    # Without the internal witness nothing is incremental.
    assert observe_increment("S", proposals, old + internal, [], index, prior_events=[])[0] == []
