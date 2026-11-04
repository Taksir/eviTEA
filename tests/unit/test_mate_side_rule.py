"""``--mate-side-rule strand`` and ``--junction-side-rule split``: the split points, and that the defaults are unchanged."""

from __future__ import annotations

import math
from types import SimpleNamespace

import pytest

from evitea.cli import build_parser
from evitea.config import DepthFloorConfig
from evitea.errors import InputContractError
from evitea.incumbent import mate_rescue as MR
from evitea.incumbent.mate_rescue import classify_side, junction_split_point, strand_split_point

INF = math.inf


def test_one_strand_places_use_infinite_points():
    assert strand_split_point([("FWD", 100.0), ("FWD", 300.0)], 200.0) == INF
    assert strand_split_point([("RC", 100.0), ("RC", 300.0)], 200.0) == -INF
    # every FWD anchor is LEFT, every RC anchor is RIGHT, whatever the midpoint of the place
    assert all(classify_side("FWD", m, INF) == "LEFT" for m in (100.0, 300.0))
    assert all(classify_side("RC", m, -INF) == "RIGHT" for m in (100.0, 300.0))
    # the midpoint rule would have dropped the FWD anchor right of the middle
    assert classify_side("FWD", 300.0, 200.0) is None


def test_both_strands_with_a_clean_gap_split_in_the_middle_of_the_gap():
    x = strand_split_point([("FWD", 100.0), ("FWD", 120.0), ("RC", 300.0), ("RC", 340.0)], 150.0)
    assert x == 210.0
    assert [classify_side("FWD", m, x) for m in (100.0, 120.0)] == ["LEFT", "LEFT"]
    assert [classify_side("RC", m, x) for m in (300.0, 340.0)] == ["RIGHT", "RIGHT"]


def test_both_strands_overlapping_keep_the_most_consistent_pairs():
    # FWD at 100 and 300, RC at 200: x in (100, 200) keeps FWD 100 and RC 200; x > 300 keeps both FWD (tie at two);
    # the bounded gap wins the tie and x is its middle
    x = strand_split_point([("FWD", 100.0), ("FWD", 300.0), ("RC", 200.0)], 200.0)
    assert x == 150.0
    # three FWD to the left of one RC beats the single RC: the clean gap
    assert strand_split_point([("FWD", 100.0), ("FWD", 110.0), ("FWD", 120.0), ("RC", 90.0)], 100.0) == INF
    # an inverted pair (RC left of FWD), one each: both outer gaps score one; the smaller x comes first
    assert strand_split_point([("RC", 100.0), ("FWD", 500.0)], 300.0) == -INF


def test_ties_take_the_widest_gap_then_the_one_closest_to_the_locus_midpoint():
    # two bounded gaps with the same best score (2): (100,110) wide 10 and (200,300) wide 100
    anchors = [("FWD", 100.0), ("RC", 110.0), ("FWD", 200.0), ("RC", 300.0)]
    # scores: x in (100,110): FWD100 + RC110/300 = 3? check by brute force below
    x = strand_split_point(anchors, 0.0)
    best = max(sum(1 for v, m in anchors if v == "FWD" and m < t) + sum(1 for v, m in anchors if v == "RC" and m > t)
               for t in (105.0, 150.0, 250.0, 99.0, 301.0))
    score = sum(1 for v, m in anchors if v == "FWD" and m < x) + sum(1 for v, m in anchors if v == "RC" and m > x)
    assert score == best
    # equal width and score: the gap closer to the locus midpoint
    anchors = [("FWD", 100.0), ("RC", 200.0), ("RC", 300.0), ("FWD", 400.0), ("RC", 500.0)]
    far = strand_split_point(anchors, 1000.0)
    near = strand_split_point(anchors, 0.0)
    assert math.isfinite(far) and math.isfinite(near)
    assert abs(near - 0.0) <= abs(far - 0.0) and abs(far - 1000.0) <= abs(near - 1000.0)
    assert strand_split_point([], 42.0) == 42.0


def test_split_point_is_a_maximum_by_brute_force():
    import random

    rng = random.Random(3)
    for _ in range(200):
        anchors = [(rng.choice(["FWD", "RC"]), float(rng.randrange(0, 40))) for _ in range(rng.randrange(1, 9))]
        locus = float(rng.randrange(0, 40))
        x = strand_split_point(anchors, locus)

        def score(t):
            return sum(1 for v, m in anchors if v == "FWD" and m < t) + sum(1 for v, m in anchors if v == "RC" and m > t)

        probes = [-INF, INF] + [t + 0.5 for t in range(-1, 41)]
        assert score(x) == max(score(t) for t in probes), (anchors, x)


def test_junction_point_is_the_present_sides_edge_at_a_one_strand_place():
    assert junction_split_point([("FWD", 100.0), ("FWD", 300.0)], 200.0) == 300.0     # LEFT present: rightmost FWD host
    assert junction_split_point([("RC", 100.0), ("RC", 300.0)], 200.0) == 100.0       # RIGHT present: leftmost RC host
    both = [("FWD", 100.0), ("RC", 300.0)]
    assert junction_split_point(both, 200.0) == strand_split_point(both, 200.0) == 200.0
    # a missing-side host (RC, right of the FWD edge) is RIGHT against the edge point, as the junction test reads it
    edge = junction_split_point([("FWD", 100.0), ("FWD", 300.0)], 200.0)
    assert 350.0 > edge and 250.0 < edge


def _locus(midpoint=200.0):
    return SimpleNamespace(locus_id="L1", chrom="chrS", family="F", start=190.0, end=210.0, midpoint=midpoint, pair_ids=(1, 2, 3))


def test_strand_rule_gives_a_side_midpoint_would_drop_and_counts_it(monkeypatch):
    anchors = {1: ("FWD", 120.0), 2: ("FWD", 250.0), 3: ("FWD", 280.0)}      # one strand; two lie right of the middle (200)
    reads = {pid: ("A" * 50, "C" * (49 + pid)) for pid in anchors}      # a different host read per pair
    reference = {"chrS": "ACGT" * 500}

    def fake_anchor(read, scope, scope_start, start, end):
        pid = len(read) - 49
        # the host mate is mate 2 of every pair; mate 1 finds no anchor
        return MR.HostAnchor(anchors[pid][0], int(anchors[pid][1]) - 25, int(anchors[pid][1]) + 25, anchors[pid][1], 0, 0.0), "OK"

    monkeypatch.setattr(MR, "best_local_host_anchor", lambda read, *a: fake_anchor(read, *a) if read.startswith("C") else (None, "HOST_MATE_NO_LOCAL_MATCH"))
    monkeypatch.setattr(MR, "best_te_mate_window", lambda *a: (SimpleNamespace(seed_hits=5, seed_fraction=0.5, te_edits=0), "OK"))
    common = dict(loci=[_locus()], reads=reads, reference=reference, family_kmers={"F": {"AAAA"}}, te_index=None)
    old, old_counts = MR.mate_rescue_evidence(**common)
    new, new_counts = MR.mate_rescue_evidence(**common, side_rule="strand")
    assert [row.pair_id for row in old] == [1] and old_counts["PAIR_GEOMETRY_INCONSISTENT"] == 2
    assert sorted(row.pair_id for row in new if isinstance(row, MR.MateEvidence)) == [1, 2, 3]
    assert {row.side for row in new if isinstance(row, MR.MateEvidence)} == {"LEFT"}
    assert new_counts["STRAND_RULE_ANCHORS_ADDED_ONE_STRAND_LOCI"] == 2 and new_counts["STRAND_RULE_LOCI_ONE_STRAND_LOCI"] == 1
    assert "PAIR_GEOMETRY_INCONSISTENT" not in new_counts
    # the points alone (junction rule) leave the sides and the counts of the midpoint rule
    pts, pts_counts = MR.mate_rescue_evidence(**common, record_split=True)
    assert [row.pair_id for row in pts if isinstance(row, MR.MateEvidence)] == [1] and pts_counts["PAIR_GEOMETRY_INCONSISTENT"] == 2
    [split] = [row for row in pts if isinstance(row, MR.LocusSplit)]
    assert split.junction_point == 280.0 and split.strand_point == INF
    # finalize takes the points out of the evidence and keeps them
    result = MR.finalize_mate_rescue("S", [_locus_as_candidate()], [], pts, pts_counts)
    assert isinstance(result, MR.MateRescueSplitResult) and result.splits["L1"].junction_point == 280.0


def _locus_as_candidate():
    from evitea.incumbent.candidates import CandidateLocus
    import dataclasses

    names = [f.name for f in dataclasses.fields(CandidateLocus)]
    values = {"locus_id": "L1", "chrom": "chrS", "family": "F", "start": 190.0, "end": 210.0, "pair_ids": (1, 2, 3)}
    return SimpleNamespace(**{name: values.get(name, 0) for name in names}, midpoint=200.0)


def test_default_calls_the_search_as_before_and_returns_the_plain_result(monkeypatch):
    monkeypatch.setattr(MR, "best_local_host_anchor", lambda *a: (None, "HOST_MATE_NO_LOCAL_MATCH"))
    ev, counts = MR.mate_rescue_evidence(loci=[_locus()], reads={1: ("A", "C"), 2: ("A", "C"), 3: ("A", "C")}, reference={"chrS": "ACGT" * 500},
                                         family_kmers={"F": {"AAAA"}}, te_index=None)
    assert ev == [] and counts == {"HOST_MATE_NO_LOCAL_MATCH": 6}
    result = MR.finalize_mate_rescue("S", [], [], [], {})
    assert type(result) is MR.MateRescueResult
    with pytest.raises(ValueError):
        MR.mate_rescue_evidence(loci=[], reads={}, reference={}, family_kmers={}, te_index=None, side_rule="x")


def test_strict_pair_and_junction_tests_use_the_points_only_when_given():
    assert MR._side_point(None, "L1", 200.0) == 200.0
    assert MR._side_point({"L1": 300.0}, "L1", 200.0) == 300.0
    assert MR._side_point({"L2": 300.0}, "L1", 200.0) == 200.0


def test_config_and_command_line():
    assert DepthFloorConfig().is_default
    assert not DepthFloorConfig(mate_side_rule="strand").is_default and not DepthFloorConfig(junction_side_rule="split").is_default
    with pytest.raises(InputContractError):
        DepthFloorConfig(mate_side_rule="x")
    with pytest.raises(InputContractError):
        DepthFloorConfig(junction_side_rule="strand")
    call = build_parser().parse_args(["call"])
    assert call.mate_side_rule is None and call.junction_side_rule is None
    call = build_parser().parse_args(["call", "--mate-side-rule", "strand", "--junction-side-rule", "split"])
    assert (call.mate_side_rule, call.junction_side_rule) == ("strand", "split")


def test_stage_fingerprints_keep_the_early_checkpoints_reusable(tmp_path):
    from evitea.incumbent.engine import StageStore

    base = {"schema": "x"}
    default = StageStore(tmp_path, base)
    assert default.fingerprint_of("mate_rescue") == base == default.fingerprint_of("candidate_acquisition")
    strand = StageStore(tmp_path, base, {"mate_side_rule": "strand"}, ("mate_rescue", "strict_pairs"))
    assert strand.fingerprint_of("candidate_acquisition") == base and strand.fingerprint_of("split_junction") == base
    assert strand.fingerprint_of("mate_rescue") == {**base, "mate_side_rule": "strand"}


def test_run_engine_gives_the_rules_to_the_library_layer_and_records_them(tmp_path, monkeypatch):
    import evitea.engine as engine_module
    from evitea.engine import EngineInputs
    from evitea.incumbent.acquisition import CandidateAcquisition
    from evitea.incumbent.engine import IncumbentResult
    from evitea.sample_geometry import FragmentSizes
    from test_default_engine import _inputs as engine_inputs
    from test_depth_floors import _results

    seen = {}
    split, mates, groups = _results(pairs=[])

    def incumbent(inputs, sample, work, progress=None):
        seen["inputs"] = inputs
        return IncumbentResult((), CandidateAcquisition((), (), {}), split, mates, None, None, groups, (), {"incumbent_calls": 0})

    def extension(anchor_inputs, internal_codebook, sample, incumbent_bytes, labels, work, progress=None):
        from evitea.shadow_layers import ShadowCallset

        layer = SimpleNamespace(layer="INTERNAL21", callset=ShadowCallset.from_bytes(incumbent_bytes), unresolved=())
        return SimpleNamespace(layers=(layer,), anchor_pe=SimpleNamespace(counts={}), full_memory=SimpleNamespace(counts={}),
                               internal21=SimpleNamespace(counts={}), labels=labels)

    monkeypatch.setattr(engine_module, "MsbwtSample", lambda path: ("msbwt", path))
    monkeypatch.setattr(engine_module, "run_incumbent_caller", incumbent)
    monkeypatch.setattr(engine_module, "run_extension", extension)
    monkeypatch.setattr(engine_module, "estimate_fragment_sizes", lambda *a, **k: FragmentSizes(300.0, 420, 5000, "estimate"))
    base = engine_inputs(tmp_path)
    fields = {f: getattr(base, f) for f in base.__slots__}

    result = engine_module.run_engine(EngineInputs(**fields), tmp_path / "w0")             # defaults: nothing is passed or recorded
    assert (seen["inputs"].mate_side_rule, seen["inputs"].junction_side_rule) == ("midpoint", "midpoint")
    assert "depth_floors" not in result.counts

    rules = DepthFloorConfig(mate_side_rule="strand", junction_side_rule="split")
    result = engine_module.run_engine(EngineInputs(**{**fields, "depth_floors": rules}), tmp_path / "w1")
    assert (seen["inputs"].mate_side_rule, seen["inputs"].junction_side_rule) == ("strand", "split")
    assert result.counts["depth_floors"]["side_rules"] == {"mate_side_rule": "strand", "junction_side_rule": "split"}
    assert "weak_site_pass" not in result.counts["depth_floors"] and "candidate_floor" not in result.counts["depth_floors"]


def test_side_test_counts_and_manifest_keys_appear_only_for_an_opt_in_rule():
    from evitea.incumbent.engine import SIDE_TEST_COUNTS, _mate_side_counts, side_test_counts

    mates = SimpleNamespace(counts={"PAIR_GEOMETRY_INCONSISTENT": 7, "STRAND_RULE_ANCHORS_ADDED_ONE_STRAND_LOCI": 4,
                                    "STRAND_RULE_LOCI_ONE_STRAND_LOCI": 2})
    strict = SimpleNamespace(counts={"PAIR_GEOMETRY_INCONSISTENT_LOCAL": 3, "REJECT_GLOBAL_HOST_GEOMETRY_INCONSISTENT": 1})
    groups = SimpleNamespace(counts={"REJECT_JUNCTION_SIDE_UNRESOLVED": 5, "REJECT_JUNCTION_NOT_MISSING_SIDE": 2})
    counts = side_test_counts(mates, strict, groups)
    assert set(counts) == {name for name, _, _ in SIDE_TEST_COUNTS}
    assert counts["side_test_rejects_main_loop"] == 7 and counts["side_test_rejects_local_geometry"] == 3
    assert counts["side_test_rejects_global_geometry"] == 1 and counts["side_test_rejects_junction_unresolved"] == 5
    assert counts["side_test_rejects_junction_not_missing_side"] == 2 and counts["side_test_rejects_local_not_missing_side"] == 0
    out = _mate_side_counts(mates, strict, groups)
    assert out["mate_side_strand_rule_anchors_added_one_strand_loci"] == 4 and out["mate_side_strand_rule_loci_one_strand_loci"] == 2


def test_a_checkpoint_of_an_early_stage_is_reused_and_a_later_one_is_refused(tmp_path):
    from evitea.errors import ResumeStateError
    from evitea.incumbent.engine import StageStore

    base = {"schema": "x"}
    later = ("mate_rescue", "host_seed_routing", "strict_pairs", "junction_groups")
    default = StageStore(tmp_path, base, None, later)
    assert default.get("candidate_acquisition", lambda: "A") == "A" and default.get("mate_rescue", lambda: "M") == "M"
    strand = StageStore(tmp_path, base, {"mate_side_rule": "strand"}, later)
    assert strand.get("candidate_acquisition", lambda: "never") == "A"                      # reused
    with pytest.raises(ResumeStateError):
        strand.get("mate_rescue", lambda: "never")                                        # computed under another rule
    # the checkpoint files of a default run are those of a store without the later identity
    plain = StageStore(tmp_path / "p", base)
    plain.get("mate_rescue", lambda: "M")
    again = StageStore(tmp_path / "p", base, None, later)
    assert again.get("mate_rescue", lambda: "never") == "M"
    assert (tmp_path / "mate_rescue.complete.json").read_bytes() == (tmp_path / "p" / "mate_rescue.complete.json").read_bytes()
