"""Opt-in depth-aware floors of the library layer: weak-site pass, candidate place minimums, their config, tiers and identity.

Every option is off by default; the default must change no output (the byte-for-byte checks of the library-free layer are in
``tools/library_free_tests/test_depth_floors.py``; the manifest and identity checks of the library layer are here).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from evitea import callset, tiers
from evitea import depth_floors as DF
from evitea.config import DepthFloorConfig, EviTEAConfig
from evitea.engine import EngineInputs
from evitea.errors import InputContractError
from evitea.incumbent import acquisition as ACQ
from evitea.incumbent import weak_sites as W
from evitea.incumbent.engine import IncumbentInputs, _fingerprint
from evitea.incumbent.events import IncumbentCall, IncumbentCallClass
from evitea.incumbent.junction_groups import JunctionLocus
from evitea.incumbent.mate_rescue import MateRescueLocus, MateRescueResult
from evitea.incumbent.split_junction import SplitJunctionResult


def _event(kind, chrom, start, end, family, ids):
    exact = kind == "two_sided"
    return IncumbentCall("s", "E", IncumbentCallClass.SPLIT_JUNCTION, chrom, start, end, family,
                         start if exact else (start if kind == "one_sided_left" else -1),
                         end if exact else (end if kind == "one_sided_right" else -1),
                         0, len(ids) if kind != "one_sided_right" else 0,
                         len(ids) if kind != "one_sided_left" else 0, kind, 1.0, tuple(sorted(ids)))


def _pair(chrom, family, a, b, left=(), right=()):
    return MateRescueLocus("L", chrom, family, float(a), float(b), float(a + b) / 2, tuple(left), tuple(right), False)


def _strict(chrom, family, a, b, left=(), right=()):
    return JunctionLocus("L", chrom, family, float(a), float(b), (), (), tuple(left), tuple(right), (), (), (), ())


def _results(events=(), pairs=(), stricts=()):
    return (SplitJunctionResult(tuple(events), (), (), {}), MateRescueResult((), tuple(pairs), (), {}),
            SimpleNamespace(ledger=tuple(stricts)))


def _run(events=(), pairs=(), stricts=(), existing=()):
    split, mates, groups = _results(events, pairs, stricts)
    return W.weak_site_calls("s", split, mates, groups, existing)


def test_pooled_one_read_each_side_makes_a_call():
    out = _run(pairs=[_pair("chr1", "roo", 1000, 1010, left=[1], right=[2])])
    assert len(out.calls) == 1
    call = out.calls[0]
    assert (call.legacy_event_id, call.chrom, call.family, call.evidence_type) == ("WEAK_000001", "chr1", "roo", W.EVIDENCE_TYPE)
    assert call.call_class is IncumbentCallClass.WEAK_SITE
    assert (call.left_support, call.right_support, call.total_support) == (1, 1, 2)
    assert call.fragment_pair_ids == (1, 2) and call.score == 2.0
    assert (call.start, call.end) == (1000, 1011)       # as the desk replay: int(low), int(high) + 1


def test_pooling_runs_across_the_stages():
    out = _run(events=[_event("one_sided_left", "chr1", 1000, 1001, "roo", [5])],
               stricts=[_strict("chr1", "roo", 1020, 1030, right=[6])])
    assert len(out.calls) == 1 and out.calls[0].fragment_pair_ids == (5, 6)


def test_one_side_only_needs_two_split_reads():
    assert _run(pairs=[_pair("chr1", "roo", 1000, 1010, left=[1, 2, 3])]).calls == ()      # pair reads on one side: no call
    one = _run(events=[_event("one_sided_left", "chr1", 1000, 1001, "roo", [5])])
    assert one.calls == ()
    two = _run(events=[_event("one_sided_right", "chr1", 1000, 1001, "roo", [5, 6])])
    assert len(two.calls) == 1 and (two.calls[0].left_support, two.calls[0].right_support) == (0, 2)
    pooled = _run(events=[_event("one_sided_left", "chr1", 1000, 1001, "roo", [5]),
                          _event("one_sided_left", "chr1", 1050, 1051, "roo", [6])])           # two split reads on one side, 2 events
    assert len(pooled.calls) == 1


def test_a_call_within_500_bp_blocks_the_site():
    existing = [IncumbentCall("s", "X", IncumbentCallClass.MATE_RESCUE, "chr1", 1400, 1401, "roo", None, None, None, 2, 2, "m", 4.0, (1, 2, 3, 4))]
    pair = _pair("chr1", "roo", 900, 910, left=[1], right=[2])
    assert _run(pairs=[pair], existing=existing).calls == ()
    far = IncumbentCall("s", "X", IncumbentCallClass.MATE_RESCUE, "chr1", 1412, 1413, "roo", None, None, None, 2, 2, "m", 4.0, (1, 2, 3, 4))
    assert len(_run(pairs=[pair], existing=[far]).calls) == 1                     # site ends 911; a call from 1411 on is within 500 bp, 1412 is not
    other = IncumbentCall("s", "X", IncumbentCallClass.MATE_RESCUE, "chr2", 900, 901, "roo", None, None, None, 2, 2, "m", 4.0, (1, 2, 3, 4))
    assert len(_run(pairs=[pair], existing=[other]).calls) == 1                   # another contig does not count


def test_site_grouping_gap_is_100_bp():
    near = _run(pairs=[_pair("chr1", "roo", 1000, 1010, left=[1]), _pair("chr1", "roo", 1110, 1120, right=[2])])
    assert len(near.calls) == 1 and near.sites == 1             # start 1110 <= end 1011 + 100
    far = _run(pairs=[_pair("chr1", "roo", 1000, 1010, left=[1]), _pair("chr1", "roo", 1112, 1120, right=[2])])
    assert far.calls == () and far.sites == 2


def test_family_is_the_one_with_most_ids_and_span_is_its_items():
    out = _run(pairs=[_pair("chr1", "roo", 1000, 1010, left=[1], right=[2]),
                      _pair("chr1", "gypsy", 1050, 1060, left=[3, 4], right=[5, 6])],
               events=[_event("two_sided", "chr1", 1100, 1105, "gypsy", [9])])
    (call,) = out.calls
    assert call.family == "gypsy" and (call.start, call.end) == (1050, 1105)
    assert call.fragment_pair_ids == (1, 2, 3, 4, 5, 6)         # the pooled ids of the rule; a two-sided event gives none
    assert (call.left_support, call.right_support) == (3, 3)
    tie = _run(pairs=[_pair("chr1", "roo", 1000, 1010, left=[1], right=[2]), _pair("chr1", "gypsy", 1000, 1010, left=[3], right=[4])])
    assert tie.calls[0].family == "gypsy"                        # ties: the name that sorts first


def test_calls_are_numbered_in_site_order():
    out = _run(pairs=[_pair("chr2", "roo", 50, 60, left=[1], right=[2]), _pair("chr1", "roo", 900, 910, left=[3], right=[4])])
    assert [(c.chrom, c.legacy_event_id) for c in out.calls] == [("chr1", "WEAK_000001"), ("chr2", "WEAK_000002")]


def test_weak_call_origin_and_tier():
    assert callset.CALL_ORIGINS[W.EVIDENCE_TYPE] == "INCUMBENT_WEAK_SITE"
    assert tiers.evidence_class(W.EVIDENCE_TYPE) == "WEAK_SITE"
    for support in (2, 4, 40):
        assert tiers.assign_tier("EVITEA", W.EVIDENCE_TYPE, support, inside_own_family=False, near_any=False) == 4
    assert tiers.assign_tier("EVITEA", W.EVIDENCE_TYPE, 9, inside_own_family=True, near_any=False) == 5
    assert tiers.assign_tier("EVITEA", "two_sided", 9, inside_own_family=False, near_any=False) == 1     # others unchanged


def test_config_validation_and_manifest_default():
    assert DepthFloorConfig().is_default and not DepthFloorConfig().needs_depth
    assert "depth_floors" not in EviTEAConfig().normalized()
    on = EviTEAConfig(depth_floors=DepthFloorConfig(weak_site_pass=True, weak_site_depth_cut=10))
    assert on.normalized()["depth_floors"]["weak_site_depth_cut"] == 10
    with pytest.raises(InputContractError):
        DepthFloorConfig(weak_site_pass=True)
    with pytest.raises(InputContractError):
        DepthFloorConfig(weak_site_depth_cut=10)
    for bad in ({"host_cluster_support": 0}, {"candidate_min_pairs": True}, {"candidate_floor_depth_cut": 0},
                {"weak_site_pass": 1}):
        with pytest.raises(InputContractError):
            DepthFloorConfig(**bad)
    assert DepthFloorConfig(candidate_floor_depth_cut=9).needs_depth and not DepthFloorConfig(candidate_min_pairs=1).needs_depth


def _inputs(**kw):
    return IncumbentInputs("s", "te.fa", "map", "ref.fa", "host", "bowtie2", **kw)


def test_incumbent_fingerprint_default_is_unchanged():
    base = _fingerprint(_inputs())
    assert "host_cluster_support" not in base and "candidate_min_pairs" not in base
    assert _fingerprint(_inputs(host_cluster_support=2, candidate_min_pairs=2)) == base
    changed = _fingerprint(_inputs(host_cluster_support=1, candidate_min_pairs=1))
    assert changed != base and changed["host_cluster_support"] == 1 and changed["candidate_min_pairs"] == 1
    with pytest.raises(ValueError):
        _inputs(candidate_min_pairs=0)


def test_engine_identity_records_floors_only_when_on():
    def engine(**kw):
        p = "x"
        return EngineInputs("s", p, p, p, p, p, p, p, p, p, **kw)

    assert "depth_floors" not in engine().identity()
    assert engine(depth_floors=DepthFloorConfig(candidate_min_pairs=1)).identity()["depth_floors"]["candidate_min_pairs"] == 1


def test_resolve_floors(monkeypatch):
    say = []
    monkeypatch.setattr(DF, "sample_depth", lambda *a, **k: 8.0)
    cfg = DepthFloorConfig(weak_site_pass=True, weak_site_depth_cut=10, host_cluster_support=1, candidate_min_pairs=1,
                           candidate_floor_depth_cut=10)
    got = DF.resolve_floors(cfg, Path("r"), Path("i"), Path("w"), say.append)
    assert (got.depth, got.host_cluster_support, got.candidate_min_pairs, got.weak_site_pass) == (8.0, 1, 1, True)
    rec = got.record(cfg)
    assert rec["single_copy_depth"] == 8.0 and rec["weak_site_pass"] == {"depth_cut": 10, "ran": True}
    assert rec["candidate_floor"]["host_cluster_support_used"] == 1 and rec["candidate_floor"]["depth_cut"] == 10
    monkeypatch.setattr(DF, "sample_depth", lambda *a, **k: 10.5)
    got = DF.resolve_floors(cfg, Path("r"), Path("i"), Path("w"), say.append)
    assert (got.host_cluster_support, got.candidate_min_pairs, got.weak_site_pass) == (2, 2, False)
    assert got.record(cfg)["weak_site_pass"]["ran"] is False
    monkeypatch.setattr(DF, "sample_depth", lambda *a, **k: pytest.fail("no depth cut: depth must not be measured"))
    got = DF.resolve_floors(DepthFloorConfig(candidate_min_pairs=1), Path("r"), Path("i"), Path("w"), say.append)
    assert (got.depth, got.host_cluster_support, got.candidate_min_pairs) == (None, 2, 1)


def test_acquisition_threads_the_candidate_minimums(monkeypatch):
    seen = {}
    ok = mock.MagicMock()
    for name in ("build_anchor_catalogue", "census_anchor_intervals", "select_sparse_panel", "route_sparse_panel",
                 "short_exact_census", "hard_gap_h1_loci", "optimize_rescue_panel", "select_pair_aware_panel",
                 "resident_pair_families", "plan_rescue_anchors", "enumerate_rescue_pairs", "rescue_host_requests",
                 "parse_rescue_host_sam", "run_host_mapping", "parse_primary_host_sam", "one_mate_host_requests"):
        monkeypatch.setattr(ACQ, name, ok)
    monkeypatch.setattr("evitea.incumbent.coverage_panel.baseline_anchor_loci", ok)
    monkeypatch.setattr(ACQ, "one_mate_host_requests", lambda *a, **k: {})
    monkeypatch.setattr(ACQ, "rescue_host_requests", lambda *a, **k: {})
    monkeypatch.setattr(ACQ, "supported_rescue_pairs", lambda *a, **k: set())
    monkeypatch.setattr(ACQ, "working_candidate_pairs", lambda *a, **k: np.array([], dtype=np.uint64))
    monkeypatch.setattr(ACQ, "supported_host_cluster_pairs", lambda primary, **k: seen.setdefault("support", k.get("support", "default")) and set())
    monkeypatch.setattr(ACQ, "reconstruct_candidate_pairs", lambda **k: ((), {}))
    monkeypatch.setattr(ACQ, "cluster_candidate_loci",
                        lambda candidates, **k: seen.setdefault("min_pairs", k.get("min_pairs", "default")) and ())
    index = mock.MagicMock()
    index.mate_dollars.return_value = {}
    index.recover.return_value = ""

    def run(**kw):
        seen.clear()
        return ACQ.acquire_candidate_loci(te_records={}, te_family_map={}, index=index, bowtie2="b", host_index_prefix="h",
                                          work_directory="w", **kw)

    default = run()
    assert seen == {"support": 2, "min_pairs": 2}
    assert "candidate_pairs" not in default.counts and "candidate_step_seconds" not in default.counts   # default counts unchanged
    low = run(host_cluster_support=1, candidate_min_pairs=1)
    assert seen == {"support": 1, "min_pairs": 1}
    assert low.counts["candidate_pairs"] == 0 and "candidate_step_seconds" in low.counts


def test_run_engine_with_the_floors_on(tmp_path, monkeypatch):
    """Floors resolved from the depth, candidate minimums reach the library layer, the weak-site rows join the layer's calls and
    the extension layers' base table, and the record lands in the engine counts. With the floors off the existing engine test
    (test_default_engine) shows an unchanged call."""

    import evitea.engine as engine_module
    from evitea.incumbent.acquisition import CandidateAcquisition
    from evitea.incumbent.engine import IncumbentResult
    from evitea.sample_geometry import FragmentSizes
    from test_default_engine import _fake_extension, _inputs as engine_inputs

    seen = {}
    pair = _pair("chr1", "roo", 1000, 1010, left=[1], right=[2])
    split, mates, groups = _results(pairs=[pair])

    def incumbent(inputs, sample, work, progress=None):
        seen["inputs"] = inputs
        acquisition = CandidateAcquisition((), (), {"working_candidate_pairs": 7, "candidate_step_seconds": 1.5})
        return IncumbentResult((), acquisition, split, mates, None, None, groups, (), {"incumbent_calls": 0})

    def extension(anchor_inputs, internal_codebook, sample, incumbent_bytes, labels, work, progress=None):
        seen["base"] = incumbent_bytes
        return _fake_extension(incumbent_bytes, labels)

    monkeypatch.setattr(engine_module, "MsbwtSample", lambda path: ("msbwt", path))
    monkeypatch.setattr(engine_module, "run_incumbent_caller", incumbent)
    monkeypatch.setattr(engine_module, "run_extension", extension)
    monkeypatch.setattr(engine_module, "estimate_fragment_sizes", lambda *a, **k: FragmentSizes(300.0, 420, 5000, "estimate"))
    monkeypatch.setattr(DF, "sample_depth", lambda *a, **k: 6.0)
    floors = DepthFloorConfig(weak_site_pass=True, weak_site_depth_cut=10, host_cluster_support=1, candidate_min_pairs=1,
                              candidate_floor_depth_cut=8)
    base = engine_inputs(tmp_path)
    inputs = EngineInputs(**{**{f: getattr(base, f) for f in base.__slots__}, "depth_floors": floors})
    result = engine_module.run_engine(inputs, tmp_path / "work")
    assert (seen["inputs"].host_cluster_support, seen["inputs"].candidate_min_pairs) == (1, 1)
    assert b"WEAK_000001" in seen["base"] and b"weak_site_pooled_floor" in seen["base"]
    record = result.counts["depth_floors"]
    assert record["single_copy_depth"] == 6.0
    assert record["weak_site_pass"] == {"depth_cut": 10, "ran": True, "sites": 1, "uncovered_sites": 1, "rows_added": 1,
                                        "calls_before": 0}
    assert record["candidate_floor"] == {
        "host_cluster_support_requested": 1, "candidate_min_pairs_requested": 1, "depth_cut": 8,
        "host_cluster_support_used": 1, "candidate_min_pairs_used": 1, "candidate_pairs": 0, "candidate_places": 0,
        "working_candidate_pairs": 7, "wall_seconds": 1.5}
    assert result.layer_rows["incumbent"] == 1

    # depth above both cuts: the minimums fall back to 2 and the pass does not run
    monkeypatch.setattr(DF, "sample_depth", lambda *a, **k: 12.0)
    seen.clear()
    result = engine_module.run_engine(inputs, tmp_path / "work2")
    assert (seen["inputs"].host_cluster_support, seen["inputs"].candidate_min_pairs) == (2, 2)
    assert b"WEAK_" not in seen["base"]
    assert result.counts["depth_floors"]["weak_site_pass"] == {"depth_cut": 10, "ran": False, "rows_added": 0, "calls_before": 0}


def test_run_command_passes_the_floors_only_when_given(tmp_path):
    from evitea.cli import build_parser
    from evitea.run_command import run

    def steps(*extra, library=True):
        seen = []
        args = build_parser().parse_args([
            "run", "--reference", "ref.fa", "--repeat-annotation", "rm.bed", "--read1", "a", "--read2", "b", "--sample-id", "s",
            "--output", str(tmp_path / f"o{len(list(tmp_path.iterdir()))}"), "--threads", "2",
            *(["--te-library", "te.fa", "--family-map", "m.tsv", "--library-free"] if library else ["--no-library"]), *extra])
        assert run(args, lambda words: seen.append(list(words)) or 0, lambda m: None) == 0
        return {w[0]: w for w in seen}

    floors = ("--extension-stop-min", "--extension-stop-min-depth-cut", "--weak-site-pass", "--weak-site-depth-cut",
              "--host-cluster-support", "--candidate-min-pairs", "--candidate-floor-depth-cut")
    for library in (True, False):
        plain = steps(library=library)
        assert not any(flag in words for words in plain.values() for flag in floors)         # nothing new by default
        explicit_default = steps("--extension-stop-min", "3", library=library)
        assert "--extension-stop-min" not in explicit_default["library-free"]                 # 3 is the default: not passed
    on = steps("--extension-stop-min", "2", "--extension-stop-min-depth-cut", "12", "--weak-site-pass",
               "--weak-site-depth-cut", "10", "--host-cluster-support", "1", "--candidate-min-pairs", "1",
               "--candidate-floor-depth-cut", "10")
    lf, call = on["library-free"], on["call"]
    assert lf[lf.index("--extension-stop-min") + 1] == "2" and float(lf[lf.index("--extension-stop-min-depth-cut") + 1]) == 12.0
    assert "--weak-site-pass" in call and call[call.index("--weak-site-depth-cut") + 1] == "10"
    assert call[call.index("--host-cluster-support") + 1] == "1" and call[call.index("--candidate-min-pairs") + 1] == "1"
    assert call[call.index("--candidate-floor-depth-cut") + 1] == "10"
    nolib = steps("--extension-stop-min", "2", library=False)["library-free"]
    assert nolib[nolib.index("--extension-stop-min") + 1] == "2"


def test_call_command_reads_the_floors_into_the_config():
    from evitea.cli import _call_config as _config_from_args, build_parser

    parser = build_parser()
    base = ["call", "--sample-assets", "s", "--reference-assets", "r", "--te-assets", "t", "--output", "o"]
    args = parser.parse_args(base)
    assert _config_from_args(args).depth_floors.is_default
    args = parser.parse_args([*base, "--weak-site-pass", "--weak-site-depth-cut", "10", "--host-cluster-support", "1",
                              "--candidate-min-pairs", "1", "--candidate-floor-depth-cut", "9.5"])
    floors = _config_from_args(args).depth_floors
    assert floors == DepthFloorConfig(True, 10.0, 1, 1, 9.5)
    with pytest.raises(InputContractError):
        _config_from_args(parser.parse_args([*base, "--weak-site-pass"]))
