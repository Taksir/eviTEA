"""Weak-site pass of the library-free step, switched by the allele depth Da = D x (1 - median e): constructed tables.

The run-level checks on the synthetic sample (default files unchanged, saved pieces reused, combined files) are in
``tools/library_free_tests/test_depth_floors.py``.
"""

from __future__ import annotations

import json
import pickle
from types import SimpleNamespace

import pytest

from evitea.errors import InputContractError
from evitea.incumbent.mate_rescue import MateRescueLocus, MateRescueResult
from evitea.incumbent.split_junction import SplitJunctionResult
from evitea.library_free import het_guard as HG
from evitea.library_free import weak_site as WS
from evitea.resident_filter import DroppedRow

HEAD = "event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\tleft_support\tright_support\ttotal_support\t" \
       "evidence_type\tscore\tfragment_ids\n"
NATIVE = (HEAD + "N1\tchr1\t10000\t10010\troo\t['roo']\t\t\t\t2\t2\t4\ttwo_sided_mate_cluster_rescue\t4.0\t1,2,3,4\n").encode()
ADDED = b"event_id\tchr\tstart\tend\tfamily\ttsd\nLF_000001\tchr1\t20000\t20050\tgypsy\t5\n"


def _pair(chrom, family, a, b, left=(), right=()):
    return MateRescueLocus("L", chrom, family, float(a), float(b), float(a + b) / 2, tuple(left), tuple(right), False)


def _tables(*pairs):
    return (SplitJunctionResult((), (), (), {}), MateRescueResult((), tuple(pairs), (), {}), SimpleNamespace(ledger=()))


def _weak(cut, depth, e, pairs):
    return WS.WeakSite(cut, depth, _tables(*pairs), list(e))


def test_allele_depth_median_e_and_the_fewer_than_ten_rule():
    ten = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.6]       # median 0.45
    got = WS.allele_depth(20.0, ten, HG.NMIN)
    assert (got["added_calls_used"], got["median_e"]) == (10, 0.45) and got["a"] == pytest.approx(0.55)
    assert got["allele_depth"] == pytest.approx(11.0) and got["single_copy_depth"] == 20.0
    nine = WS.allele_depth(20.0, ten[:9], HG.NMIN)
    assert nine["a"] == 1.0 and nine["allele_depth"] == 20.0 and nine["added_calls_used"] == 9       # fewer than 10: a = 1
    none = WS.allele_depth(20.0, [], HG.NMIN)
    assert none["median_e"] is None and none["a"] == 1.0
    assert HG.NMIN == 10
    assert got["e_source"] == "added library-free calls"


def test_the_e_values_are_those_of_the_added_library_free_calls():
    """runtime collects e from the ADDED rows of both passes that have a measured e, and from no other row."""

    import math

    from evitea.library_free import runtime as R

    rows = [SimpleNamespace(status="ADDED", e=0.4, scan="first"), SimpleNamespace(status="ADDED", e=0.2, scan="second"),
            SimpleNamespace(status="ADDED", e=math.nan, scan="first"), SimpleNamespace(status="NEAR_EVITEA", e=0.9, scan="first"),
            SimpleNamespace(status="DROPPED_HOM_REGION", e=0.5, scan="first")]
    weak = _weak(100.0, 10.0, [], [])
    src = open(R.__file__).read()
    assert 'weak_site.added_e = [r.e for r in rows if r.status == "ADDED" and not math.isnan(r.e)]' in src
    weak.added_e = [r.e for r in rows if r.status == "ADDED" and not math.isnan(r.e)]       # the same expression
    assert weak.added_e == [0.4, 0.2]


def test_skipped_above_the_cut_and_run_at_or_under_it():
    site = _pair("chr1", "roo", 5000, 5010, left=[1], right=[2])
    e = [0.5] * 10                                                    # a = 0.5 (ten added calls with e = 0.5)
    shown, dropped, info = WS.apply(_weak(10.0, 20.0, e, [site]), NATIVE, ADDED, None, "s")      # Da = 10: at the cut
    assert info["ran"] and info["rows_added"] == 1 and b"WEAK_000001" in shown and shown.startswith(NATIVE)
    row = shown.decode().splitlines()[-1].split("\t")       # a native table with the evidence columns keeps them for the new row
    assert row[0] == "WEAK_000001" and row[12] == "weak_site_pooled_floor" and (row[9], row[10], row[11]) == ("1", "1", "2")
    shown, dropped, info = WS.apply(_weak(9.99, 20.0, e, [site]), NATIVE, ADDED, None, "s")      # Da = 10 > 9.99
    assert not info["ran"] and shown == NATIVE and "rows_added" not in info and info["cut"] == 9.99
    assert info["allele_depth"] == pytest.approx(10.0) and info["median_e"] == 0.5 and info["a"] == 0.5
    shown, _d, info = WS.apply(_weak(10.0, 20.0, [0.1] * 10, [site]), NATIVE, ADDED, None, "s")  # a = 0.9: Da = 18 > 10
    assert not info["ran"] and shown == NATIVE


def test_library_rows_and_library_free_rows_both_cover_a_site():
    near_library = _pair("chr1", "roo", 10400, 10410, left=[1], right=[2])      # 390 bp from N1 (10000-10010)
    near_added = _pair("chr1", "roo", 20450, 20460, left=[3], right=[4])        # within 500 bp of LF_000001 (20000-20050)
    free = _pair("chr1", "roo", 30000, 30010, left=[5], right=[6])
    shown, _d, info = WS.apply(_weak(100.0, 5.0, [], [near_library, near_added, free]), NATIVE, ADDED, None, "s")
    assert info["sites"] == 3 and info["uncovered_sites"] == 1 and info["rows_added"] == 1 and info["cover_rows"] == 2
    assert shown.decode().count("WEAK_") == 1 and "30000" in shown.decode().split("WEAK_")[1]
    only_added = WS.apply(_weak(100.0, 5.0, [], [near_added]), NATIVE, ADDED, None, "s")[2]
    assert only_added["rows_added"] == 0                 # the added library-free call alone covers it


def test_own_family_copy_filter_drops_a_row_inside_a_copy_of_its_own_family():
    inside = _pair("chr1", "roo", 40000, 40010, left=[1], right=[2])
    outside = _pair("chr1", "roo", 50000, 50010, left=[3], right=[4])

    class Filter:
        def __call__(self, table):
            lines = table.decode().splitlines()
            kept, dropped = [lines[0]], []
            for line in lines[1:]:
                f = line.split("\t")
                if f[1] == "chr1" and int(f[2]) == 40000:
                    dropped.append(DroppedRow(f[0], f[1], f[2], f[3], f[4], "FBte1", "chr1", 39000, 41000, 0.1))
                else:
                    kept.append(line)
            return SimpleNamespace(filtered=("\n".join(kept) + "\n").encode(), dropped=tuple(dropped))

    shown, dropped, info = WS.apply(_weak(100.0, 5.0, [], [inside, outside]), NATIVE, ADDED, Filter(), "s")
    assert (info["passing_rule"], info["dropped_by_copy_filter"], info["rows_added"]) == (2, 1, 1)
    assert [d.event_id for d in dropped] == ["WEAK_000001"]
    text = shown.decode()
    assert "\t50000\t" in text and "\t40000\t" not in text and text.count("WEAK_") == 1


def _stage_tree(root, tables, work_in_manifest=True):
    run = root / "run"
    run.mkdir()
    work = root / "callwork"
    folder = work / "engine" / "work" / "incumbent" / "checkpoints"
    folder.mkdir(parents=True)
    for stage, table in zip(WS.STAGES, tables):
        payload = pickle.dumps(table)
        (folder / f"{stage}.pkl").write_bytes(payload)
        import hashlib
        (folder / f"{stage}.complete.json").write_text(json.dumps({"fingerprint": {}, "sha256": hashlib.sha256(payload).hexdigest()}))
    runtime = {"work_dir": str(work) if work_in_manifest else None}
    (run / "run_manifest.json").write_text(json.dumps({"values": {"normalized_config": {"runtime": runtime}}}))
    return run, folder


def test_stage_tables_are_found_from_the_manifest_and_missing_ones_stop_the_run(tmp_path):
    run, folder = _stage_tree(tmp_path, _tables(_pair("chr1", "roo", 1, 2, left=[1], right=[2])))
    assert WS.stage_directory(run) == folder
    split, mates, groups = WS.load_stage_tables(run)
    assert mates.ledger[0].left_pair_ids == (1,) and groups.ledger == ()
    (folder / "mate_rescue.pkl").unlink()
    with pytest.raises(InputContractError, match="mate_rescue"):
        WS.load_stage_tables(run)
    (folder / "mate_rescue.pkl").write_bytes(b"x")
    with pytest.raises(InputContractError, match="does not match"):
        WS.load_stage_tables(run)
    with pytest.raises(InputContractError, match="run manifest"):
        WS.load_stage_tables(tmp_path / "nowhere")
    run2 = tmp_path / "r2"
    run2.mkdir()
    (run2 / "run_manifest.json").write_text(json.dumps({"values": {"normalized_config": {"runtime": {"work_dir": None}}}}))
    assert WS.stage_directory(run2) == run2 / ".evitea-work" / "engine" / "work" / "incumbent" / "checkpoints"


def test_run_command_routes_the_pass(tmp_path):
    from evitea.cli import build_parser
    from evitea.run_command import run

    def steps(*extra, library=True):
        seen = []
        args = build_parser().parse_args([
            "run", "--reference", "ref.fa", "--repeat-annotation", "rm.bed", "--read1", "a", "--read2", "b", "--sample-id", "s",
            "--output", str(tmp_path / f"o{len(list(tmp_path.iterdir()))}"), "--threads", "2",
            *(["--te-library", "te.fa", "--family-map", "m.tsv"] if library else ["--no-library"]), *extra])
        assert run(args, lambda words: seen.append(list(words)) or 0, lambda m: None) == 0
        return {w[0]: w for w in seen}

    got = steps("--library-free", "--weak-site-pass", "--weak-site-allele-depth-cut", "10")
    assert "--weak-site-pass" not in got["call"]                      # the pass belongs to the library-free step
    lf = got["library-free"]
    assert lf[lf.index("--weak-site-allele-depth-cut") + 1] == "10.0" and "--weak-site-pass" in lf
    call_route = steps("--library-free", "--weak-site-pass", "--weak-site-depth-cut", "10")
    assert "--weak-site-pass" in call_route["call"] and "--weak-site-pass" not in call_route["library-free"]
    for bad in (("--weak-site-allele-depth-cut", "10"),                                   # without the pass
                ("--weak-site-pass", "--weak-site-allele-depth-cut", "10"),               # without --library-free
                ("--library-free", "--weak-site-pass", "--weak-site-allele-depth-cut", "10", "--weak-site-depth-cut", "10")):
        with pytest.raises(InputContractError):
            steps(*bad)
    with pytest.raises(InputContractError, match="no-library"):
        steps("--weak-site-pass", "--weak-site-allele-depth-cut", "10", library=False)
    plain = steps("--library-free")
    assert not any(word.startswith("--weak") for words in plain.values() for word in words)


def test_cli_reads_the_options():
    from evitea.cli import build_parser

    args = build_parser().parse_args(["library-free", "--run-dir", "r", "--sample-assets", "s", "--reference-assets", "x", "--output", "o"])
    assert args.weak_site_pass is False and args.weak_site_allele_depth_cut is None
    args = build_parser().parse_args(["library-free", "--run-dir", "r", "--sample-assets", "s", "--reference-assets", "x", "--output", "o",
                                      "--weak-site-pass", "--weak-site-allele-depth-cut", "8.5"])
    assert args.weak_site_pass is True and args.weak_site_allele_depth_cut == 8.5
