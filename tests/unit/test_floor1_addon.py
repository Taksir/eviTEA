"""``library-free --floor1-addon``: the add rule, the merge into the main output, the switch, and agreement with the desk script."""

from __future__ import annotations

import csv
import hashlib
import subprocess
import sys
import tarfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from evitea import vcf as V
from evitea.errors import InputContractError
from evitea.library_free import floor1_addon as F1
from evitea.library_free import runtime as R


def row(eid, chrom, start, end, source="LIBRARY_FREE", family="LF_cluster_0"):
    return {"event_id": eid, "chr": chrom, "start": str(start), "end": str(end), "family": family, "source": source}


def test_the_500_bp_edge_on_both_sides():
    main = [row("M1", "chr2L", 10_000, 10_010, "EVITEA")]
    # a gap of exactly 500 is near, 501 is not (right of the main row, and left of it)
    assert F1.added_rows(main, [row("A", "chr2L", 10_510, 10_511)]) == []
    assert [r["event_id"] for r in F1.added_rows(main, [row("A", "chr2L", 10_511, 10_512)])] == ["F1_A"]
    assert F1.added_rows(main, [row("A", "chr2L", 9_499, 9_500)]) == []
    assert [r["event_id"] for r in F1.added_rows(main, [row("A", "chr2L", 9_498, 9_499)])] == ["F1_A"]
    # an add-on interval that spans the main row is near
    assert F1.added_rows(main, [row("A", "chr2L", 5_000, 20_000)]) == []


def test_other_chromosome_other_sources_prefix_order_and_no_comparison_among_add_on_rows():
    main = [row("M1", "chr2L", 10_000, 10_010, "EVITEA"), row("M2", "chr3R", 500, 501, "LIBRARY_FREE")]
    addon = [row("B", "chr3R", 10_000, 10_001),                        # same position, other chromosome: added
             row("C", "chr2L", 50_000, 50_001, "EVITEA"),             # a weak-site / library row of the add-on: not added
             row("D", "chr3R", 20_000, 20_001),
             row("E", "chr3R", 20_100, 20_101),                       # near D, but add-on rows are not compared with each other
             row("F", "chrX", 1, 2),                                   # a chromosome the main table does not have
             row("G", "chr3R", 1_000, 1_001)]                         # within 500 of M2 (a LIBRARY_FREE main row counts)
    got = F1.added_rows(main, addon)
    assert [r["event_id"] for r in got] == ["F1_B", "F1_D", "F1_E", "F1_F"]
    assert all(r["source"] == "LIBRARY_FREE" for r in got) and got[0]["chr"] == "chr3R" and got[0]["start"] == "10000"
    assert F1.added_rows([], []) == [] and F1.added_rows(main, []) == []


def _fasta(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr2L\n" + "ACGT" * 5000 + "\n>chr3R\n" + "TGCA" * 5000 + "\n")
    return fasta


def _vcf(fasta, rows):
    calls = V.library_free_vcf_calls(_lf_bytes(rows))
    contigs = V.fasta_contigs(fasta)
    return V.calls_vcf_bytes(calls, contigs, sample_id="S", bases=V.reference_bases(fasta, calls), reference_name=fasta.name)


LF_HEAD = "event_id\tchr\tstart\tend\tfamily\tcluster\tbL\tbR\ttsd\te\tc"


def _lf_line(eid, chrom, pos, cluster):
    return f"{eid}\t{chrom}\t{pos}\t{pos + 1}\tLF_cluster_{cluster}\t{cluster}\t{pos}\t{pos - 5}\t5\t0.0\t1.0"


def _lf_bytes(rows):
    return ("\n".join([LF_HEAD, *rows]) + "\n").encode()


def test_vcf_merge_equals_a_vcf_written_from_all_records(tmp_path):
    fasta = _fasta(tmp_path)
    main = [_lf_line("LF_chr3R_9000", "chr3R", 9000, 1), _lf_line("LF_chr2L_2000", "chr2L", 2000, 0)]
    new = [_lf_line("F1_LF_chr2L_1500", "chr2L", 1500, 7), _lf_line("F1_LF_chr3R_12000", "chr3R", 12000, 0)]
    merged = F1.merge_vcf(_vcf(fasta, main), _vcf(fasta, new))
    assert merged == _vcf(fasta, [*main, *new])


def _write_main(tmp_path, depth, cut):
    out = tmp_path / "out"
    out.mkdir()
    fasta = _fasta(tmp_path)
    main_lf = [_lf_line("LF_chr2L_2000", "chr2L", 2000, 0)]
    (out / "calls.library_free.tsv").write_bytes(_lf_bytes(main_lf))
    (out / "calls.combined.tsv").write_text(
        "event_id\tchr\tstart\tend\tfamily\tsource\nN1\tchr2L\t6000\t6005\troo\tEVITEA\nLF_chr2L_2000\tchr2L\t2000\t2001\tLF_cluster_0\tLIBRARY_FREE\n")
    (out / "calls.combined.vcf").write_bytes(_vcf(fasta, main_lf))
    summary = {"weak_site_pass": {"allele_depth": depth, "cut": cut, "ran": depth <= cut}, "output_hashes": {
        n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in ("calls.combined.tsv", "calls.library_free.tsv", "calls.combined.vcf")}}
    return out, fasta, summary


def test_floor1_addon_step_merges_rows_into_every_main_file_and_updates_the_summary(tmp_path, monkeypatch):
    out, fasta, summary = _write_main(tmp_path, depth=6.0, cut=10.0)
    seen = {}

    def fake_run(**kw):
        seen.update(kw)
        folder = Path(kw["output"])
        folder.mkdir(parents=True)
        (folder / "calls.library_free.tsv").write_bytes(_lf_bytes([
            _lf_line("LF_chr2L_2000", "chr2L", 2000, 0), _lf_line("LF_chr2L_2300", "chr2L", 2300, 1),      # both within 500 of main rows
            _lf_line("LF_chr2L_9000", "chr2L", 9000, 2), _lf_line("LF_chr3R_400", "chr3R", 400, 3)]))
        (folder / "calls.combined.tsv").write_text(
            "event_id\tchr\tstart\tend\tfamily\tsource\nW1\tchr2L\t9500\t9501\troo\tEVITEA\n"
            "LF_chr2L_2000\tchr2L\t2000\t2001\tLF_cluster_0\tLIBRARY_FREE\nLF_chr2L_2300\tchr2L\t2300\t2301\tLF_cluster_1\tLIBRARY_FREE\n"
            "LF_chr2L_9000\tchr2L\t9000\t9001\tLF_cluster_2\tLIBRARY_FREE\nLF_chr3R_400\tchr3R\t400\t401\tLF_cluster_3\tLIBRARY_FREE\n")

    monkeypatch.setattr(R, "run_library_free", fake_run)
    args = {"run_dir": tmp_path, "msbwt_index": tmp_path, "contigs": None, "workers": 2}
    reference = SimpleNamespace(reference_fasta_path=fasta)
    R._floor1_addon(summary, out, tmp_path / "work", args, reference, "S", lambda m: None)
    # the second step: same arguments, floor 1, a scan of its own, add-on off, its own folders
    assert seen["extension_stop_min"] == 1 and seen["extension_stop_min_depth_cut"] is None and seen["detection_dir"] is None
    assert seen["floor1_addon"] is False and seen["weak_site_pass"] is False and seen["weak_site_allele_depth_cut"] is None
    assert seen["output"] == out / "floor1_addon" and seen["work_dir"] == tmp_path / "work" / "floor1_addon" and seen["workers"] == 2
    rows = F1.read_rows((out / "calls.combined.tsv").read_bytes())
    assert [r["event_id"] for r in rows] == ["N1", "LF_chr2L_2000", "F1_LF_chr2L_9000", "F1_LF_chr3R_400"]     # W1 (EVITEA) not added
    assert rows[-1]["source"] == "LIBRARY_FREE"
    lf = (out / "calls.library_free.tsv").read_text().splitlines()
    assert [x.split("\t")[0] for x in lf[1:]] == ["LF_chr2L_2000", "F1_LF_chr2L_9000", "F1_LF_chr3R_400"] and lf[0] == LF_HEAD
    vcf = [x.split("\t") for x in (out / "calls.combined.vcf").read_text().splitlines() if not x.startswith("#")]
    assert [(x[0], int(x[1]), x[2]) for x in vcf] == [("chr2L", 2000, "LF_chr2L_2000"), ("chr2L", 9000, "F1_LF_chr2L_9000"),
                                                       ("chr3R", 400, "F1_LF_chr3R_400")]
    assert "INS:ME:LF_cluster_3" in (out / "calls.combined.vcf").read_text()
    record = summary["floor1_addon"]
    assert record == {"requested": True, "ran": True, "allele_depth": 6.0, "cut": 10.0, "addon_output": str(out / "floor1_addon"),
                      "addon_library_free_rows": 4, "rows_added": 2, "rule": F1.RULE}
    assert all(summary["output_hashes"][n] == hashlib.sha256((out / n).read_bytes()).hexdigest() for n in summary["output_hashes"])
    assert (out / "library_free_summary.json").is_file()


def test_floor1_addon_does_not_run_over_the_cut(tmp_path, monkeypatch):
    out, fasta, summary = _write_main(tmp_path, depth=12.0, cut=10.0)
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    monkeypatch.setattr(R, "run_library_free", lambda **kw: pytest.fail("the second step must not run"))
    R._floor1_addon(summary, out, tmp_path / "work", {}, SimpleNamespace(reference_fasta_path=fasta), "S", lambda m: None)
    assert summary["floor1_addon"] == {"requested": True, "ran": False, "allele_depth": 12.0, "cut": 10.0, "addon_output": None,
                                       "addon_library_free_rows": None, "rows_added": 0, "rule": F1.RULE}
    assert not (out / "floor1_addon").exists()
    assert {p.name: p.read_bytes() for p in out.iterdir() if p.name != "library_free_summary.json"} == \
        {k: v for k, v in before.items() if k != "library_free_summary.json"}
    # at the cut exactly: runs
    out2, fasta2, summary2 = _write_main(tmp_path / "x", depth=10.0, cut=10.0) if (tmp_path / "x").mkdir() is None else None
    called = []
    monkeypatch.setattr(R, "run_library_free", lambda **kw: (called.append(1), Path(kw["output"]).mkdir(parents=True),
                                                             (Path(kw["output"]) / "calls.combined.tsv").write_text("event_id\tchr\tstart\tend\tfamily\tsource\n")))
    R._floor1_addon(summary2, out2, tmp_path / "work", {}, SimpleNamespace(reference_fasta_path=fasta2), "S", lambda m: None)
    assert called == [1] and summary2["floor1_addon"]["ran"] is True and summary2["floor1_addon"]["rows_added"] == 0


def test_refusals():
    ref = SimpleNamespace(reference_fasta_path="x")
    common = dict(run_dir=Path("/nonexistent"), reference=ref, output=Path("/nonexistent/o"), work_dir=None, contigs=None, rule_c=True,
                  detection_dir=None)
    with pytest.raises(InputContractError, match="--weak-site-pass"):
        R.run_library_free(**common, msbwt_index=Path("m"), floor1_addon=True)
    with pytest.raises(InputContractError, match="msBWT"):
        R.run_library_free(**common, msbwt_index=None, weak_site_pass=True, weak_site_allele_depth_cut=10.0, floor1_addon=True)


def test_command_line_and_run_passthrough(tmp_path):
    from evitea.cli import build_parser
    from evitea.run_command import library_steps, no_library_steps

    parse = lambda *a: build_parser().parse_args(["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o", *a])
    assert parse().floor1_addon is False and parse("--floor1-addon").floor1_addon is True
    base = ["run", "--reference", "ref.fa", "--repeat-annotation", "rm.bed", "--read1", "a", "--read2", "b", "--sample-id", "s",
            "--output", str(tmp_path), "--te-library", "te.fa", "--family-map", "m.tsv", "--library-free"]
    plain = dict(library_steps(build_parser().parse_args(base), tmp_path, "s"))
    assert "--floor1-addon" not in plain["library_free"]
    on = dict(library_steps(build_parser().parse_args([*base, "--weak-site-pass", "--weak-site-allele-depth-cut", "10", "--floor1-addon"]),
                            tmp_path, "s"))
    assert on["library_free"][-1] == "--floor1-addon" and "--weak-site-allele-depth-cut" in on["library_free"]
    with pytest.raises(InputContractError, match="floor1-addon"):
        library_steps(build_parser().parse_args([*base, "--floor1-addon"]), tmp_path, "s")
    nolib = ["run", "--reference", "ref.fa", "--read1", "a", "--read2", "b", "--sample-id", "s", "--output", str(tmp_path), "--no-library"]
    with pytest.raises(InputContractError, match="floor1-addon"):
        no_library_steps(build_parser().parse_args([*nolib, "--weak-site-pass", "--weak-site-allele-depth-cut", "10", "--floor1-addon"]),
                         tmp_path, "s")


# --- agreement with the desk script on the four training runs (archives at the project root, read only) ---

ROOT = Path("/mnt/project-files")
DESK = ROOT / "paper_audit_2026-11-03" / "e4_het_mix" / "dev2" / "scripts" / "union_desk.py"
RUNS = (
    ("10x", "cand_C2_10x_2026-11-03", "dev2_lf4_2026-11-04", "X1_10x", 40),
    ("15x", "cand_C2_15x_2026-11-03", "dev2_lf4_2026-11-04", "X1_15x", 16),
    ("D50_20", "cand_C2_D50_20_2026-11-03", "dev2_lf4b_2026-11-04", "X2_D50_20", 20),
    ("MIX25", "cand_C2_MIX25_2026-11-03", "dev2_lf4b_2026-11-04", "X2_MIX25", 31),
)


def _member(archive, name, tmp_path):
    with tarfile.open(archive) as handle:
        target = next(m for m in handle.getmembers() if m.name.endswith(name) and m.isfile())
        return handle.extractfile(target).read()


@pytest.mark.skipif(not DESK.is_file() or not all((ROOT / f"{run[1]}.tar.gz").is_file() for run in RUNS), reason="project archives not here")
@pytest.mark.parametrize("label,main_archive,floor_archive,folder,expected", RUNS)
def test_the_rule_agrees_with_the_desk_script(label, main_archive, floor_archive, folder, expected, tmp_path):
    main = _member(ROOT / f"{main_archive}.tar.gz", "calls.combined.tsv", tmp_path)
    floor = _member(ROOT / f"{floor_archive}.tar.gz", f"{folder}/calls.combined.tsv", tmp_path)
    (tmp_path / "main").mkdir()
    (tmp_path / "f1").mkdir()
    (tmp_path / "main" / "calls.combined.tsv").write_bytes(main)
    (tmp_path / "f1" / "calls.combined.tsv").write_bytes(floor)
    done = subprocess.run([sys.executable, str(DESK), str(tmp_path / "main"), str(tmp_path / "f1"), str(tmp_path / "out")],
                          capture_output=True, text=True, check=True)
    desk = list(csv.DictReader(open(tmp_path / "out" / "calls.combined.tsv"), delimiter="\t"))
    main_rows, floor_rows = F1.read_rows(main), F1.read_rows(floor)
    added = F1.added_rows(main_rows, floor_rows)
    desk_added = [(r["chr"], r["start"], r["end"]) for r in desk[len(main_rows):]]
    assert [(r["chr"], r["start"], r["end"]) for r in added] == desk_added
    assert len(added) == expected and f"added: {expected}" in done.stdout
    assert len(desk) == len(main_rows) + expected
