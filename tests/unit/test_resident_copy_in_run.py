"""The own-family reference copy filter that ``call``, ``library-free`` and ``run`` apply by default (BED6 copies)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from evitea.callset import CallsetState
from evitea.cli import build_parser
from evitea.config import EviTEAConfig, InputConfig, OutputConfig, RuntimeConfig
from evitea.constants import ARCHITECTURE_ID
from evitea.library_free.runtime import combined_calls_tsv_bytes, combined_calls_vcf_bytes
from evitea.pipeline import PreparedAssets, run_sample
from evitea.provenance import RunManifestContext
from evitea.resident_filter import (
    OwnFamilyCopyFilter,
    dropped_tsv_bytes,
    filter_calls,
    load_bed6_copies,
    load_copy_table,
    parse_bed6_copies,
    parse_copy_table,
    CopyIndex,
    FamilyIds,
)
from evitea.run_command import library_steps, no_library_steps, run

HEADER = ("event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\tleft_support\tright_support\t"
          "total_support\tevidence_type\tscore\tfragment_ids\n")


def _native_row(n: int, chrom: str, start: int, end: int, family: str) -> str:
    return (f"EV_{n:03d}\t{chrom}\t{start}\t{end}\t{family}\t['{family}']\t{start}\t-1\t0\t2\t0\t2\tone_sided_left\t2.0\t{n}\n")


# E1 inside a copy of its own family, E2 touches the copy's end, E3 inside a copy of another family, E4 one base away,
# E5 a family the map does not know
CALLS = [(1, "chr1", 150, 150, "B1"), (2, "chr1", 200, 200, "B1"), (3, "chr1", 150, 150, "L1"), (4, "chr1", 201, 201, "B1"),
         (5, "chr1", 700, 700, "Odd")]
NATIVE = (HEADER + "".join(_native_row(*c) for c in CALLS)).encode()
BED6 = "chr1\t100\t200\trec_B1\t500\t+\nchr1\t500\t600\trec_L1#LINE/L1\t300\t-\n# comment\n"
FAMILY_MAP = "rec_B1\tB1\nrec_L1\tL1\n"
ADDED = ("event_id\tchr\tstart\tend\tfamily\tcluster\tbL\tbR\ttsd\te\tc\n"
         "LF_chr1_150\tchr1\t150\t151\tLF_cluster_1\t1\t150\t156\t6\t0.01\t\n").encode()


@pytest.fixture()
def own_filter(tmp_path: Path) -> OwnFamilyCopyFilter:
    (tmp_path / "ann.bed").write_text(BED6, encoding="utf-8")
    (tmp_path / "map.tsv").write_text(FAMILY_MAP, encoding="utf-8")
    return OwnFamilyCopyFilter(tmp_path / "ann.bed", tmp_path / "map.tsv", ["chr1"])


def test_bed6_loader_gives_the_copy_tuples_of_the_out_loader() -> None:
    out = [
        "   SW   perc perc perc  query     position in query              matching    repeat          position in repeat",
        " score   div. del. ins.  sequence  begin    end          (left)   repeat      class/family begin   end    (left)    ID",
        "",
        "  51729    0.0  0.0  0.0  chrX             1     5707 (23536564) + FBte0000148 Unspecified    4757  10463     (0)    1  ",
        "    751   16.2  8.4  7.2  chrX          5685     6029 (23536242) C FBte0000143 Unspecified     269    617  (5466)    2 *",
        "    100   16.2  8.4  7.2  chr3L        24549286 24549611 (1) + FBte000XXXX Unspecified     269    617  (5466)    3",
    ]
    # the awk command of docs/REPEAT_ANNOTATION.md
    bed = ["chrX\t0\t5707\tFBte0000148\t51729\t+", "chrX\t5684\t6029\tFBte0000143\t751\t-",
           "chr3L\t24549285\t24549611\tFBte000XXXX\t100\t+"]
    from_out, from_bed = parse_copy_table(out), parse_bed6_copies(bed)
    assert [c[:4] for c in from_bed] == [c[:4] for c in from_out]
    assert all(c[4] == "0" for c in from_bed)


def test_bed6_names_read_as_records() -> None:
    copies = parse_bed6_copies(["track name=x", "chr1\t1\t9\tL1MdA_I#LINE/L1\t0\t+", "chr1\t1\t9\tdm_FBte0000001_roo\t0\t.",
                                "chr1\t1\t9\tFBte000XXXX\t0\t+", ""])
    assert [c[3] for c in copies] == ["L1MdA_I", "FBte0000001", "FBte000XXXX"]


def test_bed6_loader_file_matches_out_loader_file(tmp_path: Path) -> None:
    out = tmp_path / "ref.fa.out"
    out.write_text("h1\nh2\n\n    751   16.2  8.4  7.2  chr1   101   200 (5) + rec_B1 Unspecified 269 617 (5466) 2\n", encoding="utf-8")
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t100\t200\trec_B1\t751\t+\n", encoding="utf-8")
    assert [c[:4] for c in load_bed6_copies(bed)] == [c[:4] for c in load_copy_table(out)]


def test_filter_on_a_tiny_table(own_filter: OwnFamilyCopyFilter) -> None:
    result = own_filter(NATIVE)
    assert [d.event_id for d in result.dropped] == ["EV_001", "EV_002"]
    assert (result.kept, result.unmapped_family) == (3, 1)
    lines = NATIVE.decode().splitlines(keepends=True)
    assert result.filtered == "".join(lines[:1] + lines[3:]).encode()  # kept rows byte for byte, in order
    assert (result.dropped[0].family_fbte, result.dropped[0].copy_begin, result.dropped[0].copy_end) == ("rec_B1", 100, 200)
    assert dropped_tsv_bytes(result.dropped).decode().splitlines()[0].startswith("event_id\tchr\tstart\tend\tfamily")
    summary = own_filter.summary(result)
    assert (summary["enabled"], summary["kept"], summary["dropped"], summary["unmapped_family"]) == (True, 3, 2, 1)


def test_bed6_filter_equals_the_out_filter(tmp_path: Path, own_filter: OwnFamilyCopyFilter) -> None:
    out_lines = ["h", "h", "", "  500  1.0 0 0 chr1 101 200 (0) + rec_B1 X 1 100 (0) 1",
                 "  300  1.0 0 0 chr1 501 600 (0) C rec_L1 X 1 100 (0) 2"]
    index = CopyIndex(parse_copy_table(out_lines), ["chr1"])
    by_out = filter_calls(NATIVE, index, FamilyIds({"B1": ["rec_B1"], "L1": ["rec_L1"]}))
    assert by_out.filtered == own_filter(NATIVE).filtered


def test_library_free_rows_pass_untouched(own_filter: OwnFamilyCopyFilter) -> None:
    combined = combined_calls_tsv_bytes(NATIVE, ADDED)
    result = own_filter(combined)
    assert b"LF_chr1_150\tchr1\t150\t151\tLF_cluster_1\tLIBRARY_FREE" in result.filtered  # inside the B1 copy, still kept
    assert result.library_free_passed == 1 and [d.event_id for d in result.dropped] == ["EV_001", "EV_002"]


def test_filtering_the_native_table_equals_filtering_the_combined_table(own_filter: OwnFamilyCopyFilter, tmp_path: Path) -> None:
    result = own_filter(NATIVE)
    assert combined_calls_tsv_bytes(result.filtered, ADDED) == own_filter(combined_calls_tsv_bytes(NATIVE, ADDED)).filtered
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 250 + "\n", encoding="ascii")
    reference = SimpleNamespace(reference_fasta_path=fasta, contigs=(("chr1", 1000),))
    vcf = combined_calls_vcf_bytes(result.filtered, ADDED, reference, "s")
    ids = [line.split("\t")[2] for line in vcf.decode().splitlines() if not line.startswith("#")]
    assert "EV_001" not in ids and "EV_002" not in ids and "EV_003" in ids and "LF_chr1_150" in ids


def _config(tmp_path: Path) -> EviTEAConfig:
    return EviTEAConfig(inputs=InputConfig(sample_id="sample"), runtime=RuntimeConfig(threads=1, work_dir=tmp_path / "work"),
                        output=OutputConfig(tmp_path / "out"))


def _run(tmp_path: Path, resident_filter) -> tuple[Path, object]:
    config = _config(tmp_path)
    assets = PreparedAssets(reference=SimpleNamespace(contigs=(("chr1", 1000),), reference_fasta_path=None), te=SimpleNamespace())
    context = RunManifestContext(normalized_config=config.normalized(), command_line=("evitea", "call"),
                                 source_identity={"architecture": ARCHITECTURE_ID}, input_assets={},
                                 prepared_identities={"reference": "r", "te": "t"}, backend_versions={"exact_query": "x"},
                                 resolved_threads=config.runtime.resolved_threads)
    result = run_sample(config, assets, lambda *_: CallsetState("sample", NATIVE), manifest_context=context,
                        resident_filter=resident_filter)
    return result.output_directory, result


def _tree(folder: Path) -> dict[str, str]:
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir()) if p.is_file()}


def test_call_filters_the_user_tables_and_keeps_the_native_table(tmp_path: Path, own_filter: OwnFamilyCopyFilter) -> None:
    out, result = _run(tmp_path / "on", own_filter)
    assert (out / "calls.native.tsv").read_bytes() == NATIVE
    ids = [line.split("\t")[1] for line in (out / "calls.tsv").read_text().splitlines()[1:]]
    assert sorted(ids) == ["EV_003", "EV_004", "EV_005"]
    assert "EV_001" not in (out / "calls.bed").read_text() and "EV_001" not in (out / "calls.vcf").read_text()
    assert "EV_003" in (out / "calls.vcf").read_text()
    assert [line.split("\t")[0] for line in (out / "resident_copy_dropped.tsv").read_text().splitlines()[1:]] == ["EV_001", "EV_002"]
    summary = json.loads((out / "run_summary.json").read_text())
    assert summary["resident_copy_filter"]["dropped"] == 2 and summary["resident_copy_filter"]["kept"] == 3
    assert summary["final_call_count"] == 5 and summary["native_callset_sha256"] == hashlib.sha256(NATIVE).hexdigest()
    manifest = result.manifest.values
    assert "resident_copy_dropped.tsv" in manifest["output_hashes"] and manifest["output_counts"]["calls"] == 3


def test_off_switch_writes_the_outputs_of_a_run_without_the_filter(tmp_path: Path, own_filter: OwnFamilyCopyFilter) -> None:
    off, _ = _run(tmp_path / "off", None)
    names = set(_tree(off))
    assert "resident_copy_dropped.tsv" not in names
    assert b"resident_copy_filter" not in (off / "run_summary.json").read_bytes()
    assert b"resident_copy_dropped" not in (off / "run_manifest.json").read_bytes()
    assert len((off / "calls.tsv").read_text().splitlines()) == 6  # all five native calls


def test_parser_switch_and_run_steps(tmp_path: Path) -> None:
    parser = build_parser()
    for command in (["call"], ["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o"],
                    ["refine-breakpoints", "--run-dir", "r", "--sample-assets", "s", "--reference-assets", "a",
                     "--te-assets", "t", "--output", "o"]):
        assert parser.parse_args(command).keep_resident_copies is False
        assert parser.parse_args([*command, "--keep-resident-copies"]).keep_resident_copies is True
    base = ["run", "--reference", "ref.fa", "--repeat-annotation", "rm.bed", "--read1", "a", "--read2", "b", "--sample-id", "s",
            "--output", str(tmp_path / "out"), "--te-library", "te.fa", "--family-map", "m.tsv", "--library-free"]
    on = dict(library_steps(parser.parse_args(base), tmp_path / "out", "s"))
    assert "--keep-resident-copies" not in on["call"] and "--keep-resident-copies" not in on["library_free"]
    assert on["library_free"][on["library_free"].index("--te-assets") + 1] == str(tmp_path / "out" / "te")
    off = dict(library_steps(parser.parse_args([*base, "--keep-resident-copies"]), tmp_path / "out", "s"))
    assert "--keep-resident-copies" in off["call"] and "--keep-resident-copies" in off["library_free"]
    assert "--te-assets" not in off["library_free"]


def test_no_library_route_never_applies_the_filter(tmp_path: Path) -> None:
    parser = build_parser()
    for extra in ([], ["--keep-resident-copies"]):
        args = parser.parse_args(["run", "--reference", "ref.fa", "--read1", "a", "--read2", "b", "--sample-id", "s", "--output",
                                  str(tmp_path / "out"), "--no-library", *extra])
        steps = dict(no_library_steps(args, tmp_path / "out", "s"))
        assert "--keep-resident-copies" in steps["record_pass_call"]
        assert "--keep-resident-copies" not in steps["library_free"] and "--te-assets" not in steps["library_free"]
        assert "--keep-resident-copies" not in steps["record_pass_added"]


def test_library_free_refuses_library_calls_without_family_map_or_switch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import evitea.cli as cli
    from evitea.errors import InputContractError

    run_dir = tmp_path / "call"
    run_dir.mkdir()
    (run_dir / "calls.native.tsv").write_bytes(NATIVE)
    monkeypatch.setattr(cli, "load_prepared_reference", lambda _path: SimpleNamespace())
    monkeypatch.setattr(cli, "_apply_depth_cap", lambda _args: None)
    args = build_parser().parse_args(["library-free", "--run-dir", str(run_dir), "--msbwt-index", str(tmp_path),
                                      "--reference-assets", str(tmp_path), "--output", str(tmp_path / "o")])
    with pytest.raises(InputContractError, match="--te-assets"):
        cli._library_free(args)


def _assemble(tmp_path: Path, resident_filter, monkeypatch: pytest.MonkeyPatch):
    from evitea.library_free import assemble as AS
    from evitea.library_free import calling as CL
    from evitea.library_free import detector as D
    from evitea.library_free.options import Options
    from evitea.library_free.runtime import _assemble_and_write
    from evitea.shadow_layers import ShadowCallset

    seen = {}

    def assemble(first, second, index, opts, **kwargs):
        seen["index"] = index
        return [], {"families_with_similarity": 0}

    monkeypatch.setattr(D, "parse_places", lambda _b: None)
    monkeypatch.setattr(D, "parse_segments", lambda _b: None)
    monkeypatch.setattr(AS, "assemble", assemble)
    monkeypatch.setattr(AS, "added_calls_tsv_bytes", lambda _rows: ADDED)
    monkeypatch.setattr(AS, "family_calls_tsv_bytes", lambda _rows, extended: b"")
    tmp_path.mkdir()
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 250 + "\n", encoding="ascii")
    out = tmp_path / "lf"
    out.mkdir()
    summary = _assemble_and_write(
        out, (b"p", b"s"), None, Options.previous(), {}, native=NATIVE, native_rows=ShadowCallset.from_bytes(NATIVE).rows,
        sample_id="s", reference=SimpleNamespace(reference_fasta_path=fasta, contigs=(("chr1", 1000),)), genome=None, copies=None,
        genome_bp=0, rule_c=False, run_dir=tmp_path, command_line=(), source_identity=None, repeat_map=None, repeat_annotation=None,
        say=lambda _m: None, started=0.0, resident_filter=resident_filter)
    return out, summary, seen, CL


def test_library_free_reads_every_library_call_and_filters_only_the_final_tables(tmp_path: Path, own_filter: OwnFamilyCopyFilter,
                                                                              monkeypatch: pytest.MonkeyPatch) -> None:
    from evitea.shadow_layers import ShadowCallset

    off, off_summary, off_seen, CL = _assemble(tmp_path / "off", None, monkeypatch)
    on, on_summary, on_seen, _ = _assemble(tmp_path / "on", own_filter, monkeypatch)
    # the "no library call within 500 bp" index is built from the full table in both runs
    full = CL.CallIndex.from_native(ShadowCallset.from_bytes(NATIVE).rows)
    assert on_seen["index"].__dict__ == off_seen["index"].__dict__ == full.__dict__
    assert (on / "calls.library_free.tsv").read_bytes() == (off / "calls.library_free.tsv").read_bytes() == ADDED
    assert len((off / "calls.combined.tsv").read_text().splitlines()) == 7  # header, five library calls, one library-free call
    assert len((on / "calls.combined.tsv").read_text().splitlines()) == 5
    assert "resident_copy_dropped.tsv" in on_summary["output_hashes"] and "resident_copy_dropped.tsv" not in off_summary["output_hashes"]
    assert on_summary["resident_copy_filter"]["dropped"] == 2 and "resident_copy_filter" not in off_summary
    assert on_summary["source_native_callset_sha256"] == off_summary["source_native_callset_sha256"]
    assert not (off / "resident_copy_dropped.tsv").exists()
