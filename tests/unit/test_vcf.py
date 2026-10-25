from __future__ import annotations

from pathlib import Path

from evitea.callset import CallsetState, calls_bed_bytes, calls_tsv_bytes, calls_vcf_bytes
from evitea.vcf import calls_vcf_bytes as vcf_bytes
from evitea.vcf import library_free_vcf_calls, native_vcf_calls, reference_bases

NATIVE = (
    "event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\tleft_support\tright_support\t"
    "total_support\tevidence_type\tscore\tfragment_ids\n"
    "V2FCRC1_EVENT_000002\tchr2\t40\t46\tL1MdA\t['L1MdA']\t46\t40\t-6\t5\t5\t10\ttwo_sided\t12.0\t1,2,3\n"
    "V2FCRC1_EVENT_000001\tchr1\t10\t15\tB1\t['B1', 'B2']\t15\t10\t-5\t8\t5\t13\ttwo_sided\t15.0\t4,5\n"
)
ADDED = (
    "event_id\tchr\tstart\tend\tfamily\tcluster\tbL\tbR\ttsd\te\tc\n"
    "LF_chr1_30\tchr1\t30\t31\tLF_cluster_7\t7\t30\t36\t6\t0.01\t\n"
)
CONTIGS = (("chr1", 60), ("chr2", 80))
REFERENCE = ">chr1\n" + "ACGTACGTAC" * 3 + "\n" + "TTTTTGGGGG" * 3 + "\n>chr2\n" + "C" * 39 + "G" + "A" * 40 + "\n"


def _records(text: bytes) -> list[list[str]]:
    return [line.split("\t") for line in text.decode().splitlines() if not line.startswith("#")]


def test_native_calls_become_symbolic_insertions(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fa"
    fasta.write_text(REFERENCE, encoding="ascii")
    state = CallsetState(sample_id="s1", native=NATIVE.encode())
    tsv, bed = calls_tsv_bytes(state), calls_bed_bytes(state, dict(CONTIGS))
    text = calls_vcf_bytes(state, CONTIGS, fasta)
    header = [line for line in text.decode().splitlines() if line.startswith("#")]
    assert header[0] == "##fileformat=VCFv4.2"
    assert "##contig=<ID=chr1,length=60>" in header and "##reference=ref.fa" in header
    assert any(line.startswith("##ALT=<ID=INS:ME:L1MdA,") for line in header)
    assert header[-1] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    first, second = _records(text)
    # reference order, then position: chr1 before chr2
    assert first[:7] == ["chr1", "10", "V2FCRC1_EVENT_000001", "C", "<INS:ME:B1>", ".", "PASS"]
    info = dict(item.split("=", 1) if "=" in item else (item, True) for item in first[7].split(";"))
    assert info["SVTYPE"] == "INS" and info["END"] == "15" and info["CIPOS"] == "0,5" and info["IMPRECISE"] is True
    assert info["FAMILY"] == "B1" and info["FAMILIES"] == "B1,B2" and info["SOURCE"] == "library"
    assert info["ORIGIN"] == "INCUMBENT_SPLIT_JUNCTION" and info["SUPPORT"] == "13" and info["LEFT_SUPPORT"] == "8"
    assert second[:5] == ["chr2", "40", "V2FCRC1_EVENT_000002", "G", "<INS:ME:L1MdA>"]
    assert "FAMILIES" not in second[7]
    # the TSV and BED views do not change by writing the VCF
    assert (calls_tsv_bytes(state), calls_bed_bytes(state, dict(CONTIGS))) == (tsv, bed)
    assert calls_vcf_bytes(state, CONTIGS, fasta) == text


def test_library_free_calls_carry_their_source_and_tsd(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fa"
    fasta.write_text(REFERENCE, encoding="ascii")
    state = CallsetState(sample_id="s1", native=NATIVE.encode())
    from evitea.callset import call_views

    calls = [*native_vcf_calls(call_views(state)), *library_free_vcf_calls(ADDED.encode())]
    text = vcf_bytes(calls, CONTIGS, sample_id="s1", bases=reference_bases(fasta, calls))
    rows = _records(text)
    assert [r[2] for r in rows] == ["V2FCRC1_EVENT_000001", "LF_chr1_30", "V2FCRC1_EVENT_000002"]
    lf = rows[1]
    assert lf[:5] == ["chr1", "30", "LF_chr1_30", "C", "<INS:ME:LF_cluster_7>"]
    assert "SOURCE=library-free" in lf[7] and "ORIGIN=LIBRARY_FREE" in lf[7] and "TSD=6" in lf[7]
    # without the reference the REF base is N
    assert _records(vcf_bytes(calls, CONTIGS, sample_id="s1"))[0][3] == "N"
