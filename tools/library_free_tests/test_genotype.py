"""The genotype step (evitea.genotype, `evitea genotype`).   Run: PYTHONPATH=src python -B tools/library_free_tests/test_genotype.py   (venv python)

The step adds genotypes.tsv and genotype_summary.json to a library-free output and changes no call. EXPECTED_CODE_HASH (a prefix of the detection
code hash) makes H1 check that the detection code is the code of the build to keep calls identical with.

G1  Step B on a homozygous sample: every library call gets ref_share <= 0.15, a count (n_ref, n_total) and hom, method ref_walk.
G2  Step B on a half-and-half sample (reads from the insertion allele and the reference): ref_share in (0.15, 0.85], het.
G3  Step A on a run's library-free calls: ref_share equals calling.e_value of the call's ADDED row, method lf_walk, n_ref and n_total are `.`.
G4  Step A on crafted tables: the min of fracL and fracR when the right walk broke, fracL alone otherwise, hom / het / low at the thresholds of
    _p4_lib (0.15, 0.85), `.` for an empty share.
G5  no call changes: the run's output (every file) is byte for byte the same before and after the step, in the folder itself and in a folder of its own;
    the same run made twice (with and without the step) writes the same call files; the input hashes are recorded.
G6  the step does not overwrite genotypes.tsv or genotype_summary.json.
G7  without an msBWT index only Step A is made: the library calls are left out, step_b_run is false.
G8  the two library-free tables of different runs (different numbers of added calls, or a different bL) are refused.
G9  walk_at_break gives no count for an anchor at a contig edge, for an anchor with a letter other than A, C, G, T, and for a site where the
    insertion is not carried by any read (no divergence).
G10 the summary records the thresholds of _p4_lib (not tuned), the counter, the wall time of the walks and the genotype counts.
G11 plain and fast counter give the same genotypes.tsv.
C1  command line: `evitea genotype` exists; --msbwt-index and --sample-assets exclude each other.
H1  genotype.py is not among the files of the detection identity; the detection code hash is unchanged when EXPECTED_CODE_HASH is set.
"""
import hashlib
import json
import os
import random
import sys
import tempfile
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))
import synth_elements as SE  # noqa: E402
import evitea.backends.msbwt as M  # noqa: E402
from evitea import genotype as G  # noqa: E402
from evitea.errors import InputContractError  # noqa: E402
from evitea.library_free import _p4_lib as P  # noqa: E402
from evitea.library_free import calling as CL  # noqa: E402
from evitea.library_free import runtime as R  # noqa: E402

FAIL, NCHK = [], [0]


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail)[:600] if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def tree_hashes(folder: Path) -> dict:
    return {str(p.relative_to(folder)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.rglob("*")) if p.is_file()}


def tsv(path: Path) -> list[dict]:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


def make_run(tmp: Path, name: str, het: bool, seed: int = 1):
    """A library-free run on the small synthetic sample; returns (output folder, reference fasta, index folder)."""
    S = SE.build(seed=seed, ref_len=36_000, n_px=3, n_lx=2, n_uq=1, n_st=0)
    reads = S["reads"]
    if het:
        r = random.Random(seed + 1)
        reads = []
        for _ in range(len(S["reads"])):
            src = S["sample"] if r.random() < 0.5 else S["ref"]
            a = r.randrange(0, len(src) - 125)
            w = src[a:a + 125]
            reads.append(w if r.random() < 0.5 else P.rc(w))
    base = tmp / name
    base.mkdir()
    fa = base / "ref.fa"
    fa.write_text(">chrS\n" + S["ref"] + "\n")
    run = base / "run"
    run.mkdir()
    rows = "".join(f"N{k}\tchrS\t{s['pos']}\t{s['pos'] + 8}\tPX\tPX\n" for k, s in enumerate(x for x in S["sites"] if x["kind"] == "PX"))
    (run / "calls.native.tsv").write_text("event_id\tchr\tstart\tend\tfamily\tfamilies\n" + rows)
    (run / "run_summary.json").write_text(json.dumps({"sample_id": "S"}))
    counter = SE.SubCounter(reads)
    M.load_msbwt = lambda path: types.SimpleNamespace(countOccurrencesOfSeq=counter._occ)
    idx = base / "nobwt"
    idx.mkdir()
    R.MIN_DEFAULT_CONTIG = 1
    out = base / "out"
    R.run_library_free(run_dir=run, reference=types.SimpleNamespace(reference_fasta_path=fa), msbwt_index=idx, output=out,
                       work_dir=base / "work", contigs=("chrS",), rule_c=True, detection_dir=None, progress=lambda m: None,
                       het_guard="off", single_copy="none")
    return out, fa, idx


def main():
    tmp = Path(tempfile.mkdtemp(prefix="gt_test_"))
    hom_out, hom_fa, hom_idx = make_run(tmp, "hom", het=False)
    hom_before = tree_hashes(hom_out)

    # G1, G3, G5, G11 on the homozygous sample
    lib_msbwt = lambda path: types.SimpleNamespace(countOccurrencesOfSeq=SE.SubCounter(
        SE.build(seed=1, ref_len=36_000, n_px=3, n_lx=2, n_uq=1, n_st=0)["reads"])._occ)
    M.load_msbwt = lib_msbwt
    own = tmp / "gt_hom"
    summary = G.run_genotype(library_free_dir=hom_out, output=own, reference_fasta=hom_fa, msbwt_index=hom_idx, counter="plain")
    rows = tsv(own / "genotypes.tsv")
    lib = [r for r in rows if r["source"] == "EVITEA"]
    lf = [r for r in rows if r["source"] == "LIBRARY_FREE"]
    check("G1 three library calls, ref_walk, counts, hom",
          len(lib) == 3 and all(r["method"] == "ref_walk" and r["genotype"] == "hom" and float(r["ref_share"]) <= P.HOM_MAX
                                and r["n_total"] != "." and int(r["n_total"]) > 0 and r["n_ref"] != "." for r in lib), lib)
    family = tsv(hom_out / "library_free_family_calls.tsv")
    added = [r for r in family if r["status"] == "ADDED"]
    want = [CL.e_value(r["fracL"], r["fracR"], r["rstat"]) for r in added]
    check("G3 library-free rows carry e_value of the ADDED rows, lf_walk, no counts",
          len(lf) == len(added) >= 1 and all(abs(float(r["ref_share"]) - w) < 5e-4 for r, w in zip(lf, want))
          and all(r["method"] == "lf_walk" and r["n_ref"] == "." and r["n_total"] == "." for r in lf), (lf, want))
    check("G3b the table has one row per row of calls.combined.tsv, in its order",
          [r["event_id"] for r in rows] == [r["event_id"] for r in tsv(hom_out / "calls.combined.tsv")])

    # G5 no call changes
    check("G5a the output folder is unchanged by a step with a folder of its own", tree_hashes(hom_out) == hom_before)
    in_place = tmp / "hom_in_place"
    import shutil
    shutil.copytree(hom_out, in_place)
    G.run_genotype(library_free_dir=in_place, output=in_place, reference_fasta=hom_fa, msbwt_index=hom_idx)
    after = tree_hashes(in_place)
    new = set(after) - set(hom_before)
    check("G5b in place: only genotypes.tsv and genotype_summary.json are new, every other file is byte for byte the same",
          new == {"genotypes.tsv", "genotype_summary.json"} and all(after[k] == v for k, v in hom_before.items()), sorted(new))
    again, _, _ = make_run(tmp, "hom_again", het=False)
    check("G5c the same run made again writes the same call files (the step is not part of the run)",
          all(tree_hashes(again)[f] == hom_before[f] for f in ("calls.combined.tsv", "calls.library_free.tsv", "library_free_family_calls.tsv")))
    check("G5d the summary records the sha256 of the three input files, equal to the files",
          summary["input_sha256"] == {n: hashlib.sha256((hom_out / n).read_bytes()).hexdigest()
                                      for n in ("calls.combined.tsv", "calls.library_free.tsv", "library_free_family_calls.tsv")})
    fast = tmp / "gt_hom_fast"
    G.run_genotype(library_free_dir=hom_out, output=fast, reference_fasta=hom_fa, msbwt_index=hom_idx, counter="fast")
    check("G11 plain and fast counter give the same genotypes.tsv", (fast / "genotypes.tsv").read_bytes() == (own / "genotypes.tsv").read_bytes())

    # G6 no overwrite
    try:
        G.run_genotype(library_free_dir=hom_out, output=own, reference_fasta=hom_fa, msbwt_index=hom_idx)
        refused = False
    except InputContractError:
        refused = True
    check("G6 an existing genotypes.tsv is not overwritten", refused)

    # G7 no index
    only_a = tmp / "gt_only_a"
    s7 = G.run_genotype(library_free_dir=hom_out, output=only_a, reference_fasta=None, msbwt_index=None)
    rows7 = tsv(only_a / "genotypes.tsv")
    check("G7 without an msBWT index: only the library-free rows, step_b_run false",
          not s7["step_b_run"] and all(r["source"] == "LIBRARY_FREE" for r in rows7) and len(rows7) == len(lf)
          and s7["genotype_counts"].get("library_calls_left_out") == 3, s7["genotype_counts"])

    # G2 half-and-half sample
    het_out, het_fa, het_idx = make_run(tmp, "het", het=True)
    het_gt = tmp / "gt_het"
    G.run_genotype(library_free_dir=het_out, output=het_gt, reference_fasta=het_fa, msbwt_index=het_idx)
    hrows = [r for r in tsv(het_gt / "genotypes.tsv") if r["source"] == "EVITEA"]
    check("G2 half-and-half sample: library calls het with ref_share in (0.15, 0.85]",
          len(hrows) == 3 and all(r["genotype"] == "het" and P.HOM_MAX < float(r["ref_share"]) <= P.HET_MAX for r in hrows), hrows)

    # G4 crafted family-call tables
    def fam(status="ADDED", fracL="", fracR="", rstat="break", bL="100", chrom="chrS"):
        return {"status": status, "fracL": fracL, "fracR": fracR, "rstat": rstat, "bL": bL, "chrom": chrom}

    def added_row(i, bL="100"):
        return {"event_id": f"LF_chrS_{bL}_{i}", "chr": "chrS", "bL": bL}

    cases = [("0.1", "0.4", "break", "0.100", "hom"), ("0.4", "0.1", "break", "0.100", "hom"), ("0.5", "", "no_break", "0.500", "het"),
             ("0.5", "0.05", "dead", "0.500", "het"), ("0.9", "0.95", "break", "0.900", "low"), ("0.15", "0.2", "break", "0.150", "hom"),
             ("0.85", "0.86", "break", "0.850", "het"), ("0.851", "0.9", "break", "0.851", "low"), ("", "", "no_break", ".", ".")]
    got = G.step_a_rows([added_row(i) for i in range(len(cases))],
                        [fam(fracL=c[0], fracR=c[1], rstat=c[2]) for c in cases] + [fam(status="DROPPED_V2", fracL="0.0")])
    ok = [got[f"LF_chrS_100_{i}"] == (c[3], ".", ".", c[4]) for i, c in enumerate(cases)]
    check("G4 shares and labels on crafted tables", all(ok), [(i, c, got[f"LF_chrS_100_{i}"]) for i, c in enumerate(cases) if not ok[i]])

    # G8 mismatching tables
    def refused_pair(added_rows, family_rows):
        try:
            G.step_a_rows(added_rows, family_rows)
            return False
        except InputContractError:
            return True
    check("G8 different counts are refused", refused_pair([added_row(1), added_row(2)], [fam(fracL="0.1")]))
    check("G8b a different bL is refused", refused_pair([added_row(1, "100")], [fam(fracL="0.1", bL="101")]))

    # G9 walk_at_break edges, on the homozygous sample's reference and reads
    S = SE.build(seed=1, ref_len=36_000, n_px=3, n_lx=2, n_uq=1, n_st=0)
    C = P.Counter(types.SimpleNamespace(countOccurrencesOfSeq=SE.SubCounter(S["reads"])._occ))
    ref = S["ref"].upper()
    share, n_ref, n_total, why = G.walk_at_break(C, ref, 5, 13)
    check("G9a an anchor beyond the contig start gives no count", share is None and why == "edge+edge" or (share is None and "edge" in why), why)
    nref = ref[:20_000] + "N" * 60 + ref[20_060:]
    share, _, _, why = G.walk_at_break(P.Counter(types.SimpleNamespace(countOccurrencesOfSeq=SE.SubCounter(S["reads"])._occ)), nref, 20_030, 20_040)
    check("G9b an N run at the anchors gives no count", share is None, why)
    quiet = next(p for p in range(2_000, 34_000, 997) if all(abs(p - s["pos"]) > 400 for s in S["sites"]))
    share, _, _, why = G.walk_at_break(C, ref, quiet, quiet + 1)
    check("G9c a site where no read carries an insertion gives no count", share is None, (quiet, why))

    # G10 summary
    s = json.loads((own / "genotype_summary.json").read_text())
    check("G10 summary: thresholds from _p4_lib, counter, walk time, counts",
          s["thresholds"]["hom_max"] == P.HOM_MAX == 0.15 and s["thresholds"]["het_max"] == P.HET_MAX == 0.85 and s["counter"] == "plain"
          and "step_b_walks" in s["seconds"] and s["genotype_counts"].get("EVITEA:hom") == 3 and s["step_b_outcomes"] == {"break": 3}, s)

    # C1 command line
    from evitea.cli import build_parser
    parser = build_parser()
    args = parser.parse_args(["genotype", "--library-free-dir", "x", "--msbwt-index", "i", "--reference-assets", "r"])
    check("C1a evitea genotype parses", args.command == "genotype" and args.counter == "fast")
    try:
        parser.parse_args(["genotype", "--library-free-dir", "x", "--msbwt-index", "i", "--sample-assets", "s"])
        excl = False
    except SystemExit:
        excl = True
    check("C1b --msbwt-index and --sample-assets exclude each other", excl)

    # H1
    check("H1a genotype.py is not a detection code file", "genotype.py" not in set(R.CODE_FILES))
    h = R._code_sha256()
    want_hash = os.environ.get("EXPECTED_CODE_HASH")
    if want_hash:
        check("H1b detection code hash", h.startswith(want_hash), h)
    else:
        print("SKIP H1b (EXPECTED_CODE_HASH not set); hash", h[:12])
    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
