"""Opt-in depth-aware floors of the library-free layer (--extension-stop-min, --extension-stop-min-depth-cut).
Run: python -B tools/library_free_tests/test_depth_floors.py   (venv python)

V1  default switches write the freeze's files byte for byte (FREEZE_SRC = the src folder of the freeze, 296e725; skipped without it):
    every output file is equal (compare_scan.py outputs mode: the summary after the timings and the source identity), and the
    current code reuses the freeze's saved scan pieces with EVITEA_PIECES_READONLY=1 (a changed piece identity would refuse the run).
V2  --extension-stop-min 3 without a cut is the default (no extension_floor entry, files equal to the default run).
V3  a floor other than 3 is recorded in the summary (requested, used, no depth) and enters the identity of the saved pieces: the run
    does not reuse default pieces (EVITEA_PIECES_READONLY=1 refuses it).
V4  the depth cut: the floor is N at or under the cut and 3 above it; the single-copy depth used is written to the summary; above the
    cut the files equal the default run's.
V5  --detection-dir with a floor switch is refused.
V6  in-process: the floor is set for the scan and restored after, also after an error; bad values are refused; the depth cut picks the
    floor from the depth of the layer's own single-copy profile and refuses a sample whose depth could not be measured.
V7  the extension itself: on a path read by 2 reads, floor 3 stops and floor 2 goes on.
V8  the weak-site pass of the library-free step (--weak-site-pass, --weak-site-allele-depth-cut), on constructed stage tables: a run
    with the pass reuses the saved pieces of the run without it (EVITEA_PIECES_READONLY=1) and changes only the combined files and the
    summary; a site no combined row covers becomes a WEAK_ row (source EVITEA, also in the VCF), a site near a library row or near an
    added library-free call does not; above the cut every file equals the run without the pass; a missing stage table stops the run.
V9  the floor-1 add-on (--floor1-addon): at or under the cut a second step runs into floor1_addon/ with the extension floor at 1, and the
    main combined table is the run without the add-on plus the rule's rows (F1_ ids, at the end), with calls.library_free.tsv, the VCF
    and the summary (floor1_addon, current output hashes) to match; over the cut there is no second step and every file equals the run
    without it; without the option nothing is written; the option without --weak-site-pass is refused.
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
FREEZE = os.environ.get("FREEZE_SRC")
FAIL, NCHK = [], [0]


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def run(out, src, *extra, env_extra=None, small=True):
    env = dict(os.environ, PYTHONPATH=str(src))
    env.update(env_extra or {})
    return subprocess.run([sys.executable, "-B", str(HERE / "rt_driver.py"), str(out), "7", "-", *(["--small"] if small else []), *extra], env=env,
                          capture_output=True, text=True)


def same_outputs(a, b):
    r = subprocess.run([sys.executable, str(ROOT / "tools" / "compare_scan.py"), "outputs", str(a), str(b)], capture_output=True, text=True)
    return r.returncode == 0, r.stdout[-400:]


def same_tables(a, b):
    """Every file but the summary (which names the floor) equal byte for byte."""
    names = sorted(x.name for x in Path(a).iterdir() if x.is_file() and x.name != "library_free_summary.json")
    return names == sorted(x.name for x in Path(b).iterdir() if x.is_file() and x.name != "library_free_summary.json") and all(
        hashlib.sha256((Path(a) / n).read_bytes()).digest() == hashlib.sha256((Path(b) / n).read_bytes()).digest() for n in names), names


def summary(d):
    return json.loads((Path(d) / "library_free_summary.json").read_text())


def main():
    t = Path(tempfile.mkdtemp(prefix="df_"))
    src = ROOT / "src"
    base = t / "base"
    r = run(base, src)
    check("base run finishes", r.returncode == 0, r.stderr[-300:])
    check("base summary has no extension_floor", "extension_floor" not in summary(base)["detection"], list(summary(base)["detection"]))

    if FREEZE:
        fz = t / "fz"
        shared = ["--index-dir", str(t / "idx"), "--work-dir", str(t / "fz_work")]    # the index path is part of a piece's identity
        a = run(fz, FREEZE, *shared)
        check("V1 freeze run finishes", a.returncode == 0, a.stderr[-300:])
        ok, text = same_outputs(fz, base)
        check("V1 default switches: every output equals the freeze's", ok, text)
        reuse = t / "reuse"
        b = run(reuse, src, *shared, env_extra={"EVITEA_PIECES_READONLY": "1"})
        check("V1 current code reuses the freeze's saved pieces (identity unchanged)", b.returncode == 0, b.stderr[-300:])
        ok, text = same_outputs(fz, reuse)
        check("V1 and gives the freeze's files", ok, text)
    else:
        print("V1 skipped: FREEZE_SRC is not set")

    v2 = t / "v2"
    r = run(v2, src, "--extension-stop-min", "3")
    ok, text = same_outputs(base, v2)
    check("V2 explicit 3 is the default", r.returncode == 0 and ok and "extension_floor" not in summary(v2)["detection"], text)

    v3 = t / "v3"
    r = run(v3, src, "--extension-stop-min", "2")
    check("V3 floor 2 finishes", r.returncode == 0, r.stderr[-300:])
    rec = summary(v3)["detection"].get("extension_floor")
    check("V3 summary records the floor", rec == {"requested": 2, "depth_cut": None, "single_copy_depth": None, "used": 2, "default": 3}, rec)
    r = run(t / "v3b", src, "--extension-stop-min", "2", "--work-dir", str(t / "base_work"), env_extra={"EVITEA_PIECES_READONLY": "1"})
    check("V3 the floor is part of the piece identity (default pieces are not reused read-only)", r.returncode != 0, r.stderr[-200:])

    hi, lo = t / "hi", t / "lo"
    r1 = run(hi, src, "--extension-stop-min", "2", "--extension-stop-min-depth-cut", "1000000")
    r2 = run(lo, src, "--extension-stop-min", "2", "--extension-stop-min-depth-cut", "0.001")
    check("V4 both cut runs finish", r1.returncode == 0 and r2.returncode == 0, (r1.stderr[-200:], r2.stderr[-200:]))
    ra, rb = summary(hi)["detection"]["extension_floor"], summary(lo)["detection"]["extension_floor"]
    check("V4 depth at or under the cut: floor N", ra["used"] == 2 and ra["single_copy_depth"] and ra["depth_cut"] == 1000000.0, ra)
    check("V4 depth above the cut: floor 3, depth written", rb["used"] == 3 and rb["single_copy_depth"] == ra["single_copy_depth"]
          and rb["requested"] == 2, rb)
    ok, text = same_tables(base, lo)
    check("V4 above the cut the tables equal the default run's", ok, text)
    ok, text = same_tables(v3, hi)
    check("V4 at or under the cut the tables equal the run with the floor and no cut", ok, text)

    r = run(t / "v5", src, "--extension-stop-min", "2", "--detection-dir", str(base))
    check("V5 --detection-dir with a floor switch is refused", r.returncode != 0 and "--detection-dir" in r.stderr, r.stderr[-200:])

    sys.path.insert(0, str(ROOT / "src"))
    from evitea.errors import InputContractError
    from evitea.library_free import _p4_lib as P
    from evitea.library_free import runtime as R

    check("V6 default STOP_MIN is 3", P.STOP_MIN == 3 == R.STOP_MIN_DEFAULT)
    with R._stop_min(1):
        inside = P.STOP_MIN
    check("V6 floor set inside and restored after", inside == 1 and P.STOP_MIN == 3, (inside, P.STOP_MIN))
    try:
        with R._stop_min(1):
            raise RuntimeError("x")
    except RuntimeError:
        pass
    check("V6 restored after an error", P.STOP_MIN == 3)
    say = lambda m: None
    check("V6 default: no record, no depth measured", R._extension_floor(3, None, {}, Path("x"), Path("y"), None, t, say) == (3, None))
    check("V6 floor without cut: no depth measured",
          R._extension_floor(2, None, {}, Path("x"), Path("y"), None, t, say) ==
          (2, {"requested": 2, "depth_cut": None, "single_copy_depth": None, "used": 2, "default": 3}))
    for bad in ((0, None), (True, None), (2.5, None), (2, 0), (2, -1.0)):
        try:
            R._extension_floor(bad[0], bad[1], {}, Path("x"), Path("y"), None, t, say)
            check(f"V6 refused {bad}", False)
        except InputContractError:
            check(f"V6 refused {bad}", True)
    from evitea.library_free import single_copy as SC
    real = R._single_copy
    try:
        for depth, cut, want in ((10.0, 10.0, 2), (10.5, 10.0, 3), (4.0, 10.0, 2)):
            R._single_copy = lambda *a, **k: SC.Decision(SC.Choice("auto"), profile={"genome_depth": depth})
            got = R._extension_floor(2, cut, {}, Path("x"), Path("y"), None, t, say)
            check(f"V6 depth {depth} cut {cut}: floor {want}", got[0] == want and got[1]["single_copy_depth"] == depth, got)
        R._single_copy = lambda *a, **k: SC.Decision(SC.Choice("auto"), reason="no msBWT index")
        try:
            R._extension_floor(2, 10.0, {}, Path("x"), Path("y"), None, t, say)
            check("V6 unmeasurable depth is refused", False)
        except InputContractError as exc:
            check("V6 unmeasurable depth is refused", "no msBWT index" in str(exc), exc)
    finally:
        R._single_copy = real
    # V8: the weak-site pass of the library-free step, on constructed stage tables (elements sample with library calls at 6 sites)
    import hashlib
    import pickle
    sys.path.insert(0, str(ROOT / "src"))
    from evitea.incumbent.mate_rescue import MateRescueLocus, MateRescueResult
    from evitea.incumbent.split_junction import SplitJunctionResult

    def locus(pos, ids):
        return MateRescueLocus("L", "chrS", "PX", float(pos), float(pos + 8), float(pos + 4), tuple(ids[:1]), tuple(ids[1:2]), False)

    callwork = t / "callwork"
    folder = callwork / "engine" / "work" / "incumbent" / "checkpoints"
    folder.mkdir(parents=True)
    tables = (SplitJunctionResult((), (), (), {}),
              MateRescueResult((), (locus(35000, [1, 2]), locus(10300, [3, 4]), locus(6400, [5, 6])), (), {}),
              types.SimpleNamespace(ledger=()))
    for stage, table in zip(("split_junction", "mate_rescue", "junction_groups"), tables):
        payload = pickle.dumps(table)
        (folder / f"{stage}.pkl").write_bytes(payload)
        (folder / f"{stage}.complete.json").write_text(json.dumps({"fingerprint": {}, "sha256": hashlib.sha256(payload).hexdigest()}))
    shared = ["--native", "px", "--index-dir", str(t / "idx8"), "--work-dir", str(t / "w8_work")]
    b8 = t / "b8"
    r = run(b8, src, *shared, small=False)
    check("V8 run without the pass finishes", r.returncode == 0, r.stderr[-300:])
    rows = (b8 / "calls.combined.tsv").read_text().splitlines()
    check("V8 the run has library rows and added library-free calls", sum("\tEVITEA" in x for x in rows) == 6 and
          sum("LIBRARY_FREE" in x for x in rows) == 2, rows)
    ro = {"EVITEA_PIECES_READONLY": "1"}
    on = t / "w8"
    r = run(on, src, *shared, "--weak-site-allele-depth-cut", "1000000", "--weak-stage-dir", str(callwork), env_extra=ro, small=False)
    check("V8 the pass reuses the saved pieces (identity unchanged) and finishes", r.returncode == 0, r.stderr[-300:])
    new = (on / "calls.combined.tsv").read_text().splitlines()
    added = [x for x in new if x not in rows]
    check("V8 exactly one new row: the site no row covers, source EVITEA", len(added) == 1 and added[0].startswith("WEAK_000001\tchrS\t35000\t")
          and added[0].endswith("\tPX\tEVITEA"), added)
    check("V8 library rows first, then the weak row, then the library-free calls",
          [x.split("\t")[5] for x in new[1:]] == ["EVITEA"] * 7 + ["LIBRARY_FREE"] * 2, [x.split("\t")[5] for x in new[1:]])
    same = all(hashlib.sha256((b8 / n).read_bytes()).digest() == hashlib.sha256((on / n).read_bytes()).digest()
               for n in ("library_free_places.tsv", "library_free_segments.fa", "library_free_family_calls.tsv", "calls.library_free.tsv"))
    check("V8 places, segments, family calls and added calls are unchanged", same)
    vcf = (on / "calls.combined.vcf").read_text()
    check("V8 the VCF has the weak row as a library call (the synthetic native table has no evidence_type column)",
          "WEAK_000001" in vcf and "SOURCE=library" in vcf.split("WEAK_000001")[1].split("\n")[0]
          and "WEAK_000001" not in (b8 / "calls.combined.vcf").read_text())
    info = summary(on)["weak_site_pass"]
    check("V8 the summary records depth, e, a, Da, cut and the counts",
          info["ran"] and info["cut"] == 1000000.0 and info["single_copy_depth"] > 0 and info["a"] == 1.0 and info["allele_depth"] == info["single_copy_depth"]
          and info["e_source"] == "added library-free calls" and info["added_calls_used"] == 2 and (info["sites"], info["uncovered_sites"], info["passing_rule"], info["dropped_by_copy_filter"],
                                                 info["rows_added"], info["cover_rows"]) == (3, 1, 1, 0, 1, 8), info)
    off = t / "w8off"
    r = run(off, src, *shared, "--weak-site-allele-depth-cut", "0.001", "--weak-stage-dir", str(callwork), env_extra=ro, small=False)
    ok, names = same_tables(b8, off)
    check("V8 above the cut every file equals the run without the pass", r.returncode == 0 and ok, r.stderr[-200:])
    check("V8 above the cut the summary says it did not run", summary(off)["weak_site_pass"]["ran"] is False and "rows_added" not in summary(off)["weak_site_pass"])
    check("V8 without the option the summary has no weak_site_pass entry", "weak_site_pass" not in summary(b8))
    # V9: the floor-1 add-on (--floor1-addon): a second step with the extension floor at 1 when the allele depth is at or under the cut
    from evitea.library_free import floor1_addon as F1

    a9 = t / "a9"
    r = run(a9, src, *shared, "--weak-site-allele-depth-cut", "1000000", "--weak-stage-dir", str(callwork), "--floor1-addon", small=False)
    check("V9 the add-on run finishes", r.returncode == 0, r.stderr[-300:])
    s9 = summary(a9)["floor1_addon"]
    addon_dir = a9 / "floor1_addon"
    f1_rows = F1.read_rows((addon_dir / "calls.combined.tsv").read_bytes())
    main_rows = F1.read_rows((on / "calls.combined.tsv").read_bytes())          # the main step without the add-on (same weak pass)
    want = F1.added_rows(main_rows, f1_rows)
    check("V9 the add-on step wrote its own output folder with a summary at floor 1", (addon_dir / "library_free_summary.json").is_file()
          and summary(addon_dir)["detection"]["extension_floor"]["used"] == 1 and "floor1_addon" not in summary(addon_dir)
          and "weak_site_pass" not in summary(addon_dir), summary(addon_dir)["detection"].get("extension_floor"))
    got = F1.read_rows((a9 / "calls.combined.tsv").read_bytes())
    check("V9 the main combined table is the run without the add-on plus the add-on rows by the rule (F1_ ids, at the end)",
          got[:len(main_rows)] == main_rows and got[len(main_rows):] == [{**w} for w in want]
          and all(w["event_id"].startswith("F1_") for w in want), (len(main_rows), len(got), len(want)))
    check("V9 the summary records the switch", s9["requested"] is True and s9["ran"] is True and s9["rows_added"] == len(want)
          and s9["addon_library_free_rows"] == sum(1 for x in f1_rows if x["source"] == "LIBRARY_FREE") and s9["cut"] == 1000000.0
          and s9["allele_depth"] == summary(a9)["weak_site_pass"]["allele_depth"] and s9["addon_output"] == str(addon_dir)
          and s9["rule"] == F1.RULE, s9)
    added_lf = F1.read_rows((a9 / "calls.library_free.tsv").read_bytes())
    check("V9 calls.library_free.tsv has the same rows added, and the output hashes in the summary are current",
          len(added_lf) == len(F1.read_rows((on / "calls.library_free.tsv").read_bytes())) + len(want)
          and all(summary(a9)["output_hashes"][n] == hashlib.sha256((a9 / n).read_bytes()).hexdigest()
                  for n in ("calls.combined.tsv", "calls.library_free.tsv", "calls.combined.vcf")))
    check("V9 the VCF has one record per call of the combined table",
          sum(1 for x in (a9 / "calls.combined.vcf").read_text().splitlines() if not x.startswith("#")) == len(got))
    same9 = all(hashlib.sha256((on / n).read_bytes()).digest() == hashlib.sha256((a9 / n).read_bytes()).digest()
                for n in ("library_free_places.tsv", "library_free_segments.fa", "library_free_family_calls.tsv"))
    check("V9 places, segments and family calls of the main step are unchanged", same9)
    o9 = t / "o9"
    r = run(o9, src, *shared, "--weak-site-allele-depth-cut", "0.001", "--weak-stage-dir", str(callwork), "--floor1-addon", small=False)
    ok, names = same_tables(off, o9)
    s9 = summary(o9)["floor1_addon"]
    check("V9 allele depth over the cut: no second step, ran false, every file equals the run without the add-on",
          r.returncode == 0 and ok and not (o9 / "floor1_addon").exists() and s9["ran"] is False and s9["rows_added"] == 0
          and s9["addon_output"] is None and s9["allele_depth"] > s9["cut"], (r.stderr[-200:], s9))
    check("V9 without the option the summary has no floor1_addon entry and no floor1_addon folder",
          "floor1_addon" not in summary(on) and not (on / "floor1_addon").exists())
    r = run(t / "o9b", src, *shared, "--floor1-addon", small=False)
    check("V9 --floor1-addon without --weak-site-pass is refused", r.returncode != 0 and "--weak-site-pass" in r.stderr, r.stderr[-200:])
    (folder / "mate_rescue.pkl").unlink()
    r = run(t / "w8miss", src, *shared, "--weak-site-allele-depth-cut", "10", "--weak-stage-dir", str(callwork), env_extra=ro, small=False)
    check("V8 a missing stage table stops the run with a clear error", r.returncode != 0 and "mate_rescue" in r.stderr, r.stderr[-200:])

    # the extension itself: a path read by 2 reads goes on with the floor 2 and stops with 3 (n_alt 10 keeps the 10% term at 1)
    K = P.K
    path = "ACGTTGCAAGGCTTACCGATCGGATTCAGCTAAGCTTGGAC" * 3

    class Path2:
        def cnt(self, w):
            return 2 if any(path[i:i + K] == w for i in range(len(path) - K + 1)) else 0

    ctx = path[:K - 1]
    long, why3 = P.extend(Path2(), ctx[:-3], ctx[-3:], 10)
    with R._stop_min(2):
        longer, why2 = P.extend(Path2(), ctx[:-3], ctx[-3:], 10)
    check("V7 floor 3 stops on a 2-read path, floor 2 goes on", len(longer) > len(long) and why3 == "support", (len(long), len(longer), why3, why2))
    print(f"{NCHK[0]} checks, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
