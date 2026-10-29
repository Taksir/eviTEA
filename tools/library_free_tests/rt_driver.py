"""Run evitea.library_free.runtime.run_library_free on the synthetic sample with the package on PYTHONPATH (previous or current).
usage: rt_driver.py <out dir> <seed> <options json, merge or -> [--het] [--annotation file] [--repeat-map file] [--native none|px]
                    [--detection-dir dir] [--repeat-anchor windows,at_cutoff] [--profile-min letters] [--sample elements|cross]
                    [--single-copy auto|all|none|LIST] [--het-guard on|off]
                    [--tsd-skip on|off]
merge: Options.merge(). The cross-species repeat cutoff needs the anchor reference's window counts; the driver sets them to a share of
1 (the cutoff stays REP_MAX) unless --repeat-anchor gives others, and --profile-min lowers the 1 Mb contig minimum of the profile.
--sample cross: the sample of synth_cross.py (one scenario per cross-species rule) instead of synth_elements.py."""
import json
import os
import random
import sys
import types
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import synth_elements as SE  # noqa: E402
from evitea.library_free import _p4_lib as P  # noqa: E402
import evitea.backends.msbwt as M  # noqa: E402
from evitea.library_free import runtime as R  # noqa: E402

out = Path(sys.argv[1])
seed = int(sys.argv[2])
opt_json = sys.argv[3]
het = "--het" in sys.argv
annotation = sys.argv[sys.argv.index("--annotation") + 1] if "--annotation" in sys.argv else None
repeat_map = sys.argv[sys.argv.index("--repeat-map") + 1] if "--repeat-map" in sys.argv else None
detection_dir = sys.argv[sys.argv.index("--detection-dir") + 1] if "--detection-dir" in sys.argv else None
native_mode = sys.argv[sys.argv.index("--native") + 1] if "--native" in sys.argv else "none"

sample_kind = sys.argv[sys.argv.index("--sample") + 1] if "--sample" in sys.argv else "elements"
if sample_kind == "cross":
    import synth_cross as SX  # noqa: E402
    S = SX.build(seed=seed)
else:
    S = SE.build(seed=seed, ref_len=80_000, n_px=6, n_lx=3, n_uq=2, n_st=1)
ref = S["ref"]
if het and sample_kind == "cross":
    raise SystemExit("--het is for the elements sample")
if het:
    r = random.Random(seed + 1)
    reads = []
    for _ in range(len(S["reads"])):
        src = S["sample"] if r.random() < 0.5 else S["ref"]
        a = r.randrange(0, len(src) - 125)
        w = src[a:a + 125]
        reads.append(w if r.random() < 0.5 else P.rc(w))
else:
    reads = S["reads"]
work = out.parent / (out.name + "_in")
work.mkdir(parents=True, exist_ok=True)
fa = work / "ref.fa"
fa.write_text(">chrS\n" + ref + "\n")
run_dir = work / "run"
run_dir.mkdir(exist_ok=True)
head = "event_id\tchr\tstart\tend\tfamily\tfamilies\n"
rows = ""
if native_mode == "px":      # library-confirmed calls at the P-element-like sites
    for k, s in enumerate(x for x in S["sites"] if x["kind"] == "PX"):
        rows += f"N{k}\tchrS\t{s['pos']}\t{s['pos'] + 8}\tPX\tPX\n"
(run_dir / "calls.native.tsv").write_text(head + rows)
(run_dir / "run_summary.json").write_text(json.dumps({"sample_id": "S"}))
C = SE.SubCounter(reads)
M.load_msbwt = lambda path: types.SimpleNamespace(countOccurrencesOfSeq=C._occ)
reference = types.SimpleNamespace(reference_fasta_path=fa)
idx = work / "nobwt"
idx.mkdir(exist_ok=True)
kw = {}
if opt_json == "merge":
    from evitea.library_free.options import Options
    kw["options"] = Options.merge()
elif opt_json != "-":
    from evitea.library_free.options import Options
    kw["options"] = Options(**json.loads(opt_json))
try:
    from evitea.library_free import cross as X
except ImportError:      # a version without the cross-species rules
    X = None
if X is not None:
    a = sys.argv[sys.argv.index("--repeat-anchor") + 1].split(",") if "--repeat-anchor" in sys.argv else ("1", "1")
    X.REPEAT_ANCHOR = {"windows": int(a[0]), "at_cutoff": int(a[1])}
    if "--profile-min" in sys.argv:
        X.MIN_PROFILE_CONTIG = int(sys.argv[sys.argv.index("--profile-min") + 1])
if annotation:
    kw["repeat_annotation"] = Path(annotation)
if repeat_map:
    kw["repeat_map"] = Path(repeat_map)
if "--single-copy" in sys.argv:
    kw["single_copy"] = sys.argv[sys.argv.index("--single-copy") + 1]
if "--het-guard" in sys.argv:
    kw["het_guard"] = sys.argv[sys.argv.index("--het-guard") + 1]
if "--tsd-skip" in sys.argv:
    from evitea.library_free.tsd_skip import Skip
    kw["tsd_skip"] = Skip(on=sys.argv[sys.argv.index("--tsd-skip") + 1] == "on")
R.MIN_DEFAULT_CONTIG = 1
R.run_library_free(run_dir=run_dir, reference=reference, msbwt_index=idx, output=out, work_dir=out.parent / (out.name + "_work"),
                   contigs=("chrS",), rule_c=True, detection_dir=Path(detection_dir) if detection_dir else None,
                   progress=lambda m: None, **kw)
print("done")
