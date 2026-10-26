"""Run the library-layer engine end to end and compare it to recorded artifacts.

Validation-only.  Production code never imports this module.  It runs
``evitea.incumbent.run_incumbent_caller`` from raw inputs (sample msBWT, TE
library, family map, reference, host Bowtie2 index) and, when a recorded report
root is given, compares candidate loci and every emitted callset stage with
the recorded stage tables.  Truth and comparator callsets are never read.

Usage (Windows PowerShell launcher: tools/windows/VALIDATE_INCUMBENT.ps1):

    python validation/incumbent_stage_checks/run_incumbent_validation.py --work <work-dir>
    python validation/incumbent_stage_checks/run_incumbent_validation.py --sample ISO1 --work <work-dir>

Per-sample recorded files (``--sample`` sets the defaults):

* AKA017: msBWT ``AKA017/msbwt_SRR9951105``; stage tables
  ``AKA017/route_b_r9e_validation_full55x/shadow_v2/reports``; callset
  of 792 rows, SHA-256 ``1d1110f9...1757``.
* ISO1: msBWT ``TEMP2_ISO1/msbwt_SRR11906527``; stage tables
  ``TE_locator_route_b_v2/reports``; callset of 192 rows, SHA-256
  ``b3e48681...6ebc5``.

Both samples share the TE library, FBte family map, dm6 reference and its
prebuilt Bowtie2 index.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path


def _require_env(name):
    """Path from environment variable ``name``; exit with a clear message if unset."""
    value = os.environ.get(name)
    if not value:
        raise SystemExit(f"Set the {name} environment variable (see validation/README.md).")
    return Path(value)

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from evitea.incumbent import IncumbentInputs, MsbwtSample, legacy_tsv_bytes, run_incumbent_caller  # noqa: E402

DATA = _require_env("EVITEA_DATA")
ASSET = DATA / "TEMP2_ISO1"
ROUTE_B_V2 = _require_env("EVITEA_ROUTE_B_V2")
FROZEN_R9E_SHA256 = {
    "AKA017": "1d1110f939996223379edf74847d626d4bc397a5cac32ad16d5d985dca051757",
    "ISO1": "b3e48681568d616142fe534995d3eeb9e5163907992e9d0893919aea42e6ebc5",
}
FROZEN_R9E_ROWS = {"AKA017": 792, "ISO1": 192}
SAMPLE_DEFAULTS = {
    "AKA017": {
        "msbwt": DATA / "AKA017" / "msbwt_SRR9951105",
        "v2_reports": DATA / "AKA017" / "route_b_r9e_validation_full55x" / "shadow_v2" / "reports",
        "rescue_read_names": "pair",
    },
    "ISO1": {
        "msbwt": ASSET / "msbwt_SRR11906527",
        "v2_reports": ROUTE_B_V2 / "reports",
        # The recorded ISO1 wide-anchor host mapping names queries {pid}.
        "rescue_read_names": "numeric",
    },
}


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", "--sample-id", dest="sample_id", default="AKA017", choices=sorted(SAMPLE_DEFAULTS))
    ap.add_argument("--msbwt", type=Path, default=None, help="default: the sample's msBWT")
    ap.add_argument("--te-fasta", type=Path, default=ASSET / "TE consensus library" / "transposon_sequence_set.fa")
    ap.add_argument("--family-map", type=Path, default=ASSET / "fbte_to_flybase_family.tsv")
    ap.add_argument("--reference", type=Path, default=ASSET / "Reference" / "dm6_TEMP2_compatible.fa")
    ap.add_argument("--host-index", type=Path, default=ASSET / "Reference" / "bowtie2_index" / "dm6")
    ap.add_argument("--bowtie2", type=Path, default=ASSET / "tools" / "bowtie2" / "bowtie2.bat")
    ap.add_argument("--bowtie2-build", type=Path, default=ASSET / "tools" / "bowtie2" / "bowtie2-build.bat")
    ap.add_argument("--te-index", type=Path, default=None, help="prebuilt TE-library Bowtie2 prefix (built if omitted)")
    ap.add_argument("--v2-reports", type=Path, default=None,
                    help="recorded stage-report root to compare against (default: the sample's; 'none' skips)")
    ap.add_argument("--work", type=Path, required=True, help="work/checkpoint directory (resumable)")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args(argv)
    defaults = SAMPLE_DEFAULTS[args.sample_id]
    args.msbwt = args.msbwt or defaults["msbwt"]
    args.v2_reports = args.v2_reports or defaults["v2_reports"]
    return args


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main(argv=None) -> int:
    args = parse_args(argv)
    args.work.mkdir(parents=True, exist_ok=True)
    started = time.time()
    log = (args.work / "validation.log").open("a", encoding="utf-8")

    def progress(message: str) -> None:
        line = f"[{time.time() - started:9.1f}s] {message}"
        print(line, flush=True)
        log.write(line + "\n")
        log.flush()

    inputs = IncumbentInputs(
        sample_id=args.sample_id, te_fasta=args.te_fasta, te_family_map=args.family_map,
        reference_fasta=args.reference, host_index_prefix=args.host_index, bowtie2=args.bowtie2,
        bowtie2_build=args.bowtie2_build, te_index_prefix=args.te_index, threads=args.threads,
        workers=args.workers, rescue_read_names=SAMPLE_DEFAULTS[args.sample_id]["rescue_read_names"],
    )
    result = run_incumbent_caller(inputs, MsbwtSample(args.msbwt), args.work / "engine", progress=progress)
    final = result.legacy_tsv()
    (args.work / "INCUMBENT_COMBINED_CALLS.tsv").write_bytes(final)
    report: dict[str, object] = {
        "sample_id": args.sample_id,
        "elapsed_seconds": round(time.time() - started, 1),
        "counts": result.counts,
        "final_sha256": hashlib.sha256(final).hexdigest(),
        "checks": {},
    }
    checks: dict[str, bool] = report["checks"]  # type: ignore[assignment]
    if str(args.v2_reports).lower() != "none" and args.v2_reports.is_dir():
        v2 = args.v2_reports
        loci = [{
            "locus_id": L.locus_id, "chr": L.chrom, "candidate_family": L.family,
            "candidate_start_midpoint": f"{L.start:.1f}", "candidate_end_midpoint": f"{L.end:.1f}",
            "candidate_midpoint": f"{L.midpoint:.1f}", "candidate_pair_count": str(len(L.pair_ids)),
        } for L in result.acquisition.loci]
        frozen = _tsv(v2 / "V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1" / "V2FCRC1_CANDIDATE_LOCI.tsv")
        checks["candidate_loci"] = loci == [{key: row[key] for key in loci[0]} for row in frozen] if loci else not frozen
        fc = v2 / "V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1"
        checks["split_junction_all_events"] = legacy_tsv_bytes(result.split_junction.all_events) == (fc / "V2FCRC1_ALL_EVENTS.tsv").read_bytes()
        checks["split_junction_calls"] = legacy_tsv_bytes(result.split_junction.calls) == (fc / "V2FCRC1_CALLS_SCORE4.tsv").read_bytes()
        r4 = v2 / "V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4" / "V2R4_COMBINED_CALLS.tsv"
        checks["mate_rescue_combined"] = legacy_tsv_bytes((*result.split_junction.calls, *result.mate_rescue.calls)) == r4.read_bytes()
        r9e = v2 / "V2_R9E_REPEAT_AWARE_DIVERGENCE" / "R9E_COMBINED_CALLS.tsv"
        checks["repeat_aware_final"] = final == r9e.read_bytes()
    expected = FROZEN_R9E_SHA256.get(args.sample_id)
    if expected:
        checks["frozen_final_sha256"] = report["final_sha256"] == expected
        checks["frozen_final_rows"] = len(result.calls) == FROZEN_R9E_ROWS[args.sample_id]
    report["verdict"] = "INCUMBENT_EQUIVALENCE_PASS" if checks and all(checks.values()) else "INCUMBENT_EQUIVALENCE_NOT_SHOWN"
    (args.work / "INCUMBENT_VALIDATION_REPORT.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    progress(json.dumps({"verdict": report["verdict"], "checks": checks, "final_sha256": report["final_sha256"]}))
    return 0 if report["verdict"] == "INCUMBENT_EQUIVALENCE_PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
