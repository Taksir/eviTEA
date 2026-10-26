"""Validation-only: run the extension engine end to end, per sample.

``evitea.extension_engine.run_extension`` runs from each prepared msBWT (plus
the recorded safe terminal and internal-only codebooks, TE library, family map,
RepeatMasker BED and host Bowtie2 index) and adds the four layers to
the library-layer bytes.  Required, per sample:

* paired-end anchor candidates equal to the recorded tables;
* full-memory proposals (full data + ten rotations), event table and
  selected events equal to the recorded tables;
* internal 21-mer exact-junction-evidence witnesses equal to the recorded
  witness table, and the incremental event instances and locus table equal to
  the recorded observability outputs;
* layers byte-identical: AKA017 801 / 820 / 835 / 838, ISO1 192 / 193 / 194 / 194.

Host read-name profiles: AKA017 ``aka-historical`` (reproduction only),
ISO1 ``pid`` (default).  Resumable.
No truth, comparator or other-sample file is read.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import os
from pathlib import Path


def _env_path(name):
    """Path from environment variable ``name``; None when it is not set."""
    value = os.environ.get(name)
    return Path(value) if value else None

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_shadow_layer_equivalence as frozen  # noqa: E402

from evitea.anchor_pe import AnchorPEInputs  # noqa: E402
from evitea.extension_engine import run_extension  # noqa: E402
from evitea.full_memory.strict_a import EVENT_LEDGER_FIELDS  # noqa: E402
from evitea.incumbent import MsbwtSample  # noqa: E402
from evitea.internal21.observability import INSTANCE_FIELDS  # noqa: E402

SAMPLES = {
    "AKA017": {"msbwt": "AKA017/msbwt_SRR9951105", "profile": "aka-historical"},
    "ISO1": {"msbwt": "TEMP2_ISO1/msbwt_SRR11906527", "profile": "pid"},
}
LAYER_NAMES = ("bilateral", "asymmetric", "full-memory", "INTERNAL21")


def _text(value) -> str:
    return ";".join(map(str, value)) if isinstance(value, (list, tuple, set)) else str(value)


def _rows(path: Path, sample: str, extra=lambda r: True) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return [r for r in csv.DictReader(handle, delimiter="\t") if r.get("sample", sample) == sample and extra(r)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=_env_path("EVITEA_DATA"), required=_env_path("EVITEA_DATA") is None)
    ap.add_argument("--route-b-v2", type=Path, required=True)
    ap.add_argument("--route-b-v3", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--bowtie2", type=Path, default=None)
    ap.add_argument("--edlib-python", type=Path, default=None,
                    help="interpreter with edlib, used only if this one lacks it")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--samples", nargs="+", default=list(SAMPLES), choices=list(SAMPLES))
    args = ap.parse_args(argv)

    shared = args.data / "TEMP2_ISO1"
    reports = args.route_b_v3 / "reports"
    bowtie2 = args.bowtie2 or shared / "tools" / "bowtie2" / "bowtie2.bat"
    codebook = reports / "V3_FINAL_SAFE_TERMINAL_CODEBOOK_20260918" / "FINAL_SAFE_TERMINAL_CODEBOOK.tsv"
    internal_codebook = reports / "V3_SAFE_INTERNAL21_CODEBOOK_GATE_20260922" / "SAFE_INTERNAL21_INTERNAL_ONLY_CODEBOOK.tsv"
    mem = reports / "V3_SAFE_PAIR_PROVISIONAL_CANDIDATE_MEMORY_CROSSFIT_FALSIFICATION_20260922" / "results"
    rebuild = reports / "V3_FULL_MEMORY_ANCHORPE_A_REBUILD_FALSIFICATION_20260922" / "results"
    hostile = reports / "V3_FULL_MEMORY_ANCHORPE_HOSTILE_SPECIFICITY_AUDIT_20260922" / "results" / "FROZEN_A_ONLY_SELECTED_ANCHORS.tsv"
    r2 = reports / "V3_SAFE_INTERNAL21_ADDITIVE_PROPOSAL_MEMORY_TO_HELDOUT_ANCHORPE_FALSIFICATION_20260922_R2"
    obs = reports / "V3_SAFE_INTERNAL21_STRICT_A_JUNCTION_OBSERVABILITY_FALSIFICATION_20260923"
    report: dict[str, object] = {"truth_opened": False, "samples": {}}
    ok = True
    for sample in args.samples:
        started = time.time()
        spec = SAMPLES[sample]
        inputs = AnchorPEInputs(
            sample_id=sample, safe_codebook=codebook,
            te_fasta=shared / "TE consensus library" / "transposon_sequence_set.fa",
            te_family_map=args.route_b_v3 / "fbte_to_flybase_family_FROZEN.tsv",
            repeatmasker_bed=shared / "dm6_repeatmasker_TE_TEMP2.bed6",
            host_index_prefix=shared / "Reference" / "bowtie2_index" / "dm6",
            bowtie2=bowtie2, host_read_names=spec["profile"], threads=args.threads,
            edlib_python=args.edlib_python,
        )
        incumbent = (args.data / frozen.AKA_R9E) if sample == "AKA017" else (args.route_b_v2 / frozen.ISO_R9E)
        run = run_extension(inputs, internal_codebook, MsbwtSample(args.data / spec["msbwt"]), incumbent.read_bytes(),
                            frozen.LABELS[sample], args.work / sample,
                            progress=lambda m, s=sample: print(f"[{s} {time.time() - started:8.1f}s] {m}", flush=True))
        key = lambda c: (c.cluster_id, c.breakpoint0, c.anchor_host_side, c.family)  # noqa: E731
        fm = run.full_memory
        wkey = lambda r: (int(r["pair_id"]), r["family"], r["chrom"], int(r["breakpoint0"]), r["anchor_host_side"])  # noqa: E731
        instances = [{k: _text(e.get(k, "OBSERVED_INTERNAL_A" if k == "analysis_arm" else "")) for k in INSTANCE_FIELDS}
                     for e in run.internal21_events]
        loci = [{k: str(v) for k, v in x.ledger_row("OBSERVED_INTERNAL_A").items()} for x in run.internal21_loci]
        observed = lambda r: r["analysis_arm"] == "OBSERVED_INTERNAL_A"  # noqa: E731
        checks = {
            "anchor_pe_candidates": sorted(run.anchor_pe.candidates, key=key)
            == sorted(frozen.anchor_candidates(reports.parent / frozen.CANDIDATES[sample]), key=key),
            "full_proposals": [{k: _text(v) for k, v in r.items()} for r in fm.full_proposals]
            == _rows(mem / "STAGE_A_FULL_DATA_CANDIDATES.tsv", sample),
            "rotation_proposals": [{k: _text(v) for k, v in r.items()} for r in fm.rotation_proposals]
            == _rows(mem / "STAGE_A_CROSSFIT_CANDIDATES.tsv", sample),
            "full_memory_event_ledger": [{k: _text(r[k]) for k in EVENT_LEDGER_FIELDS} for r in fm.event_ledger]
            == _rows(rebuild / "FULL_MEMORY_ANCHORPE_EVENT_LEDGER.tsv", sample),
            "full_memory_selected_events": fm.events == frozen.full_memory_events(hostile, sample),
            "internal21_witnesses": [wkey(r) for r in run.internal21.witnesses]
            == [wkey(r) for r in _rows(r2 / f"{sample}_INTERNAL21_STRICT_A_WITNESSES.tsv", sample)],
            "internal21_instances": instances == _rows(obs / "INCREMENTAL_EVENT_INSTANCES.tsv", sample, observed),
            "internal21_loci": loci == _rows(obs / "INCREMENTAL_LOCUS_LEDGER.tsv", sample, observed),
        }
        layers = {}
        for name, layer, (rows, sha) in zip(LAYER_NAMES, run.layers, frozen.EXPECTED[sample][1:]):
            got = hashlib.sha256(layer.callset.data).hexdigest()
            layers[name] = {"rows": layer.callset.row_count, "sha256": got, "pass": layer.callset.row_count == rows and got == sha}
        passed = all(checks.values()) and all(v["pass"] for v in layers.values())
        ok &= passed
        report["samples"][sample] = {
            "host_read_names": spec["profile"], "anchor_pe_counts": run.anchor_pe.counts,
            "full_memory_counts": fm.counts, "internal21_counts": run.internal21.counts,
            "checks": checks, "layers": layers, "pass": passed, "seconds": round(time.time() - started, 1),
        }
        print(f"{sample}: checks {checks}", flush=True)
        print(f"{sample}: " + " ".join(f"{k}={v['rows']}:{v['pass']}" for k, v in layers.items()), flush=True)
    report["verdict"] = "V3_EXTENSION_ENGINE_EQUIVALENCE_PASS" if ok else "V3_EXTENSION_ENGINE_EQUIVALENCE_FAIL"
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "V3_EXTENSION_ENGINE_EQUIVALENCE.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(report["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
