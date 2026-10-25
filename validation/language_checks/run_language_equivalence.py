"""Validation-only: check the codebook builders against the recorded tables.

Each stage runs in isolation on its recorded inputs and is compared byte for
byte with its recorded outputs. Then the whole chain runs from the regenerated
anchor catalogue, and its two codebooks are compared with the recorded ones:

    final safe terminal codebook   685 rows     d631f16a...3e00
    internal-only safe codebook    12,268 rows  44b98960...62f7

Steps cache their results (``--cache``), so the harness can run in short
invocations (``--step``) or all at once (default).  No sample data, truth or
comparator file is read.

    python run_language_equivalence.py [--step g1a|rebuild|gapfill|final|internal|chain-...|report]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import pickle
import sys
import time
import zipfile
from pathlib import Path


def _env_path(name):
    """Path from environment variable ``name``; None when it is not set."""
    value = os.environ.get(name)
    return Path(value) if value else None

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from evitea.languages import g1a_tables  # noqa: E402
from evitea.languages.internal import build_internal21  # noqa: E402
from evitea.languages.terminal import fill_underpowered, finalize_terminal, rebuild_terminal21  # noqa: E402

FROZEN_INPUT_SHA = {
    "g1a_zip": "f08f1c1936566b0ca469a1c583ac9a96164b2fbd94e6d82f4235b88c49c6b695",
    "family_map": "6616b6305e6eac2f97336620376859151083a4b9ed30ba3bd64d7eead7c70714",
    "te_fasta": "8dff7f6c25b357d64184f232eff7e533058492ca7b118a06dac7d54bce13a922",
    "reference": "eeb896ecb227ec5b2c779faa46f064cf781adddaae509c7f50b9d6503c194397",
    "repeatmasker": "29ddba016d1afa6fa9b7c0fd45882ac5a3931408779580a946053f9bbab11e04",
}
FINAL_SHA = {
    "terminal_final": "d631f16a0252e4820fdfbe4f2822165595ebf726018618468602fea5dbfb3e00",
    "internal_only": "44b9896015eb1dd52b11d3a4c56c802f001af66d9d21d65cfd2f1a5ee6c662f7",
}
REBUILD = "reports/V3_REBUILT_SAFE_TERMINAL21_CODEBOOK_20260918"
GAPFILL = "reports/V3_UNDERPOWERED_TERMINAL_MIXED_K_GAP_FILL_20260918"
FINAL = "reports/V3_FINAL_SAFE_TERMINAL_CODEBOOK_20260918"
INTERNAL = "reports/V3_SAFE_INTERNAL21_CODEBOOK_GATE_20260922"
STEPS = ("g1a", "rebuild", "gapfill", "final", "internal", "chain-rebuild", "chain-gapfill", "chain-final",
         "chain-internal", "report")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


class Harness:
    def __init__(self, args):
        self.v3 = args.route_b_v3
        data = args.data / "TEMP2_ISO1"
        self.assets = {
            "te_fasta": data / "TE consensus library" / "transposon_sequence_set.fa",
            "family_map": self.v3 / "fbte_to_flybase_family_FROZEN.tsv",
            "reference_fasta": data / "Reference" / "dm6_TEMP2_compatible.fa",
            "repeatmasker_bed": data / "dm6_repeatmasker_TE_TEMP2.bed6",
        }
        self.cache = args.cache
        self.cache.mkdir(parents=True, exist_ok=True)
        self.report_path = self.cache / "LANGUAGE_EQUIVALENCE.json"
        self.report = json.loads(self.report_path.read_text()) if self.report_path.exists() else {"truth_opened": False}

    def record(self, step, checks):
        self.report[step] = checks
        self.report_path.write_text(json.dumps(self.report, indent=2, sort_keys=True) + "\n")
        print(step, json.dumps(checks), flush=True)

    def save(self, name, value):
        (self.cache / f"{name}.pkl").write_bytes(pickle.dumps(value))

    def load(self, name):
        return pickle.loads((self.cache / f"{name}.pkl").read_bytes())

    def frozen(self, rel, name) -> bytes:
        return (self.v3 / rel / name).read_bytes()

    def compare(self, produced: dict[str, bytes], rel: str) -> dict[str, object]:
        out = {}
        for name, data in produced.items():
            frozen = self.frozen(rel, name)
            out[name] = {"bytes_equal": data == frozen, "sha256": sha(data)[:16], "rows": data.count(b"\n") - 1}
        return out

    def zip_rows(self, member):
        with zipfile.ZipFile(self.v3 / "G1A_FROZEN_RESULTS_FOR_CODEBOOK.zip") as archive, archive.open(member) as raw:
            return list(csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8-sig", newline=""), delimiter="\t"))

    def frozen_g1a(self):
        path = self.cache / "frozen_g1a.pkl"
        if path.exists():
            return pickle.loads(path.read_bytes())
        value = (self.zip_rows("G1A_ANCHOR_CATALOGUE.tsv"), self.zip_rows("G1A_ANCHOR_PROVENANCE.tsv"))
        path.write_bytes(pickle.dumps(value))
        return value

    # -- steps -------------------------------------------------------------

    def step_g1a(self):
        inputs = {"g1a_zip": self.v3 / "G1A_FROZEN_RESULTS_FOR_CODEBOOK.zip", "family_map": self.assets["family_map"],
                  "te_fasta": self.assets["te_fasta"], "reference": self.assets["reference_fasta"],
                  "repeatmasker": self.assets["repeatmasker_bed"]}
        input_ok = {k: sha_file(p) == FROZEN_INPUT_SHA[k] for k, p in inputs.items()}
        catalogue, provenance = g1a_tables(self.assets["te_fasta"], self.assets["family_map"])
        frozen_cat, frozen_prov = self.frozen_g1a()
        cat_cols = list(catalogue[0])
        prov_cols = list(provenance[0])
        cat_eq = [tuple(r[c] for c in cat_cols) for r in catalogue] == [tuple(r[c] for c in cat_cols) for r in frozen_cat]
        prov_eq = [tuple(r[c] for c in prov_cols) for r in provenance] == [tuple(r[c] for c in prov_cols) for r in frozen_prov]
        self.save("g1a", (catalogue, provenance))
        self.record("g1a", {"frozen_inputs_sha256_ok": input_ok, "catalogue_rows": len(catalogue),
                            "catalogue_equal_on_consumed_columns": cat_eq, "provenance_rows": len(provenance),
                            "provenance_equal_on_consumed_columns": prov_eq})

    def step_rebuild(self, chain=False):
        catalogue, provenance = self.load("g1a") if chain else self.frozen_g1a()
        result = rebuild_terminal21(catalogue=catalogue, provenance=provenance, **self.assets,
                                    progress=lambda m: print(m, flush=True))
        self.save(("chain_" if chain else "") + "rebuild", result)
        self.record(("chain-" if chain else "") + "rebuild", {**self.compare({
            "SAFE_TERMINAL21_CANDIDATE_AUDIT.tsv": result.candidate_audit,
            "SAFE_TERMINAL21_SOURCE_AUDIT.tsv": result.source_audit,
            "REBUILT_SAFE_TERMINAL21_CODEBOOK.tsv": result.codebook,
            "REBUILT_SAFE_TERMINAL21_TE_RECORD_SUMMARY.tsv": result.te_record_summary,
            "REBUILT_SAFE_TERMINAL21_FAMILY_SUMMARY.tsv": result.family_summary,
        }, REBUILD), "counts": result.counts})

    def step_gapfill(self, chain=False):
        if chain:
            catalogue, provenance = self.load("g1a")
            prior = self.load("chain_rebuild")
            codebook, summary = prior.codebook, prior.te_record_summary
        else:
            catalogue, provenance = self.frozen_g1a()
            codebook = self.frozen(REBUILD, "REBUILT_SAFE_TERMINAL21_CODEBOOK.tsv")
            summary = self.frozen(REBUILD, "REBUILT_SAFE_TERMINAL21_TE_RECORD_SUMMARY.tsv")
        result = fill_underpowered(rebuilt_codebook=codebook, rebuilt_te_summary=summary, catalogue=catalogue,
                                   provenance=provenance, **self.assets, progress=lambda m: print(m, flush=True))
        self.save(("chain_" if chain else "") + "gapfill", result)
        self.record(("chain-" if chain else "") + "gapfill", {**self.compare({
            "UNDERPOWERED_MIXED_K_SOURCE_AUDIT.tsv": result.source_audit,
            "UNDERPOWERED_TARGET_SUMMARY.tsv": result.target_summary,
            "HYBRID_SAFE_TERMINAL_CODEBOOK.tsv": result.hybrid_codebook,
        }, GAPFILL), "counts": result.counts})

    def step_final(self, chain=False):
        if chain:
            rebuild, gapfill = self.load("chain_rebuild"), self.load("chain_gapfill")
            kwargs = dict(hybrid_codebook=gapfill.hybrid_codebook, rebuilt_source_audit=rebuild.source_audit,
                          gapfill_source_audit=gapfill.source_audit, rebuilt_te_summary=rebuild.te_record_summary)
        else:
            kwargs = dict(hybrid_codebook=self.frozen(GAPFILL, "HYBRID_SAFE_TERMINAL_CODEBOOK.tsv"),
                          rebuilt_source_audit=self.frozen(REBUILD, "SAFE_TERMINAL21_SOURCE_AUDIT.tsv"),
                          gapfill_source_audit=self.frozen(GAPFILL, "UNDERPOWERED_MIXED_K_SOURCE_AUDIT.tsv"),
                          rebuilt_te_summary=self.frozen(REBUILD, "REBUILT_SAFE_TERMINAL21_TE_RECORD_SUMMARY.tsv"))
        result = finalize_terminal(**kwargs)
        self.save(("chain_" if chain else "") + "final", result)
        self.record(("chain-" if chain else "") + "final", {**self.compare({
            "FINAL_SAFE_TERMINAL_CODEBOOK.tsv": result.codebook,
            "FINAL_SAFE_TERMINAL_RECORD_STATUS.tsv": result.record_status,
            "FINAL_SAFE_TERMINAL_REPAIR_LOG.tsv": result.repair_log,
        }, FINAL), "final_sha256_ok": sha(result.codebook) == FINAL_SHA["terminal_final"],
            "counts": result.counts})

    def step_internal(self, chain=False):
        if chain:
            catalogue, provenance = self.load("g1a")
            terminal = self.load("chain_final").codebook
        else:
            catalogue, provenance = self.frozen_g1a()
            terminal = self.frozen(FINAL, "FINAL_SAFE_TERMINAL_CODEBOOK.tsv")
        result = build_internal21(catalogue=catalogue, provenance=provenance, terminal_codebook=terminal,
                                  **self.assets, progress=lambda m: print(m, flush=True))
        self.save(("chain_" if chain else "") + "internal", result)
        self.record(("chain-" if chain else "") + "internal", {**self.compare({
            "SAFE_INTERNAL21_CANDIDATE_AUDIT.tsv": result.candidate_audit,
            "SAFE_INTERNAL21_SOURCE_AUDIT.tsv": result.source_audit,
            "SAFE_INTERNAL21_FULL_CODEBOOK.tsv": result.full_codebook,
            "SAFE_INTERNAL21_INTERNAL_ONLY_CODEBOOK.tsv": result.internal_only_codebook,
            "SAFE_INTERNAL21_TE_RECORD_SUMMARY.tsv": result.te_record_summary,
            "SAFE_INTERNAL21_FAMILY_SUMMARY.tsv": result.family_summary,
        }, INTERNAL), "internal_only_sha256_ok": sha(result.internal_only_codebook) == FINAL_SHA["internal_only"],
            "counts": result.counts})

    def step_report(self):
        def all_true(value):
            if isinstance(value, dict):
                return all(all_true(v) for k, v in value.items() if k not in ("sha256", "rows", "counts"))
            return value is True if isinstance(value, bool) else True

        needed = [s for s in STEPS if s != "report"]
        missing = [s for s in needed if s not in self.report]
        ok = not missing and all(all_true(self.report[s]) for s in needed)
        self.report["verdict"] = "LANGUAGE_BUILDER_EQUIVALENCE_PASS" if ok else "LANGUAGE_BUILDER_EQUIVALENCE_NOT_SHOWN"
        self.report["missing_steps"] = missing
        self.report_path.write_text(json.dumps(self.report, indent=2, sort_keys=True) + "\n")
        print(self.report["verdict"], missing)
        return 0 if ok else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=_env_path("EVITEA_DATA"), required=_env_path("EVITEA_DATA") is None)
    ap.add_argument("--route-b-v3", type=Path,
                    default=_env_path("EVITEA_ROUTE_B_V3"), required=_env_path("EVITEA_ROUTE_B_V3") is None)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--step", choices=STEPS, default=None, help="run one step (default: all)")
    args = ap.parse_args(argv)
    harness = Harness(args)
    steps = [args.step] if args.step else list(STEPS)
    for step in steps:
        started = time.time()
        if step == "report":
            return harness.step_report()
        chain = step.startswith("chain-")
        getattr(harness, "step_" + step.removeprefix("chain-"))(**({"chain": True} if chain else {}))
        print(f"[{step}] {time.time() - started:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
