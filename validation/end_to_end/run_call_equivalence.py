"""Validation-only: run the public ``evitea call`` pipeline end to end, per sample.

Prepares the reference and TE assets once, with hash checks, as a user would.
``prepare-reference`` adopts the existing dm6 Bowtie2 index. ``prepare-te``
either builds both codebooks (``--build-languages``) or binds the recorded
codebook files (``--languages supplied``). The script checks the codebooks
against their recorded SHA-256 values. Then it runs

    python -m evitea call --msbwt-index <sample msBWT> --sample-id S --profile S ...

and compares the SHA-256 of each layer (recorded in ``run_summary.json``) and of
``calls.native.tsv`` with the recorded callsets. The call counts per layer, from
the library layer to the internal 21-mer layer, and the final SHA-256 are:

    AKA017  792, 801, 820, 835, 838  (6432c8e1...f701)
    ISO1    192, 192, 193, 194, 194  (9e538d74...92b0)

The first ISO1 layer has the recorded SHA-256 b3e48681...
Resumable: an existing work directory is resumed with ``--resume``.  No truth
file, comparator file or file of any other sample is read.  Production code
never imports this module.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path


def _env_path(name):
    """Path from environment variable ``name``; None when it is not set."""
    value = os.environ.get(name)
    return Path(value) if value else None

ROOT = Path(__file__).resolve().parents[2]
DATA = _env_path("EVITEA_DATA")
ROUTE_B_V3 = _env_path("EVITEA_ROUTE_B_V3")

EXPECTED = {
    "AKA017": {
        "incumbent": (792, "1d1110f939996223379edf74847d626d4bc397a5cac32ad16d5d985dca051757"),
        "bilateral": (801, "ac352333be3912e54bfe8cbd6aaf9664231e34ba023d90c0a519d51c6a32a0f3"),
        "asymmetric": (820, "8c12a5619d56dc603dc49c9a0b80b716056959fc39a12137d9c1b66bd883858a"),
        "full_memory": (835, "5efade8fdb64a41028419037a022a17a3380173cf87a1f7ac111c750f968b7c0"),
        "internal21": (838, "6432c8e12db75f922c6fe44ff721b35be272c1eed34193364bbc4dc51160f701"),
    },
    "ISO1": {
        "incumbent": (192, "b3e48681568d616142fe534995d3eeb9e5163907992e9d0893919aea42e6ebc5"),
        "bilateral": (192, "b3e48681568d616142fe534995d3eeb9e5163907992e9d0893919aea42e6ebc5"),
        "asymmetric": (193, "31e15da3929500226351e4da86ca8afe1edc4faa2ef6b8fd5925c8aa349f0378"),
        "full_memory": (194, "9e538d744a6ea9c6cd2cd391740c742b75ef7ab1b4bd7a30a22f2736cc6092b0"),
        "internal21": (194, "9e538d744a6ea9c6cd2cd391740c742b75ef7ab1b4bd7a30a22f2736cc6092b0"),
    },
}
MSBWT = {"AKA017": "AKA017/msbwt_SRR9951105", "ISO1": "TEMP2_ISO1/msbwt_SRR11906527"}
FROZEN_LANGUAGES = {
    "terminal": ("reports/V3_FINAL_SAFE_TERMINAL_CODEBOOK_20260918/FINAL_SAFE_TERMINAL_CODEBOOK.tsv",
                 "d631f16a0252e4820fdfbe4f2822165595ebf726018618468602fea5dbfb3e00"),
    "internal": ("reports/V3_SAFE_INTERNAL21_CODEBOOK_GATE_20260922/SAFE_INTERNAL21_INTERNAL_ONLY_CODEBOOK.tsv",
                 "44b9896015eb1dd52b11d3a4c56c802f001af66d9d21d65cfd2f1a5ee6c662f7"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def evitea(python: str, *args: str) -> None:
    command = [python, "-m", "evitea", *args]
    print(">", " ".join(f'"{a}"' if " " in a else a for a in command), flush=True)
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(ROOT / "src"), env.get("PYTHONPATH"))))
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(command, cwd=ROOT, env=env)
    if completed.returncode != 0:
        raise SystemExit(f"evitea exited with {completed.returncode}")


def prepare_assets(args, python: str) -> tuple[Path, Path]:
    shared = args.data / "TEMP2_ISO1"
    reference = args.work / "assets" / "reference"
    te = args.work / "assets" / "te"
    if not (reference / "reference_manifest.json").is_file():
        evitea(python, "prepare-reference",
               "--reference", str(shared / "Reference" / "dm6_TEMP2_compatible.fa"),
               "--repeat-annotation", str(shared / "dm6_repeatmasker_TE_TEMP2.bed6"),
               "--bowtie2-index", str(shared / "Reference" / "bowtie2_index" / "dm6"),
               "--output", str(reference))
    if not (te / "te_manifest.json").is_file():
        common = ["prepare-te",
                  "--te-library", str(shared / "TE consensus library" / "transposon_sequence_set.fa"),
                  "--family-map", str(shared / "fbte_to_flybase_family.tsv"),
                  "--reference-assets", str(reference), "--output", str(te)]
        if args.languages == "built":
            evitea(python, *common, "--build-languages")
        else:
            languages = {}
            for key, (relative, expected) in FROZEN_LANGUAGES.items():
                path = args.route_b_v3 / relative
                if sha256(path) != expected:
                    raise SystemExit(f"{key} codebook hash does not match the recorded hash: {path}")
                languages[key] = path
            evitea(python, *common, "--terminal-codebook", str(languages["terminal"]),
                   "--internal-codebook", str(languages["internal"]))
    evitea(python, "validate", "reference", "--reference-assets", str(reference))
    evitea(python, "validate", "te", "--te-assets", str(te))
    manifest = json.loads((te / "te_manifest.json").read_text(encoding="utf-8"))
    bound = {key: entry["sha256"] for key, entry in manifest["languages"].items()}
    expected = {"terminal_final": FROZEN_LANGUAGES["terminal"][1], "internal_only": FROZEN_LANGUAGES["internal"][1]}
    if bound != expected:
        raise SystemExit(f"prepared codebooks do not match the recorded hashes: {bound}")
    print("languages:", {k: v["source"] for k, v in manifest["languages"].items()}, "= frozen hashes", flush=True)
    return reference, te


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=DATA, required=DATA is None)
    ap.add_argument("--route-b-v3", type=Path, default=ROUTE_B_V3, required=ROUTE_B_V3 is None)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--bowtie2", type=Path, default=None)
    ap.add_argument("--edlib-python", type=Path, default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--samples", nargs="+", default=["ISO1", "AKA017"], choices=sorted(EXPECTED))
    ap.add_argument("--languages", choices=("built", "supplied"), default="built",
                    help="build both codebooks with EviTEA (default) or bind the recorded files")
    args = ap.parse_args(argv)
    args.work.mkdir(parents=True, exist_ok=True)
    bowtie2 = args.bowtie2 or args.data / "TEMP2_ISO1" / "tools" / "bowtie2" / "bowtie2.bat"
    reference, te = prepare_assets(args, args.python)
    report: dict[str, object] = {"truth_opened": False, "languages": args.languages, "samples": {}}
    ok = True
    for sample in args.samples:
        started = time.time()
        out, work = args.work / sample / "out", args.work / sample / "work"
        command = [
            "call", "--msbwt-index", str(args.data / MSBWT[sample]), "--sample-id", sample,
            "--profile", sample, "--reference-assets", str(reference), "--te-assets", str(te),
            "--output", str(out), "--work-dir", str(work), "--bowtie2", str(bowtie2),
            "--threads", str(args.threads),
        ]
        if args.workers:
            command += ["--workers", str(args.workers)]
        if args.edlib_python:
            command += ["--edlib-python", str(args.edlib_python)]
        if (out / "run_manifest.json").is_file():
            print(f"{sample}: completed run found in {out}; re-checking it", flush=True)
        else:
            if work.exists():
                command.append("--resume")
            evitea(args.python, *command)
        evitea(args.python, "validate", "run", "--output", str(out))
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        native_sha = sha256(out / "calls.native.tsv")
        layers = {}
        for name, (rows, sha) in EXPECTED[sample].items():
            got_rows = summary["layer_rows"].get(name)
            got_sha = summary["layer_sha256"].get(name)
            layers[name] = {"rows": got_rows, "sha256": got_sha, "pass": got_rows == rows and got_sha == sha}
        final_rows, final_sha = EXPECTED[sample]["internal21"]
        checks = {
            "native_callset_sha256": native_sha == final_sha == summary["native_callset_sha256"],
            "final_call_count": summary["final_call_count"] == final_rows,
            "profile": summary["profile"] == sample,
        }
        passed = all(checks.values()) and all(v["pass"] for v in layers.values())
        ok &= passed
        report["samples"][sample] = {
            "checks": checks, "layers": layers, "pass": passed,
            "calls_by_origin": summary["calls_by_origin"], "unresolved_count": summary["unresolved_count"],
            "seconds_this_invocation": round(time.time() - started, 1),
        }
        print(f"{sample}: " + " ".join(f"{k}={v['rows']}:{v['pass']}" for k, v in layers.items()), flush=True)
    report["verdict"] = "EVITEA_CALL_END_TO_END_EQUIVALENCE_PASS" if ok else "EVITEA_CALL_END_TO_END_EQUIVALENCE_FAIL"
    (args.work / "EVITEA_CALL_END_TO_END_EQUIVALENCE.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(report["verdict"], flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
