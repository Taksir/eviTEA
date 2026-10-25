"""Validation-only: merge the recorded layer tables into the callset and check each layer hash.

Feeds the recorded candidate tables (the inputs of the four layer merges) into
``evitea.shadow_layers`` on top of the recorded library-layer callset bytes, and
checks the byte identity of every layer:

    AKA017  792, 801, 820, 835, 838
    ISO1    192, 192, 193, 194, 194

This isolates the merge of the layers into the callset from the upstream
evidence stages.  It reads no truth file, comparator file or file of any other
sample.  Production never imports it.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
import os
from pathlib import Path


def _env_path(name):
    """Path from environment variable ``name``; None when it is not set."""
    value = os.environ.get(name)
    return Path(value) if value else None

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from evitea.shadow_layers import (  # noqa: E402
    AnchorPECandidate,
    ExtensionLayerInputs,
    FullMemoryEvent,
    Internal21Instance,
    Internal21Locus,
    ShadowLabels,
    compose_frozen_shadow,
)

V3 = "reports"
AKA_R9E = Path("AKA017/route_b_r9e_validation_full55x/shadow_v2/reports/V2_R9E_REPEAT_AWARE_DIVERGENCE/R9E_COMBINED_CALLS.tsv")
ISO_R9E = Path("reports/DUAL_STRAIN_FORENSIC_20260917/ISO1_R9E_COMBINED_CALLS.frozen.tsv")
CANDIDATES = {
    "AKA017": f"{V3}/V3_OBSERVED_ANCHOR_BILATERAL_PE_COMPLETION_FALSIFICATION_20260920/ANCHOR_PE_CANDIDATES.tsv",
    "ISO1": f"{V3}/V3_ANCHOR_PE_ISO1_TRANSFER_20260920/ISO1_ANCHOR_PE_CANDIDATES.tsv",
}
FULL_MEMORY = f"{V3}/V3_FULL_MEMORY_ANCHORPE_HOSTILE_SPECIFICITY_AUDIT_20260922/results/FROZEN_A_ONLY_SELECTED_ANCHORS.tsv"
INT21_DIR = f"{V3}/V3_SAFE_INTERNAL21_STRICT_A_JUNCTION_OBSERVABILITY_FALSIFICATION_20260923"
# Recorded per-layer integration tables (decision and nearest existing call).
LEDGERS = {
    ("AKA017", "BILATERAL"): f"{V3}/V3_ANCHOR_PE_PRODUCTION_SHADOW_FREEZE_20260920/ANCHOR_PE_PRODUCTION_EVENTS.tsv",
    ("ISO1", "BILATERAL"): f"{V3}/V3_ANCHOR_PE_ISO1_TRANSFER_20260920/ISO1_ANCHOR_PE_PRODUCTION_EVENTS.tsv",
    ("AKA017", "ASYMMETRIC"): f"{V3}/V3_DUAL_DATASET_ASYMMETRIC_PRODUCTION_SHADOW_FREEZE_20260920/AKA017_ASYMMETRIC_ANCHOR_PE_EVENTS.tsv",
    ("ISO1", "ASYMMETRIC"): f"{V3}/V3_DUAL_DATASET_ASYMMETRIC_PRODUCTION_SHADOW_FREEZE_20260920/ISO1_ASYMMETRIC_ANCHOR_PE_EVENTS.tsv",
    ("AKA017", "FULL_MEMORY"): f"{V3}/V3_FULL_MEMORY_ANCHORPE_PRODUCTION_SHADOW_FREEZE_20260922/results/AKA017_FULL_MEMORY_ANCHOR_PE_EVENTS.tsv",
    ("ISO1", "FULL_MEMORY"): f"{V3}/V3_FULL_MEMORY_ANCHORPE_PRODUCTION_SHADOW_FREEZE_20260922/results/ISO1_FULL_MEMORY_ANCHOR_PE_EVENTS.tsv",
}
INT21_LEDGER = f"{V3}/V3_INTERNAL21_STRICT_A_INTEGRATION_SHADOW_FREEZE_20260923/INTEGRATION_LEDGER.tsv"
FORBIDDEN = ("truth", "teforest", "temp2_calls", "gim024", "gim-024", "quarantined")
LABELS = {
    "AKA017": ShadowLabels(bilateral_prefix="", asymmetric="AKA", full_memory="AKA", internal21="AKA"),
    "ISO1": ShadowLabels(bilateral_prefix="ISO1_", asymmetric="ISO1", full_memory="ISO", internal21="ISO"),
}
EXPECTED = {
    "AKA017": [
        (792, "1d1110f939996223379edf74847d626d4bc397a5cac32ad16d5d985dca051757"),
        (801, "ac352333be3912e54bfe8cbd6aaf9664231e34ba023d90c0a519d51c6a32a0f3"),
        (820, "8c12a5619d56dc603dc49c9a0b80b716056959fc39a12137d9c1b66bd883858a"),
        (835, "5efade8fdb64a41028419037a022a17a3380173cf87a1f7ac111c750f968b7c0"),
        (838, "6432c8e12db75f922c6fe44ff721b35be272c1eed34193364bbc4dc51160f701"),
    ],
    "ISO1": [
        (192, "b3e48681568d616142fe534995d3eeb9e5163907992e9d0893919aea42e6ebc5"),
        (192, "b3e48681568d616142fe534995d3eeb9e5163907992e9d0893919aea42e6ebc5"),
        (193, "31e15da3929500226351e4da86ca8afe1edc4faa2ef6b8fd5925c8aa349f0378"),
        (194, "9e538d744a6ea9c6cd2cd391740c742b75ef7ab1b4bd7a30a22f2736cc6092b0"),
        (194, "9e538d744a6ea9c6cd2cd391740c742b75ef7ab1b4bd7a30a22f2736cc6092b0"),
    ],
}
LAYERS = ("R9E", "bilateral", "asymmetric", "full-memory", "final INTERNAL21")


def _guard(path: Path) -> Path:
    lowered = str(path).lower()
    if any(token in lowered for token in FORBIDDEN):
        raise RuntimeError(f"blocked {path}: truth and comparator files must not be read")
    return path


def _tsv(path: Path) -> list[dict[str, str]]:
    with _guard(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _ids(raw: str) -> tuple[int, ...]:
    return tuple(int(float(x)) for x in str(raw or "").split(";") if x.strip())


def _check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def anchor_candidates(path: Path) -> tuple[AnchorPECandidate, ...]:
    out = []
    for row in _tsv(path):
        cand = AnchorPECandidate(
            cluster_id=int(row["cluster_id"]), family=row["family"], chrom=row["chrom"],
            breakpoint0=int(row["breakpoint0"]), anchor_host_side=row["anchor_host_side"],
            a_pair_ids=_ids(row["channel_A_pair_ids"]),
            left_b_pair_ids=_ids(row["left_B_pair_ids"]),
            right_b_pair_ids=_ids(row["right_B_pair_ids"]),
        )
        # The recorded table's summary columns must agree with the parsed ids.
        _check(int(row["channel_A_pairs"]) == len(cand.a_pair_ids), f"A count drift {row['cluster_id']}")
        _check(int(row["left_B_pairs"]) == len(cand.left_b_pair_ids), f"left count drift {row['cluster_id']}")
        _check(int(row["right_B_pairs"]) == len(cand.right_b_pair_ids), f"right count drift {row['cluster_id']}")
        _check(int(row["distinct_B_pairs"]) == len(cand.distinct_b_pair_ids), f"distinct drift {row['cluster_id']}")
        _check(int(row["bilateral_B_2p2fold"]) == int(cand.b_state == "BILATERAL"), f"B flag drift {row['cluster_id']}")
        out.append(cand)
    return tuple(out)


def full_memory_events(path: Path, sample: str) -> tuple[FullMemoryEvent, ...]:
    out = []
    for row in _tsv(path):
        if row["sample"] != sample:
            continue
        event = FullMemoryEvent(
            candidate_id=row["candidate_id"], family=row["family"], chrom=row["chrom"],
            breakpoint0=int(row["breakpoint0"]), anchor_host_side=row["anchor_host_side"],
            a_pair_ids=_ids(row["A_pair_ids"]),
            left_b_pair_ids=_ids(row["target_left_pair_ids"]),
            right_b_pair_ids=_ids(row["target_right_pair_ids"]),
            b_state=row["target_B_state"],
        )
        _check(int(row["target_pair_mentions"]) == event.b_pair_mentions,
               f"target_pair_mentions != |L|+|R| for {row['candidate_id']}")
        out.append(event)
    return tuple(out)


def internal21_loci(directory: Path, sample: str) -> tuple[Internal21Locus, ...]:
    instances: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in _tsv(directory / "INCREMENTAL_EVENT_INSTANCES.tsv"):
        if row["sample"] == sample and row["analysis_arm"] == "OBSERVED_INTERNAL_A":
            instances[(row["family"], row["chrom"])].append(row)
    out = []
    for row in _tsv(directory / "INCREMENTAL_LOCUS_LEDGER.tsv"):
        if (row["sample"], row["analysis_arm"], row["classification"]) != (
            sample, "OBSERVED_INTERNAL_A", "NOVEL_UNCONFOUNDED"
        ):
            continue
        members = tuple(
            Internal21Instance(
                rotation=r["rotation"], breakpoint0=int(r["breakpoint0"]),
                anchor_host_side=r["anchor_host_side"], a_pair_ids=_ids(r["A_pair_ids"]),
                left_b_pair_ids=_ids(r["left_B_pair_ids"]), right_b_pair_ids=_ids(r["right_B_pair_ids"]),
            )
            for (family, chrom), rows in instances.items()
            if chrom == row["chrom"] and family == row["family"]
            for r in rows
        )
        out.append(Internal21Locus(row["family"], row["chrom"], int(row["locus_lo"]), int(row["locus_hi"]), members))
    return tuple(out)


def _status(raw: str) -> str:
    if raw.startswith("ADD"):
        return "ADDED"
    if raw.startswith("UNRESOLVED"):
        return "UNRESOLVED"
    return "ABSORBED"


def _dist(raw) -> str:
    return "" if raw in (None, "") else str(int(raw))


def ledger_mismatches(v3: Path, sample: str, results) -> list[str]:
    """Compare every integration decision and nearest-call attribution."""

    problems: list[str] = []
    by_layer = {r.layer: r for r in results}
    for layer in ("BILATERAL", "ASYMMETRIC", "FULL_MEMORY"):
        frozen = {}
        for row in _tsv(v3 / LEDGERS[(sample, layer)]):
            eid = row.get("anchor_pe_event_id") or row["event_id"]
            dist = row.get("nearest_incumbent_same_family_distance", row.get("nearest_existing_same_family_distance"))
            near = row.get("nearest_incumbent_event_id", row.get("nearest_existing_event_id"))
            frozen[eid] = (_status(row["integration_status"]), _dist(dist), near)
        mine = {
            r.event_id: (_status(r.status), _dist(r.nearest_same_family_distance), r.nearest_same_family_event_id)
            for r in by_layer[layer].ledger
        }
        if mine != frozen:
            diff = sorted(set(mine.items()) ^ set(frozen.items()))[:5]
            problems.append(f"{sample} {layer} ledger differs: {diff}")
    frozen = {}
    for row in _tsv(v3 / INT21_LEDGER):
        if row["sample"] != sample:
            continue
        key = f"{row['family']}:{row['chrom']}:{row['locus_lo']}-{row['locus_hi']}"
        frozen[key] = (
            _status(row["integration_status"]),
            _dist(row["nearest_same_family_base_distance"]), row["nearest_same_family_base_event_id"],
            _dist(row["nearest_different_family_base_distance"]), row["nearest_different_family_base_event_id"],
            row["chosen_breakpoint0"],
        )
    mine = {
        r.event_key: (
            _status(r.status), _dist(r.nearest_same_family_distance), r.nearest_same_family_event_id,
            _dist(r.nearest_different_family_distance), r.nearest_different_family_event_id,
            "" if r.breakpoint0 is None else str(r.breakpoint0),
        )
        for r in by_layer["INTERNAL21"].ledger
    }
    if mine != frozen:
        diff = sorted(set(mine.items()) ^ set(frozen.items()))[:5]
        problems.append(f"{sample} INTERNAL21 ledger differs: {diff}")
    return problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=_env_path("EVITEA_DATA"), required=_env_path("EVITEA_DATA") is None)
    ap.add_argument("--route-b-v2", type=Path, required=True)
    ap.add_argument("--route-b-v3", type=Path, required=True)
    ap.add_argument("--incumbent-aka017", type=Path, default=None,
                    help="library-layer callset to build on (default: the recorded callset; pass the engine output to chain it)")
    ap.add_argument("--output", type=Path, default=None, help="optional JSON report path")
    args = ap.parse_args(argv)

    incumbents = {
        "AKA017": args.incumbent_aka017 or args.data / AKA_R9E,
        "ISO1": args.route_b_v2 / ISO_R9E,
    }
    report: dict[str, object] = {"datasets": {}, "truth_opened": False}
    ok = True
    for sample in ("AKA017", "ISO1"):
        inputs = ExtensionLayerInputs(
            bilateral_candidates=anchor_candidates(args.route_b_v3 / CANDIDATES[sample]),
            full_memory_events=full_memory_events(args.route_b_v3 / FULL_MEMORY, sample),
            internal21_loci=internal21_loci(args.route_b_v3 / INT21_DIR, sample),
        )
        incumbent = _guard(incumbents[sample]).read_bytes()
        results = compose_frozen_shadow(incumbent, inputs, LABELS[sample], sample_id=sample)
        callsets = [incumbent] + [r.callset.data for r in results]
        layers = []
        for name, data, (rows, sha) in zip(LAYERS, callsets, EXPECTED[sample]):
            got_rows = data.count(b"\n") - 1
            got_sha = hashlib.sha256(data).hexdigest()
            passed = got_rows == rows and got_sha == sha
            ok &= passed
            layers.append({"layer": name, "rows": got_rows, "sha256": got_sha, "expected_rows": rows, "pass": passed})
        problems = ledger_mismatches(args.route_b_v3, sample, results)
        ok &= not problems
        report["datasets"][sample] = {
            "layers": layers,
            "ledger_mismatches": problems,
            "unresolved": [
                {"layer": u.layer, "key": u.event_key, "reason": u.reason}
                for r in results for u in r.unresolved
            ],
            "appended": {r.layer: list(r.appended_event_ids) for r in results},
        }
    report["verdict"] = "SHADOW_LAYER_INTEGRATION_EQUIVALENCE_PASS" if ok else "SHADOW_LAYER_INTEGRATION_EQUIVALENCE_FAIL"
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    for sample, entry in report["datasets"].items():
        for layer in entry["layers"]:
            print(f"{sample:7s} {layer['layer']:17s} rows={layer['rows']:4d} pass={layer['pass']} {layer['sha256'][:16]}")
        for problem in entry["ledger_mismatches"]:
            print(f"{sample:7s} LEDGER MISMATCH {problem}")
        print(f"{sample:7s} integration ledgers: {'exact' if not entry['ledger_mismatches'] else 'MISMATCH'}")
        for u in entry["unresolved"]:
            print(f"{sample:7s} unresolved {u}")
    print(report["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
