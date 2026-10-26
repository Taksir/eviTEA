"""Validation-only readers for the recorded stage tables (TSV, NPY, zip).

They let each library-layer stage be compared with recorded outputs.
Production code never imports this module.
"""

from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path

import numpy as np

from evitea.incumbent.candidates import HostMateRecord, PairPanelEvidence


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def zip_tsv(archive: Path, member: str) -> list[dict[str, str]]:
    with zipfile.ZipFile(archive) as bundle, bundle.open(member) as raw:
        return list(csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8", newline=""), delimiter="\t"))


def zip_npy(archive: Path, member: str) -> np.ndarray:
    with zipfile.ZipFile(archive) as bundle:
        return np.load(io.BytesIO(bundle.read(member)), allow_pickle=False)


def _truthy(value: str | None) -> bool:
    return str(value if value is not None else "0").strip() in {"1", "True", "true"}


def _host(rows: list[dict[str, str]]) -> dict[int, HostMateRecord]:
    return {
        int(row["pair_id"]): HostMateRecord(
            int(row["pair_id"]),
            row.get("status", ""),
            row.get("chrom", "") or "",
            row.get("midpoint", "") or "",
            row.get("strand", "") or "",
        )
        for row in rows
    }


def _evidence(rows: list[dict[str, str]]) -> dict[int, PairPanelEvidence]:
    return {
        int(row["pair_id"]): PairPanelEvidence(
            int(row["pair_id"]), row.get("families", "") or "", _truthy(row.get("mate1_hit")), _truthy(row.get("mate2_hit"))
        )
        for row in rows
    }


def load_candidate_universe_inputs(v1_reports: Path) -> dict[str, object]:
    """Load the recorded candidate inputs from a report root."""

    g1b = v1_reports / "G1B_sparse_panel_routing" / "TE_LOCATOR_ROUTE_B_G1B_RESULTS_FOR_CHATGPT.zip"
    g1f = v1_reports / "G1F_host_mate_locus_gate" / "TE_LOCATOR_ROUTE_B_G1F_RESULTS_FOR_CHATGPT.zip"
    g1g = v1_reports / "G1G_reference_resident_compatibility" / "G1G_HOST_CLUSTER_PAIR_FAMILIES.tsv"
    if not g1g.is_file():
        g1g_zip = v1_reports / "G1G_reference_resident_compatibility" / "TE_LOCATOR_ROUTE_B_G1G_RESULTS_FOR_CHATGPT.zip"
        families_rows = zip_tsv(g1g_zip, "G1G_HOST_CLUSTER_PAIR_FAMILIES.tsv")
    else:
        families_rows = read_tsv(g1g)
    g1hr = v1_reports / "G1HR_candidate_recall_repair" / "TE_LOCATOR_ROUTE_B_G1HR_RESULTS_FOR_CHATGPT.zip"
    return {
        "final_pair_ids": zip_npy(g1hr, "G1HR_REPAIRED_CANDIDATE_PAIR_IDS.npy").tolist(),
        "working_pair_ids": zip_npy(g1f, "G1F_WORKING_CANDIDATE_PAIR_IDS.npy").tolist(),
        "primary_host": _host(zip_tsv(g1f, "G1F_HOST_MATE_MAPPING.tsv")),
        "rescue_host": _host(zip_tsv(g1hr, "G1HR_NEW_HOST_MAPPING.tsv")),
        "resident_families": {int(row["pair_id"]): row.get("families", "") for row in families_rows},
        "rescue_evidence": _evidence(zip_tsv(g1hr, "G1HR_NEW_PAIR_EVIDENCE.tsv")),
        "panel_evidence": _evidence(zip_tsv(g1b, "G1B_CANDIDATE_PAIRS.tsv")),
    }


def legacy_event_rows(path: Path) -> list[dict[str, str]]:
    return read_tsv(path)
