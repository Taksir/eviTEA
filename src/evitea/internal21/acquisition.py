"""Internal 21-mer acquisition: internal-only safe 21-mers give internal exact-junction witnesses."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from ..anchor_pe.alignments import (
    PRODUCTION_HOST_READ_NAMES,
    map_host_mates,
    map_reference_pairs,
    write_host_fasta,
    write_pair_fastas,
)
from ..anchor_pe.codebook import revcomp
from ..anchor_pe.engine import AnchorPEInputs, _fingerprint
from ..constants import N_FOLDS
from ..errors import AssetCompatibilityError, BackendCompatibilityError
from ..full_memory.common import as_int, normfam
from ..full_memory.finite_state import finite_host_alternatives, finite_reference_alternatives
from ..full_memory.strict_a import allowed_pairs, derive_witnesses
from ..incumbent.engine import StageStore

__all__ = [
    "PER_ORIENTED_QUERY_ROW_CAP", "TOTAL_INTERVAL_ROW_CAP", "load_internal_codebook", "project_internal",
    "classify_internal_pairs", "internal_accounting", "Internal21Acquisition", "run_internal21_acquisition",
]

# Resource limits that stop the run (never a biological decision).
PER_ORIENTED_QUERY_ROW_CAP = 100_000
TOTAL_INTERVAL_ROW_CAP = 30_000_000


def load_internal_codebook(path: str | Path) -> dict[str, str]:
    """Map each canonical probe to its family, from the internal-only safe codebook."""

    probe_family: dict[str, str] = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if as_int(row.get("acquisition_safe")) != 1 or as_int(row.get("internal_only_query_unit")) != 1:
                raise AssetCompatibilityError("internal-only codebook contains a non-safe / non-internal row")
            probe = str(row["canonical_strand_probe"]).upper()
            if len(probe) != 21 or set(probe) - set("ACGT"):
                raise AssetCompatibilityError(f"malformed internal probe {probe!r}")
            family = normfam(row["family"])
            if probe_family.setdefault(probe, family) != family:
                raise AssetCompatibilityError(f"internal probe maps to several families: {probe}")
    return probe_family


def project_internal(sample, probe_family: Mapping[str, str],
                     progress: Callable[[str], None] | None = None) -> dict[int, dict[int, set[str]]]:
    """Map pair_id, then mate, to the observed canonical probes (both strands queried)."""

    intervals: dict[str, tuple[int, int]] = {}
    total = 0
    over = 0
    probes = sorted(probe_family)
    for probe in probes:
        for q in dict.fromkeys((probe, revcomp(probe))):
            lo, hi = sample.interval(q, cache=False)
            intervals[q] = (lo, hi)
            total += hi - lo
            over += int(hi - lo > PER_ORIENTED_QUERY_ROW_CAP)
    if over or total > TOTAL_INTERVAL_ROW_CAP:
        raise BackendCompatibilityError(
            f"internal 21-mer msBWT interval budget exceeded (total={total:,}, over-cap queries={over}); the run stops")
    hits: dict[int, dict[int, set[str]]] = defaultdict(lambda: defaultdict(set))
    for n, probe in enumerate(probes, 1):
        for q in dict.fromkeys((probe, revcomp(probe))):
            if intervals[q][1] == intervals[q][0]:
                continue
            files, pairs = sample.file_and_pair(sample.dollar_ids(q, cache=False))
            for file_id, pid in zip(files.tolist(), pairs.tolist()):
                hits[int(pid)][int(file_id) + 1].add(probe)
        if progress and (n % 2000 == 0 or n == len(probes)):
            progress(f"INTERNAL21 projection {n:,}/{len(probes):,}; pairs={len(hits):,}")
    return {pid: {1: set(m.get(1, ())), 2: set(m.get(2, ()))} for pid, m in hits.items()}


def classify_internal_pairs(hits: Mapping[int, Mapping[int, set[str]]], probe_family: Mapping[str, str]) -> tuple[dict[int, dict], Counter]:
    meta: dict[int, dict] = {}
    states: Counter = Counter()
    for pid in sorted(hits):
        m1, m2 = hits[pid][1], hits[pid][2]
        families = {probe_family[p] for p in m1 | m2}
        if len(families) != 1:
            states["FAMILY_CONFLICT"] += 1
        elif m1 and m2:
            states["BOTH_MATES_TE_POSITIVE"] += 1
        else:
            states["ONE_MATE_FAMILY_RESOLVED"] += 1
            te_mate = 1 if m1 else 2
            meta[pid] = {"pair_id": pid, "fold": pid % N_FOLDS, "family": next(iter(families)), "te_mate": te_mate,
                         "host_mate": 3 - te_mate, "observed_probes": set(m1 if te_mate == 1 else m2)}
    return meta, states


def pair_probe_sequences(meta_row: Mapping) -> tuple[str, ...]:
    return tuple(sorted({s for p in meta_row["observed_probes"] for s in (p, revcomp(p))}))


def internal_accounting(meta: Mapping[int, Mapping], host_summary: Mapping[int, dict], ref_summary: Mapping[int, dict]) -> list[dict]:
    rows = []
    for pid in sorted(meta):
        m, hs, rs = meta[pid], host_summary[pid], ref_summary[pid]
        state, reason = "ELIGIBLE", ""
        if hs["host_alt_count"] == 0:
            state, reason = "UNMAPPED", "HOST_UNMAPPED"
        elif hs["host_class"] == "MULTI_OR_LOW" and hs["host_alt_count"] < 2:
            state, reason = "HOST_UNRESOLVED", "LOW_CONFIDENCE_SINGLE_HOST_PLACEMENT"
        rows.append({
            "pair_id": pid, "fold": pid % N_FOLDS, "family": m["family"], "families": m["family"],
            "te_mate": m["te_mate"], "host_mate": m["host_mate"], "signature_state": "INTERNAL21_ONE_MATE_FAMILY_RESOLVED",
            "experiment_state": state, "exclusion_reason": reason, **rs, **hs,
            "observed_internal_probe_count": len(m["observed_probes"]),
            "observed_internal_probes": ";".join(sorted(m["observed_probes"])),
        })
    return rows


@dataclass(slots=True)
class Internal21Acquisition:
    witnesses: list[dict]
    accounting: list[dict]
    counts: dict[str, int] = field(default_factory=dict)


def run_internal21_acquisition(
    inputs: AnchorPEInputs,
    internal_codebook: Path,
    sample,
    work_directory: Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> Internal21Acquisition:
    def say(message: str) -> None:
        if progress is not None:
            progress(message)

    work = Path(work_directory)
    work.mkdir(parents=True, exist_ok=True)
    fingerprint = {**_fingerprint(inputs), "schema": "evitea.internal21.checkpoint.v1", "internal_codebook": str(internal_codebook)}
    store = StageStore(work / "checkpoints", fingerprint)
    probe_family = load_internal_codebook(internal_codebook)
    hits = store.get("projection", lambda: project_internal(sample, probe_family, say))
    meta, states = classify_internal_pairs(hits, probe_family)
    say(f"INTERNAL21 pairs={len(hits):,}; one-mate family-resolved={len(meta):,}")
    reads = store.get("reads", lambda: sample.pair_reads(meta))

    def state():
        stem = f"host_mates_{PRODUCTION_HOST_READ_NAMES}"
        host_fa, host_sam = work / f"{stem}.fa", work / f"{stem}.sam"
        write_host_fasta(host_fa, {p: m["host_mate"] for p, m in meta.items()}, reads, PRODUCTION_HOST_READ_NAMES)
        map_host_mates(inputs.bowtie2, inputs.host_index_prefix, host_fa, host_sam, inputs.threads)
        _, host_summary = finite_host_alternatives(host_sam, meta)
        r1, r2, ref_sam = work / "pairs_R1.fa", work / "pairs_R2.fa", work / "local_reference.sam"
        write_pair_fastas(r1, r2, meta, reads)
        map_reference_pairs(inputs.bowtie2, inputs.host_index_prefix, r1, r2, ref_sam, inputs.threads, inputs.fragment_q99)
        ref_rows, ref_summary = finite_reference_alternatives(
            ref_sam, meta, {}, reads, pair_probes={p: pair_probe_sequences(m) for p, m in meta.items()})
        accounting = internal_accounting(meta, host_summary, ref_summary)
        return accounting, [r for r in ref_rows if not r["evidence_complete"] and not r["reference_cap_censored"]]

    accounting, ref_rows = store.get("finite_state", state)
    witnesses = derive_witnesses(allowed_pairs(accounting), ref_rows)
    counts = {"internal_hit_pairs": len(hits), "one_mate_family_resolved": len(meta),
              "eligible_reference_unexplained": sum(r["experiment_state"] == "ELIGIBLE" and r["reference_state"] == "REFERENCE_UNEXPLAINED" for r in accounting),
              "internal_strict_a_witnesses": len(witnesses), **{f"signature_{k.lower()}": v for k, v in states.items()}}
    say(f"INTERNAL21 counts: {counts}")
    return Internal21Acquisition(witnesses, accounting, counts)
