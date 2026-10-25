"""AnchorPE candidate engine over a prepared sample (S20-S70 for layers 1-2)."""

from __future__ import annotations

import csv
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from ..constants import Q99
from ..errors import InputContractError
from ..incumbent.engine import StageStore
from ..shadow_layers import AnchorPECandidate
from .alignments import (
    HOST_READ_NAME_PROFILES,
    PRODUCTION_HOST_READ_NAMES,
    map_host_mates,
    map_reference_pairs,
    parse_host_sam,
    parse_reference_sam,
    write_host_fasta,
    write_pair_fastas,
)
from .anchors import channel_a, channel_b_universe, score_candidates
from .clusters import (
    annotate_corrected_clusters,
    cluster_mappings,
    load_family_map,
    parse_repeatmasker,
    relaxed_provisional_clusters,
)
from .codebook import classify_pairs, load_safe_codebook, one_mate_universe, project_codebook
from .decoder import (
    decode_fragments,
    decode_fragments_external,
    edlib_available,
    load_record_families,
    load_te_library,
)

__all__ = ["AnchorPEInputs", "AnchorPEResult", "HOST_READ_NAME_PROFILES", "derive_candidates", "run_anchor_pe"]


@dataclass(frozen=True, slots=True)
class AnchorPEInputs:
    sample_id: str
    safe_codebook: Path
    te_fasta: Path
    te_family_map: Path
    repeatmasker_bed: Path
    host_index_prefix: Path
    bowtie2: Path
    host_read_names: str = PRODUCTION_HOST_READ_NAMES
    threads: int = 1
    edlib_python: Path | None = None
    # 99th percentile of the sample's fragment size (bp); see ``sample_geometry``.
    fragment_q99: int = Q99

    def __post_init__(self) -> None:
        if self.host_read_names not in HOST_READ_NAME_PROFILES:
            raise InputContractError(f"host_read_names must be one of {HOST_READ_NAME_PROFILES}")


@dataclass(slots=True)
class AnchorPEResult:
    candidates: tuple[AnchorPECandidate, ...]
    candidate_rows: list[dict]
    corrected_clusters: list[dict]
    b_universe: list[dict]
    witnesses: list[dict]
    counts: dict[str, int] = field(default_factory=dict)
    # Shared evidence handed to the downstream (full-memory) branch.
    signature_rows: list[dict] = field(default_factory=list)
    eligible: dict[int, dict] = field(default_factory=dict)
    reads: dict[int, tuple[str, str]] = field(default_factory=dict)
    reference_sam: Path | None = None
    host_sams: dict[str, Path] = field(default_factory=dict)


def _fingerprint(inputs: AnchorPEInputs) -> dict[str, object]:
    def stat(path: Path) -> object:
        p = Path(path)
        return [str(p), p.stat().st_size, int(p.stat().st_mtime)] if p.is_file() else str(p)

    return {
        "schema": "evitea.anchor_pe.checkpoint.v1",
        "sample_id": inputs.sample_id,
        "safe_codebook": stat(inputs.safe_codebook),
        "te_fasta": stat(inputs.te_fasta),
        "te_family_map": stat(inputs.te_family_map),
        "repeatmasker_bed": stat(inputs.repeatmasker_bed),
        "host_index_prefix": str(inputs.host_index_prefix),
        "host_read_names": inputs.host_read_names,
        # recorded only when not the default, so earlier work directories stay resumable
        **({"fragment_q99": inputs.fragment_q99} if inputs.fragment_q99 != Q99 else {}),
    }


def _tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_anchor_pe(
    inputs: AnchorPEInputs,
    sample,
    work_directory: Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> AnchorPEResult:
    """Derive the AnchorPE candidates of one sample.

    ``sample`` is an ``incumbent.MsbwtSample`` over the prepared sample index;
    reads are recovered from it (no FASTQ is needed after sample preparation).
    """

    def say(message: str) -> None:
        if progress is not None:
            progress(message)

    work = Path(work_directory)
    work.mkdir(parents=True, exist_ok=True)
    store = StageStore(work / "checkpoints", _fingerprint(inputs))
    profile = inputs.host_read_names
    codebook = load_safe_codebook(inputs.safe_codebook)
    by_family = codebook.by_family()

    pair_hits = store.get("projection", lambda: project_codebook(sample, codebook, progress=say)[0])
    signature_rows, eligible = classify_pairs(pair_hits, codebook)
    universe = one_mate_universe(pair_hits, eligible)
    say(f"hit pairs={len(pair_hits):,}; one-mate family-resolved={len(eligible):,}; decoder universe={len(universe):,}")

    reads = store.get("reads", lambda: sample.pair_reads(universe))

    def decode() -> list[dict]:
        rows = [{"pair_id": pid, **{k: universe[pid][k] for k in ("te_mate", "host_mate", "unit_ids", "old_primary")},
                 "te_read": reads[pid][universe[pid]["te_mate"] - 1]} for pid in sorted(universe)]
        if not edlib_available() and inputs.edlib_python is not None:
            say(f"decoding {len(rows):,} fragments under {inputs.edlib_python}")
            return decode_fragments_external(
                Path(inputs.edlib_python), rows, codebook_path=Path(inputs.safe_codebook),
                te_fasta=Path(inputs.te_fasta), family_map=Path(inputs.te_family_map), work=work / "decoder")
        families = load_record_families(codebook, _tsv(inputs.te_family_map))
        return decode_fragments(rows, codebook, load_te_library(inputs.te_fasta), families)

    decoded = store.get("decoder", decode)
    strict_new = {int(r["pair_id"]): r for r in decoded
                  if int(r["strict_fragment_decodable"]) == 1 and not r["old_primary"]}
    say(f"strict whole-read decoded outside the old primary set: {len(strict_new):,}")

    def host_mapping(stem: str, host_mates: dict[int, int]) -> dict[int, dict]:
        fasta, sam = work / f"{stem}.fa", work / f"{stem}.sam"
        write_host_fasta(fasta, host_mates, reads, profile)
        map_host_mates(inputs.bowtie2, inputs.host_index_prefix, fasta, sam, inputs.threads)
        return parse_host_sam(sam)

    host = store.get("host", lambda: host_mapping(
        f"host_mates_{profile}", {pid: m["host_mate"] for pid, m in eligible.items()}))

    def reference():
        r1, r2, sam = work / "pairs_R1.fa", work / "pairs_R2.fa", work / "local_reference.sam"
        write_pair_fastas(r1, r2, eligible, reads)
        map_reference_pairs(inputs.bowtie2, inputs.host_index_prefix, r1, r2, sam, inputs.threads, inputs.fragment_q99)
        ledger, te_records = parse_reference_sam(sam, eligible, reads, by_family)
        # Only LOCAL_PROPER_TE_UNEXPLAINED records can witness Channel A.
        kept = {pid: recs for pid, recs in te_records.items()
                if ledger[pid]["reference_state"] == "LOCAL_PROPER_TE_UNEXPLAINED"}
        return ledger, kept

    ledger, te_records = store.get("reference", reference)
    decoded_host = store.get("decoded_host", lambda: host_mapping(
        f"decoded_host_mates_{profile}", {pid: int(r["host_mate"]) for pid, r in strict_new.items()})
        if strict_new else {})

    result = derive_candidates(
        ledger=ledger, te_records=te_records, host=host, decoded_host=decoded_host, strict_new=strict_new,
        by_family=by_family, family_map_rows=_tsv(inputs.te_family_map), repeatmasker_bed=inputs.repeatmasker_bed,
        q99=inputs.fragment_q99,
    )
    result.counts = {
        "codebook_hit_pairs": len(pair_hits), "one_mate_family_resolved": len(eligible),
        "decoder_universe": len(universe), **result.counts,
    }
    result.signature_rows, result.eligible, result.reads = signature_rows, eligible, reads
    result.reference_sam = work / "local_reference.sam"
    result.host_sams = {profile: work / f"host_mates_{profile}.sam"}
    say(f"AnchorPE counts: {result.counts}")
    return result


def derive_candidates(
    *,
    ledger: dict[int, dict],
    te_records: dict[int, list[dict]],
    host: dict[int, dict],
    decoded_host: dict[int, dict],
    strict_new: dict[int, dict],
    by_family,
    family_map_rows: list[dict[str, str]],
    repeatmasker_bed: Path,
    q99: int = Q99,
) -> AnchorPEResult:
    """Clusters, novelty/residency, Channel A, Channel B and candidates."""

    unique = [{"pair_id": pid, "family": m["family"], **{k: host[pid][k] for k in ("chrom", "start", "end", "strand")},
               "primary": m["primary"], "reference_state": m["reference_state"]}
              for pid, m in ledger.items() if host.get(pid, {}).get("host_class") == "UNIQUE_HIGH"]
    old_primary_unique = [r for r in unique if r["primary"]]
    primary_clusters = cluster_mappings(old_primary_unique, require_folds=False, prefix="PRIMARY")
    decoded_unique = [{"pair_id": pid, "family": r.get("family", ""),
                       **{k: decoded_host[pid][k] for k in ("chrom", "start", "end", "strand")},
                       "primary": False, "reference_state": "FRAGMENT_DECODED_EXPANSION"}
                      for pid, r in strict_new.items() if decoded_host.get(pid, {}).get("host_class") == "UNIQUE_HIGH"]
    relaxed = relaxed_provisional_clusters(primary_clusters, old_primary_unique, decoded_unique, strict_new, decoded_host)
    corrected = cluster_mappings([r for r in unique if r["reference_state"] != "EVIDENCE_COMPLETE_CLEAN_REFERENCE"],
                                 require_folds=True, prefix="CORRECTED", keep_all=True)
    fbmap, known = load_family_map(family_map_rows)
    corrected = annotate_corrected_clusters(corrected, primary_clusters + relaxed,
                                            parse_repeatmasker(repeatmasker_bed, fbmap, known))
    targets = [c for c in corrected if c["novel_nonresident"] and int(c.get("event_2p2fold", 0)) == 1]
    pair_to_cluster = {int(p): int(c["cluster_id"]) for c in corrected for p in str(c["pair_ids"]).split(";") if p}
    witnesses, anchors = channel_a(targets, pair_to_cluster, ledger, te_records, by_family, q99)
    b_rows = channel_b_universe(ledger, host)
    rows, candidates = score_candidates(anchors, b_rows, q99)
    counts = {
        "strict_decoded_outside_primary": len(strict_new),
        "primary_clusters": len(primary_clusters), "relaxed_provisional_clusters": len(relaxed),
        "corrected_clusters": len(corrected), "novel_nonresident_targets": len(targets),
        "replicated_A_anchors": len(anchors), "B_universe_unique_high": len(b_rows),
        "bilateral": sum(c.b_state == "BILATERAL" for c in candidates),
        "one_sided": sum(c.b_state in {"LEFT_ONLY", "RIGHT_ONLY"} for c in candidates),
    }
    return AnchorPEResult(tuple(candidates), rows, corrected, b_rows, witnesses, counts)

