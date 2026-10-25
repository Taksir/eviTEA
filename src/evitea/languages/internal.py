"""Safe internal-21 codebooks (full and internal-only).

Selects the anchor catalogue's k=21 internal candidates
(``has_internal_provenance``, not low complexity) under the terminal
rebuild's exhaustive TE-family and RepeatMasker-aware reference safety.
Units in the final terminal codebook are excluded from the
internal-only codebook.  No sample data are read and no selection threshold
is introduced.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from .common import (
    build_interval_index,
    extract_te_id_rebuild,
    load_repeatmasker_rebuild,
    oriented_code_map,
    overlap_query,
    parse_family_map_rebuild,
    parse_fasta,
    read_tsv_bytes,
    reference_hit_class,
    scan_exact,
    tsv_bytes,
)

__all__ = ["InternalResult", "build_internal21"]

K = 21
# Default block size of the ``read_length_block_125`` column (the column keeps its name; ``read_length`` sets the size).
READ_LENGTH = 125


@dataclass(slots=True)
class InternalResult:
    candidate_audit: bytes
    source_audit: bytes
    full_codebook: bytes
    internal_only_codebook: bytes
    te_record_summary: bytes
    family_summary: bytes
    counts: dict[str, int] = field(default_factory=dict)


def build_internal21(
    *,
    catalogue: Sequence[Mapping[str, str]],
    provenance: Sequence[Mapping[str, str]],
    te_fasta: str | Path,
    family_map: str | Path,
    reference_fasta: str | Path,
    repeatmasker_bed: str | Path,
    terminal_codebook: bytes,
    progress=None,
    read_length: int = READ_LENGTH,
) -> InternalResult:
    say = progress or (lambda message: None)
    family_by_te: dict[str, str] = {}
    for r in provenance:
        old = family_by_te.setdefault(r["te_id"], r["family"])
        if old != r["family"]:
            raise RuntimeError(f"Inconsistent anchor family for {r['te_id']}")
    cat_by_qid, candidates = {}, set()
    for r in catalogue:
        if int(r["k"]) != K or int(r["has_internal_provenance"]) != 1 or int(r["low_complexity"]) != 0:
            continue
        cat_by_qid[r["query_id"]] = r
        candidates.add(r["canonical_strand_key"])
    source_rows, seen = [], set()
    for p in provenance:
        if int(p["k"]) != K or p["region"] not in {"internal", "short_internal_fallback"}:
            continue
        c = cat_by_qid.get(p["query_id"])
        if c is None:
            continue
        row = {"canonical_strand_probe": c["canonical_strand_key"], "te_id": p["te_id"], "family": p["family"],
               "source_region": p["region"], "source_start": int(p["start"]), "source_end": int(p["start"]) + K}
        key = tuple(row.values())
        if key not in seen:
            seen.add(key)
            source_rows.append(row)
    say(f"internal k21 candidates={len(candidates):,}; source rows={len(source_rows):,}")
    source_by_probe = defaultdict(list)
    for r in source_rows:
        source_by_probe[r["canonical_strand_probe"]].append(r)

    qmap = oriented_code_map(candidates)
    te_hits, te_len = defaultdict(list), {}
    for header, seq in parse_fasta(te_fasta):
        te = extract_te_id_rebuild(header)
        te_len[te] = len(seq)
        family = family_by_te.get(te, te)
        for st, en, probe in scan_exact(seq, qmap, K):
            te_hits[probe].append((te, family, st, en))
    missing = sorted(set(family_by_te) - set(te_len))
    if missing:
        raise RuntimeError(f"TE library missing anchor catalogue records: {missing[:20]}")

    probe_state, family_safe, failures = {}, set(), Counter()
    for probe in sorted(candidates):
        hits = te_hits.get(probe, [])
        hit_families = sorted({fam for _, fam, _, _ in hits})
        src = source_by_probe.get(probe, [])
        src_families = sorted({r["family"] for r in src})
        single = len(hit_families) == 1
        consistent = len(src_families) == 1 and single and src_families[0] == hit_families[0]
        if not hits:
            failures["NO_TE_LIBRARY_HIT"] += 1
        elif not single:
            failures["MULTI_FAMILY_TE_LIBRARY"] += 1
        elif not consistent:
            failures["SOURCE_FAMILY_MISMATCH"] += 1
        else:
            family_safe.add(probe)
        probe_state[probe] = {
            "canonical_strand_probe": probe, "source_families": src_families,
            "source_te_ids": sorted({r["te_id"] for r in src}), "exhaustive_hit_families": hit_families,
            "exhaustive_hit_te_ids": sorted({te for te, _, _, _ in hits}),
            "exhaustive_single_family": int(single), "source_family_consistent": int(consistent),
            "exhaustive_te_occurrence_count": len(hits),
        }

    rm_index = build_interval_index(load_repeatmasker_rebuild(repeatmasker_bed, parse_family_map_rebuild(family_map)))
    qmap_ref = oriented_code_map(family_safe)
    target = {}
    for probe in family_safe:
        src = source_by_probe[probe]
        families = sorted({r["family"] for r in src})
        if len(families) != 1:
            raise RuntimeError("Internal source-family invariant failed")
        target[probe] = (families[0], sorted({r["te_id"] for r in src})[0])
    ref_count, context = Counter(), defaultdict(Counter)
    for chrom, seq in parse_fasta(reference_fasta):
        for st, en, probe in scan_exact(seq, qmap_ref, K):
            family, te = target[probe]
            ref_count[probe] += 1
            context[probe][reference_hit_class(family, te, st, en, overlap_query(rm_index, chrom, st, en))] += 1
        say(f"reference scanned: {chrom}")

    terminal_rows = read_tsv_bytes(terminal_codebook)
    terminal_probes = {r["canonical_strand_probe"] for r in terminal_rows}
    terminal_families = {r["family"] for r in terminal_rows}
    terminal_tes = {r["te_id"] for r in terminal_rows}

    candidate_audit, safe_units = [], set()
    for probe in sorted(candidates):
        s = probe_state[probe]
        total, c = ref_count[probe], context[probe]
        if probe not in family_safe:
            ref_safe, ref_class = False, "NOT_EVALUATED_SOURCE_UNSAFE"
        elif total == 0:
            ref_safe, ref_class = True, "REFERENCE_ABSENT"
        else:
            ref_safe = total - c["TARGET_CONTEXT_ANNOTATED_TE"] == 0
            ref_class = "TARGET_CONTEXT_ANNOTATED_TE_ONLY" if ref_safe else "UNSAFE_OR_MIXED_REFERENCE_CONTEXT"
        safe = bool(s["source_family_consistent"] and ref_safe)
        if safe:
            safe_units.add(probe)
        row = dict(s)
        row.update({
            "reference_exact_occurrence_count": total,
            "reference_target_context_hits": c["TARGET_CONTEXT_ANNOTATED_TE"],
            "reference_unsafe_context_hits": total - c["TARGET_CONTEXT_ANNOTATED_TE"],
            "reference_class": ref_class, "acquisition_safe": int(safe),
            "present_in_final_terminal_codebook": int(probe in terminal_probes),
            "internal_only_query_unit": int(safe and probe not in terminal_probes),
        })
        candidate_audit.append(row)
    by_probe = {r["canonical_strand_probe"]: r for r in candidate_audit}

    source_audit, safe_by_te = [], defaultdict(list)
    for s in source_rows:
        probe = s["canonical_strand_probe"]
        ps = by_probe[probe]
        positions = sorted({(h[2], h[3]) for h in te_hits.get(probe, []) if h[0] == s["te_id"]})
        match = (s["source_start"], s["source_end"]) in positions
        length = te_len[s["te_id"]]
        row = dict(s)
        row.update({
            "te_record_length": length,
            "nearest_terminal_distance_bp": min(s["source_start"], length - s["source_end"]),
            "read_length_block_125": s["source_start"] // read_length,
            "source_exact_position_match": int(match),
            "source_te_distinct_occurrence_count": len(positions),
            "reference_exact_occurrence_count": ps["reference_exact_occurrence_count"],
            "reference_class": ps["reference_class"], "acquisition_safe": ps["acquisition_safe"],
            "present_in_final_terminal_codebook": ps["present_in_final_terminal_codebook"],
            "internal_only_query_unit": int(ps["acquisition_safe"] and probe not in terminal_probes),
        })
        source_audit.append(row)
        if int(ps["acquisition_safe"]) == 1 and match:
            safe_by_te[s["te_id"]].append(row)
    full = [r for r in source_audit if int(r["acquisition_safe"]) == 1 and int(r["source_exact_position_match"]) == 1]
    internal_only = [r for r in full if int(r["internal_only_query_unit"]) == 1]

    te_rows = []
    for te in sorted(family_by_te):
        xs = safe_by_te.get(te, [])
        internal_units = {r["canonical_strand_probe"] for r in xs if int(r["internal_only_query_unit"]) == 1}
        te_rows.append({
            "te_id": te, "family": family_by_te[te], "te_length": te_len.get(te, 0),
            "terminal_codebook_units": len({r["canonical_strand_probe"] for r in terminal_rows if r["te_id"] == te}),
            "safe_internal_units": len({r["canonical_strand_probe"] for r in xs}),
            "safe_internal_only_units": len(internal_units),
            "safe_internal_only_read_length_blocks": len({int(r["read_length_block_125"]) for r in xs
                                                          if int(r["internal_only_query_unit"]) == 1}),
            "terminal_uncovered_record": int(te not in terminal_tes),
            "terminal_uncovered_but_internal_safe": int(te not in terminal_tes and bool(internal_units)),
        })
    fam_rows = []
    for fam in sorted(set(family_by_te.values())):
        trs = [r for r in te_rows if r["family"] == fam]
        fam_rows.append({
            "family": fam, "te_record_count": len(trs),
            "terminal_codebook_records": sum(r["terminal_codebook_units"] > 0 for r in trs),
            "records_with_safe_internal": sum(r["safe_internal_units"] > 0 for r in trs),
            "records_with_internal_only": sum(r["safe_internal_only_units"] > 0 for r in trs),
            "total_safe_internal_units": sum(r["safe_internal_units"] for r in trs),
            "total_internal_only_units": sum(r["safe_internal_only_units"] for r in trs),
            "terminal_uncovered_family": int(fam not in terminal_families),
            "terminal_uncovered_but_internal_safe": int(fam not in terminal_families
                                                        and any(r["safe_internal_only_units"] > 0 for r in trs)),
        })
    counts = {
        "frozen_clean_internal21_candidates": len(candidates),
        "internal_source_rows": len(source_rows),
        "exhaustive_source_family_safe_candidates": len(family_safe),
        "reference_safe_internal_query_units": len(safe_units),
        "safe_internal_only_query_units": len(safe_units - terminal_probes),
        "safe_internal_source_rows": len(full),
        "safe_internal_only_source_rows": len(internal_only),
    }
    return InternalResult(
        tsv_bytes(candidate_audit, style="internal"), tsv_bytes(source_audit, style="internal"),
        tsv_bytes(full, style="internal"), tsv_bytes(internal_only, style="internal"),
        tsv_bytes(te_rows, style="internal"), tsv_bytes(fam_rows, style="internal"), counts,
    )
