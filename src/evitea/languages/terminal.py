"""Final safe terminal codebook in three steps: rebuild, mixed-k gap fill, finalizer.

* :func:`rebuild_terminal21`: clean terminal k=21 anchor catalogue
  candidates; exhaustive TE-library family compatibility;
  RepeatMasker-aware reference safety; up to four non-overlapping probes per
  TE record, two per side first.
* :func:`fill_underpowered`: for records with fewer than two selected
  probes, terminal anchor catalogue candidates at k in {21, 25, 31, 35};
  cardinality/side-balanced interval selection.
* :func:`finalize_terminal`: deterministic repair of duplicate canonical
  probes within a record.

The steps exchange TSV bytes exactly as the original scripts exchanged files,
so every intermediate table can be compared byte for byte with the original.
No sample data are read.
"""

from __future__ import annotations

import bisect
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from .common import (
    build_interval_index,
    encode,
    extract_te_id_gapfill,
    extract_te_id_rebuild,
    load_repeatmasker_gapfill,
    load_repeatmasker_rebuild,
    oriented_code_map,
    overlap_query,
    parse_family_map_gapfill,
    parse_family_map_rebuild,
    parse_fasta,
    read_tsv_bytes,
    reference_hit_class,
    revcomp,
    scan_exact,
    tsv_bytes,
)

__all__ = [
    "CODEBOOK_FIELDS",
    "RebuildResult",
    "GapFillResult",
    "FinalTerminalResult",
    "rebuild_terminal21",
    "fill_underpowered",
    "finalize_terminal",
]

K = 21
ALLOWED_K = (21, 25, 31, 35)
REBUILD_CODE_FIELDS = (
    "family", "te_id", "probe_rank_within_te_record", "terminal_side", "canonical_strand_probe",
    "query_sequences_both_strands", "source_start", "source_end", "terminal_distance_bp",
    "positionally_interpretable", "source_te_distinct_occurrence_count",
    "reference_exact_occurrence_count", "reference_class", "library_role",
)
CODEBOOK_FIELDS = (
    "family", "te_id", "probe_rank_within_te_record", "terminal_side", "k",
    "canonical_strand_probe", "query_sequences_both_strands",
    "source_start", "source_end", "terminal_distance_bp",
    "positionally_interpretable", "source_te_distinct_occurrence_count",
    "reference_exact_occurrence_count", "reference_class", "library_role",
)


def _g1a_family_by_te(provenance: Sequence[Mapping[str, str]], *, terminal21_only: bool) -> dict[str, str]:
    family_by_te: dict[str, str] = {}
    for row in provenance:
        if terminal21_only and (int(row["k"]) != K or row["region"] not in {"5p_terminal", "3p_terminal"}):
            continue
        old = family_by_te.setdefault(row["te_id"], row["family"])
        if old != row["family"]:
            raise RuntimeError(f"Inconsistent anchor family for {row['te_id']}: {old} vs {row['family']}")
    return family_by_te


# ---------------------------------------------------------------------------
# Rebuild the safe terminal k=21 codebook
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class RebuildResult:
    candidate_audit: bytes
    source_audit: bytes
    codebook: bytes
    te_record_summary: bytes
    family_summary: bytes
    counts: dict[str, int] = field(default_factory=dict)


def _choose_record_codebook(candidates: list[dict]) -> list[dict]:
    def key(r):
        return (-int(r["positionally_interpretable"]), int(r["reference_exact_occurrence_count"]),
                int(r["terminal_distance_bp"]), r["canonical_strand_probe"])

    def nonoverlap(start, end, chosen):
        return all(not (max(start, r["source_start"]) < min(end, r["source_end"])) for r in chosen)

    by_side = {side: sorted([r for r in candidates if r["terminal_side"] == side], key=key) for side in ("5p", "3p")}
    chosen: list[dict] = []
    for side in ("5p", "3p"):
        for r in by_side[side]:
            if sum(x["terminal_side"] == side for x in chosen) >= 2:
                break
            if nonoverlap(r["source_start"], r["source_end"], chosen):
                chosen.append(r)
    for r in sorted(candidates, key=key):
        if len(chosen) >= 4:
            break
        if r in chosen:
            continue
        if nonoverlap(r["source_start"], r["source_end"], chosen):
            chosen.append(r)
    chosen = sorted(chosen, key=lambda r: (r["terminal_side"], r["source_start"], r["canonical_strand_probe"]))
    for rank, r in enumerate(chosen, 1):
        r["probe_rank_within_te_record"] = rank
    return chosen


def rebuild_terminal21(
    *,
    catalogue: Sequence[Mapping[str, str]],
    provenance: Sequence[Mapping[str, str]],
    te_fasta: str | Path,
    family_map: str | Path,
    reference_fasta: str | Path,
    repeatmasker_bed: str | Path,
    progress=None,
) -> RebuildResult:
    say = progress or (lambda message: None)
    cat_by_qid: dict[str, Mapping[str, str]] = {}
    canonical_candidates: set[str] = set()
    for r in catalogue:
        if int(r["k"]) != K or int(r["has_terminal_provenance"]) != 1 or int(r["low_complexity"]) != 0:
            continue
        cat_by_qid[r["query_id"]] = r
        canonical_candidates.add(r["canonical_strand_key"])

    source_rows: list[dict] = []
    family_by_te: dict[str, str] = {}
    for p in provenance:
        if int(p["k"]) != K or p["region"] not in {"5p_terminal", "3p_terminal"}:
            continue
        c = cat_by_qid.get(p["query_id"])
        if c is None:
            continue
        old = family_by_te.setdefault(p["te_id"], p["family"])
        if old != p["family"]:
            raise RuntimeError(f"Inconsistent anchor family for {p['te_id']}: {old} vs {p['family']}")
        source_rows.append({
            "canonical_strand_probe": c["canonical_strand_key"], "te_id": p["te_id"], "family": p["family"],
            "terminal_side": "5p" if p["region"] == "5p_terminal" else "3p",
            "source_start": int(p["start"]), "source_end": int(p["start"]) + K,
        })
    seen, dedup = set(), []
    for r in source_rows:
        key = (r["canonical_strand_probe"], r["te_id"], r["family"], r["terminal_side"], r["source_start"], r["source_end"])
        if key not in seen:
            seen.add(key)
            dedup.append(r)
    source_rows = dedup
    source_by_probe: dict[str, list[dict]] = defaultdict(list)
    for r in source_rows:
        source_by_probe[r["canonical_strand_probe"]].append(r)
    say(f"terminal 21-mer candidates={len(canonical_candidates):,}; records={len(family_by_te)}")

    qmap = oriented_code_map(canonical_candidates)
    te_hits: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    te_len: dict[str, int] = {}
    for header, seq in parse_fasta(te_fasta):
        te = extract_te_id_rebuild(header)
        te_len[te] = len(seq)
        family = family_by_te.get(te, te)
        for start, end, probe in scan_exact(seq, qmap, K):
            te_hits[probe].append((te, family, start, end))
    missing = sorted(set(family_by_te) - set(te_len))
    if missing:
        raise RuntimeError(f"TE library missing anchor catalogue records: {missing[:20]}")

    rm_index = build_interval_index(load_repeatmasker_rebuild(repeatmasker_bed, parse_family_map_rebuild(family_map)))
    context: dict[str, Counter] = defaultdict(Counter)
    ref_hits: Counter = Counter()
    priority = ("OTHER_FAMILY_ANNOTATED_TE", "NON_TE_OR_UNANNOTATED_HOST", "UNMAPPED_ANNOTATED_TE",
                "TE_BOUNDARY_OR_PARTIAL", "MIXED_REFERENCE_CONTEXT")
    for header, seq in parse_fasta(reference_fasta):
        chrom = header.split()[0]
        for start, end, probe in scan_exact(seq, qmap, K):
            ref_hits[probe] += 1
            labels = sorted({(x["family"], x["te_id"]) for x in source_by_probe.get(probe, [])})
            overlaps = overlap_query(rm_index, chrom, start, end)
            per_label = [reference_hit_class(fam, te, start, end, overlaps) for fam, te in labels]
            if not per_label:
                cls = "NO_SOURCE_PROVENANCE"
            elif all(x == "TARGET_CONTEXT_ANNOTATED_TE" for x in per_label):
                cls = "TARGET_CONTEXT_ANNOTATED_TE"
            else:
                cls = next((x for x in priority if x in per_label), per_label[0])
            context[probe][cls] += 1
        say(f"reference scanned: {chrom}")

    candidate_rows, probe_global = [], {}
    for probe in sorted(canonical_candidates):
        hits = te_hits.get(probe, [])
        hit_families = sorted({fam for _, fam, _, _ in hits})
        hit_tes = sorted({te for te, _, _, _ in hits})
        sources = source_by_probe.get(probe, [])
        source_families = sorted({x["family"] for x in sources})
        single = len(hit_families) == 1
        consistent = len(source_families) == 1 and single and source_families[0] == hit_families[0]
        counts = context.get(probe, Counter())
        total = ref_hits[probe]
        if total == 0:
            ref_safe, ref_class = True, "REFERENCE_ABSENT"
        else:
            ref_safe = total - counts["TARGET_CONTEXT_ANNOTATED_TE"] == 0
            ref_class = "TARGET_CONTEXT_ANNOTATED_TE_ONLY" if ref_safe else "UNSAFE_OR_MIXED_REFERENCE_CONTEXT"
        probe_global[probe] = {
            "canonical_strand_probe": probe, "source_families": source_families,
            "source_te_ids": sorted({x["te_id"] for x in sources}),
            "exhaustive_hit_families": hit_families, "exhaustive_hit_te_ids": hit_tes,
            "exhaustive_single_family": int(single), "source_family_consistent": int(consistent),
            "reference_exact_occurrence_count": total,
            "reference_target_context_hits": counts["TARGET_CONTEXT_ANNOTATED_TE"],
            "reference_unsafe_context_hits": total - counts["TARGET_CONTEXT_ANNOTATED_TE"],
            "reference_class": ref_class, "acquisition_safe": int(bool(consistent and ref_safe)),
        }
        candidate_rows.append(probe_global[probe])

    source_audit, safe_by_te = [], defaultdict(list)
    for s in source_rows:
        g = probe_global[s["canonical_strand_probe"]]
        positions = sorted({(h[2], h[3]) for h in te_hits[s["canonical_strand_probe"]] if h[0] == s["te_id"]})
        length = te_len[s["te_id"]]
        distance = s["source_start"] if s["terminal_side"] == "5p" else length - s["source_end"]
        opposite = 0
        for st, en in positions:
            if (st, en) == (s["source_start"], s["source_end"]):
                continue
            if s["terminal_side"] == "5p":
                opposite += length - en <= distance
            else:
                opposite += st <= distance
        match = (s["source_start"], s["source_end"]) in positions
        positional = int(g["acquisition_safe"]) == 1 and match and len(positions) == 1 and opposite == 0
        row = dict(s)
        row.update({
            "te_record_length": length, "terminal_distance_bp": distance,
            "source_exact_position_match": int(match), "source_te_distinct_occurrence_count": len(positions),
            "opposite_terminal_recurrence_count": opposite,
            "exhaustive_hit_family_count": len(g["exhaustive_hit_families"]),
            "exhaustive_hit_te_count": len(g["exhaustive_hit_te_ids"]),
            "reference_exact_occurrence_count": g["reference_exact_occurrence_count"],
            "reference_class": g["reference_class"], "acquisition_safe": g["acquisition_safe"],
            "positionally_interpretable": int(positional),
        })
        source_audit.append(row)
        if int(g["acquisition_safe"]) == 1 and match:
            safe_by_te[s["te_id"]].append(row)

    selected = []
    for te in sorted(family_by_te):
        unique = {}
        for r in safe_by_te.get(te, []):
            unique[(r["canonical_strand_probe"], r["terminal_side"], r["source_start"], r["source_end"])] = dict(r)
        selected.extend(_choose_record_codebook(list(unique.values())))
    selected.sort(key=lambda r: (r["family"], r["te_id"], r["probe_rank_within_te_record"]))
    for r in selected:
        r["query_sequences_both_strands"] = f"{r['canonical_strand_probe']};{revcomp(r['canonical_strand_probe'])}"
        r["library_role"] = "REBUILT_SAFE_TERMINAL21_ACQUISITION"

    by_te = defaultdict(list)
    for r in selected:
        by_te[r["te_id"]].append(r)
    te_rows = []
    for te in sorted(family_by_te):
        safe = safe_by_te.get(te, [])
        chosen = by_te.get(te, [])
        te_rows.append({
            "te_id": te, "family": family_by_te[te],
            "available_safe_terminal21_probes": len({r["canonical_strand_probe"] for r in safe}),
            "available_safe_5p": len({r["canonical_strand_probe"] for r in safe if r["terminal_side"] == "5p"}),
            "available_safe_3p": len({r["canonical_strand_probe"] for r in safe if r["terminal_side"] == "3p"}),
            "selected_probe_count": len(chosen),
            "selected_5p": sum(r["terminal_side"] == "5p" for r in chosen),
            "selected_3p": sum(r["terminal_side"] == "3p" for r in chosen),
            "selected_positionally_interpretable": sum(int(r["positionally_interpretable"]) for r in chosen),
        })
    family_to_tes = defaultdict(list)
    for te, fam in family_by_te.items():
        family_to_tes[fam].append(te)
    fam_rows = []
    for fam in sorted(family_to_tes):
        tes = sorted(family_to_tes[fam])
        rows = [r for r in te_rows if r["te_id"] in tes]
        fam_rows.append({
            "family": fam, "te_record_count": len(tes),
            "records_with_ge1_safe_terminal21": sum(r["available_safe_terminal21_probes"] >= 1 for r in rows),
            "records_with_ge2_safe_terminal21": sum(r["available_safe_terminal21_probes"] >= 2 for r in rows),
            "records_with_ge4_safe_terminal21": sum(r["available_safe_terminal21_probes"] >= 4 for r in rows),
            "family_total_available_safe_terminal21": sum(r["available_safe_terminal21_probes"] for r in rows),
            "selected_probe_count": sum(r["selected_probe_count"] for r in rows),
            "selected_positionally_interpretable": sum(r["selected_positionally_interpretable"] for r in rows),
        })
    counts = {
        "clean_terminal21_canonical_probes": len(canonical_candidates),
        "te_records": len(family_by_te),
        "source_family_consistent_probes": sum(int(r["source_family_consistent"]) for r in candidate_rows),
        "reference_context_safe_probes": sum(r["reference_class"] in {"REFERENCE_ABSENT", "TARGET_CONTEXT_ANNOTATED_TE_ONLY"}
                                             for r in candidate_rows),
        "acquisition_safe_probes": sum(int(r["acquisition_safe"]) for r in candidate_rows),
        "selected_rows": len(selected),
    }
    return RebuildResult(
        tsv_bytes(candidate_rows, list(candidate_rows[0])), tsv_bytes(source_audit, list(source_audit[0])),
        tsv_bytes(selected, REBUILD_CODE_FIELDS), tsv_bytes(te_rows, list(te_rows[0])),
        tsv_bytes(fam_rows, list(fam_rows[0])), counts,
    )


# ---------------------------------------------------------------------------
# Mixed-k gap fill of underpowered records
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class GapFillResult:
    source_audit: bytes
    target_summary: bytes
    hybrid_codebook: bytes
    counts: dict[str, int] = field(default_factory=dict)


def _select_best_nonoverlapping(cands: list[dict], max_count: int = 4) -> list[dict]:
    if not cands:
        return []
    cands = sorted(cands, key=lambda r: (r["source_end"], r["source_start"], r["k"], r["canonical_strand_probe"]))
    ends = [r["source_end"] for r in cands]
    prev = [bisect.bisect_right(ends, r["source_start"], 0, i) - 1 for i, r in enumerate(cands)]
    dp: list[dict] = [dict() for _ in range(len(cands) + 1)]
    dp[0][(0, 0, 0)] = ((0, 0, 0, 0), tuple())

    def additive(r):
        return (int(r["positionally_interpretable"]), -int(r["k"]), -int(r["reference_exact_occurrence_count"]),
                -int(r["terminal_distance_bp"]))

    def better(x, y):
        if y is None:
            return True
        if x[0] != y[0]:
            return x[0] > y[0]
        sig = lambda s: tuple((cands[j]["k"], cands[j]["canonical_strand_probe"], cands[j]["source_start"]) for j in s[1])  # noqa: E731
        return sig(x) < sig(y)

    for i in range(1, len(cands) + 1):
        r = cands[i - 1]
        dp[i] = dict(dp[i - 1])
        for (n, c5, c3), state in dp[prev[i - 1] + 1].items():
            if n >= max_count:
                continue
            key = (n + 1, c5 + (r["terminal_side"] == "5p"), c3 + (r["terminal_side"] == "3p"))
            candidate = (tuple(a + b for a, b in zip(state[0], additive(r))), state[1] + (i - 1,))
            if better(candidate, dp[i].get(key)):
                dp[i][key] = candidate
    best_rank = best_state = None
    for (n, c5, c3), state in dp[-1].items():
        add = state[0]
        rank = (n, -abs(c5 - c3), min(c5, c3), add[0], add[1], add[2], add[3])
        if best_rank is None or rank > best_rank:
            best_rank, best_state = rank, state
    selected = [dict(cands[i]) for i in best_state[1]]
    selected.sort(key=lambda r: (0 if r["terminal_side"] == "5p" else 1, r["source_start"], r["k"], r["canonical_strand_probe"]))
    for rank, r in enumerate(selected, 1):
        r["probe_rank_within_te_record"] = rank
    return selected


def fill_underpowered(
    *,
    rebuilt_codebook: bytes,
    rebuilt_te_summary: bytes,
    catalogue: Sequence[Mapping[str, str]],
    provenance: Sequence[Mapping[str, str]],
    te_fasta: str | Path,
    family_map: str | Path,
    reference_fasta: str | Path,
    repeatmasker_bed: str | Path,
    progress=None,
) -> GapFillResult:
    say = progress or (lambda message: None)
    prior_codebook = read_tsv_bytes(rebuilt_codebook)
    targets = [r for r in read_tsv_bytes(rebuilt_te_summary) if int(r["selected_probe_count"]) < 2]
    target_ids = {r["te_id"] for r in targets}
    if not target_ids:
        raise RuntimeError("No underpowered records found; nothing to do.")
    family_by_te = _g1a_family_by_te(provenance, terminal21_only=False)
    cat = {(r["query_id"], int(r["k"])): r for r in catalogue}
    source_rows, by_k, seen = [], defaultdict(set), set()
    for p in provenance:
        if p["te_id"] not in target_ids or p["region"] not in {"5p_terminal", "3p_terminal"}:
            continue
        if p.get("queryable", "0") != "1":
            continue
        k = int(p["k"])
        if k not in ALLOWED_K:
            continue
        c = cat.get((p["query_id"], k))
        if c is None:
            raise RuntimeError(f"Missing catalogue row for {p['query_id']} k={k}")
        if int(c["low_complexity"]) != 0 or int(c["has_terminal_provenance"]) != 1:
            continue
        probe, start = c["canonical_strand_key"], int(p["start"])
        side = "5p" if p["region"] == "5p_terminal" else "3p"
        key = (probe, p["te_id"], p["family"], side, start, start + k, k)
        if key in seen:
            continue
        seen.add(key)
        source_rows.append({"canonical_strand_probe": probe, "k": k, "te_id": p["te_id"], "family": p["family"],
                            "terminal_side": side, "source_start": start, "source_end": start + k})
        by_k[k].add(probe)
    say(f"gap-fill targets={sorted(target_ids)}; candidates by k={ {k: len(v) for k, v in by_k.items()} }")
    qmaps = {}
    for k, probes in by_k.items():
        qmap = {}
        for probe in probes:
            if len(probe) != k:
                raise RuntimeError(f"Probe length mismatch: k={k}")
            for query in (probe, revcomp(probe)):
                code = encode(query)
                if qmap.get(code, probe) != probe:
                    raise RuntimeError(f"Oriented query collision at k={k}")
                qmap[code] = probe
        qmaps[k] = qmap

    def scan(seq):
        for k, qmap in qmaps.items():
            for start, end, probe in scan_exact(seq, qmap, k):
                yield k, start, end, probe

    te_hits, te_len, headers = defaultdict(list), {}, set()
    for header, seq in parse_fasta(te_fasta):
        te = extract_te_id_gapfill(header)
        headers.add(te)
        te_len[te] = len(seq)
        family = family_by_te.get(te, te)
        for k, start, end, probe in scan(seq):
            te_hits[(k, probe)].append((te, family, start, end))
    if target_ids - headers:
        raise RuntimeError(f"Target records missing from TE FASTA: {sorted(target_ids - headers)}")

    rm_index = build_interval_index(load_repeatmasker_gapfill(repeatmasker_bed, parse_family_map_gapfill(family_map)))
    source_by_probe = defaultdict(list)
    for r in source_rows:
        source_by_probe[(r["k"], r["canonical_strand_probe"])].append(r)
    context, ref_count = defaultdict(Counter), Counter()
    for header, seq in parse_fasta(reference_fasta):
        chrom = header.split()[0]
        for k, start, end, probe in scan(seq):
            overlaps = overlap_query(rm_index, chrom, start, end)
            for src in source_by_probe[(k, probe)]:
                skey = (src["te_id"], k, probe, src["terminal_side"], src["source_start"], src["source_end"])
                context[skey][reference_hit_class(src["family"], src["te_id"], start, end, overlaps)] += 1
                ref_count[skey] += 1

    source_audit, safe_by_te = [], defaultdict(list)
    for src in source_rows:
        k, probe, te, fam = src["k"], src["canonical_strand_probe"], src["te_id"], src["family"]
        key = (te, k, probe, src["terminal_side"], src["source_start"], src["source_end"])
        hits = te_hits.get((k, probe), [])
        hit_families = sorted({h[1] for h in hits})
        family_safe = len(hit_families) == 1 and hit_families[0] == fam
        positions = sorted({(h[2], h[3]) for h in hits if h[0] == te})
        match = (src["source_start"], src["source_end"]) in positions
        length = te_len[te]
        distance = src["source_start"] if src["terminal_side"] == "5p" else length - src["source_end"]
        opposite = 0
        for st, en in positions:
            if (st, en) == (src["source_start"], src["source_end"]):
                continue
            opposite += (length - en <= distance) if src["terminal_side"] == "5p" else (st <= distance)
        counts = context.get(key, Counter())
        n_ref = ref_count[key]
        unsafe = n_ref - counts["TARGET_CONTEXT_ANNOTATED_TE"]
        ref_safe = unsafe == 0
        ref_class = "REFERENCE_ABSENT" if n_ref == 0 else (
            "TARGET_CONTEXT_ANNOTATED_TE_ONLY" if ref_safe else "UNSAFE_OR_MIXED_REFERENCE_CONTEXT")
        safe = family_safe and ref_safe and match
        row = dict(src)
        row.update({
            "te_record_length": length, "terminal_distance_bp": distance,
            "exhaustive_hit_family_count": len(hit_families), "exhaustive_hit_families": hit_families,
            "exhaustive_hit_te_count": len({h[0] for h in hits}),
            "source_exact_position_match": int(match), "source_te_distinct_occurrence_count": len(positions),
            "opposite_terminal_recurrence_count": opposite, "reference_exact_occurrence_count": n_ref,
            "reference_target_context_hits": counts["TARGET_CONTEXT_ANNOTATED_TE"],
            "reference_unsafe_context_hits": unsafe, "reference_class": ref_class,
            "family_safe": int(family_safe), "reference_safe": int(ref_safe), "acquisition_safe": int(safe),
            "positionally_interpretable": int(safe and len(positions) == 1 and opposite == 0),
        })
        source_audit.append(row)
        if safe:
            safe_by_te[te].append(row)
    source_audit.sort(key=lambda r: (r["te_id"], r["k"], r["terminal_side"], r["source_start"], r["canonical_strand_probe"]))

    replacement, target_summary = [], []
    for t in sorted(targets, key=lambda r: r["te_id"]):
        unique = {}
        for r in safe_by_te.get(t["te_id"], []):
            unique[(r["k"], r["canonical_strand_probe"], r["terminal_side"], r["source_start"], r["source_end"])] = dict(r)
        safe = list(unique.values())
        chosen = _select_best_nonoverlapping(safe, max_count=4)
        for r in chosen:
            r["query_sequences_both_strands"] = f"{r['canonical_strand_probe']};{revcomp(r['canonical_strand_probe'])}"
            r["library_role"] = "SAFE_TERMINAL_MIXED_K_GAP_FILL"
        replacement.extend(chosen)
        target_summary.append({
            "te_id": t["te_id"], "family": t["family"], "prior_selected_probe_count": int(t["selected_probe_count"]),
            "safe_candidates_total": len(safe),
            **{f"safe_k{k}": sum(int(r["k"]) == k for r in safe) for k in ALLOWED_K},
            "selected_probe_count": len(chosen),
            **{f"selected_k{k}": sum(int(r["k"]) == k for r in chosen) for k in ALLOWED_K},
            "selected_5p": sum(r["terminal_side"] == "5p" for r in chosen),
            "selected_3p": sum(r["terminal_side"] == "3p" for r in chosen),
            "selected_positionally_interpretable": sum(int(r["positionally_interpretable"]) for r in chosen),
        })
    hybrid = []
    for r in prior_codebook:
        if r["te_id"] in target_ids:
            continue
        row = dict(r)
        row["k"] = 21
        row["library_role"] = row.get("library_role", "REBUILT_SAFE_TERMINAL21_ACQUISITION")
        hybrid.append(row)
    hybrid.extend(replacement)
    hybrid.sort(key=lambda r: (r["family"], r["te_id"], int(r["probe_rank_within_te_record"]), int(r["k"])))
    return GapFillResult(
        tsv_bytes(source_audit, list(source_audit[0])), tsv_bytes(target_summary, list(target_summary[0])),
        tsv_bytes(hybrid, CODEBOOK_FIELDS),
        {"targets": len(target_ids), "replacement_rows": len(replacement), "hybrid_rows": len(hybrid)},
    )


# ---------------------------------------------------------------------------
# Finalizer (duplicate-sequence repair)
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class FinalTerminalResult:
    codebook: bytes
    record_status: bytes
    repair_log: bytes
    counts: dict[str, object] = field(default_factory=dict)


def _overlaps(a, b) -> bool:
    return max(int(a["source_start"]), int(b["source_start"])) < min(int(a["source_end"]), int(b["source_end"]))


def _choose_repair(cands, target_n: int):
    def score(r):
        return (-int(r.get("positionally_interpretable", 0)), int(r.get("reference_exact_occurrence_count", 10**18)),
                int(r.get("terminal_distance_bp", 10**18)), int(r.get("k", 10**18)), r["canonical_strand_probe"],
                int(r["source_start"]))

    by_side = {side: sorted([dict(r) for r in cands if r["terminal_side"] == side], key=score) for side in ("5p", "3p")}
    chosen, used = [], set()

    def admissible(r):
        return r["canonical_strand_probe"] not in used and not any(_overlaps(r, x) for x in chosen)

    for side in ("5p", "3p"):
        for r in by_side[side]:
            if len([x for x in chosen if x["terminal_side"] == side]) >= 2:
                break
            if admissible(r):
                chosen.append(r)
                used.add(r["canonical_strand_probe"])
    for r in sorted([dict(r) for r in cands], key=score):
        if len(chosen) >= target_n:
            break
        if admissible(r):
            chosen.append(r)
            used.add(r["canonical_strand_probe"])
    chosen.sort(key=lambda r: (r["terminal_side"], int(r["source_start"]), int(r["k"]), r["canonical_strand_probe"]))
    for rank, r in enumerate(chosen, 1):
        q = r["canonical_strand_probe"]
        r["probe_rank_within_te_record"] = str(rank)
        r["query_sequences_both_strands"] = q if revcomp(q) == q else f"{q};{revcomp(q)}"
        r["library_role"] = "FINAL_SAFE_TERMINAL_ACQUISITION"
    return chosen


def finalize_terminal(
    *,
    hybrid_codebook: bytes,
    rebuilt_source_audit: bytes,
    gapfill_source_audit: bytes,
    rebuilt_te_summary: bytes,
) -> FinalTerminalResult:
    hybrid = read_tsv_bytes(hybrid_codebook)
    source21 = read_tsv_bytes(rebuilt_source_audit)
    if any(r.get("k", "") not in {"", "21"} for r in source21):
        raise RuntimeError("Terminal-21 source audit contains a non-21 k value.")
    for r in source21:
        r["k"] = "21"
    source_rows = source21 + read_tsv_bytes(gapfill_source_audit)
    all_te = {r["te_id"]: r["family"] for r in read_tsv_bytes(rebuilt_te_summary)}
    by_te = defaultdict(list)
    for r in hybrid:
        by_te[r["te_id"]].append(dict(r))
    for te, rows in by_te.items():
        for i, a in enumerate(rows):
            for b in rows[i + 1:]:
                if _overlaps(a, b):
                    raise RuntimeError(f"Selected source intervals overlap for {te}")
    dup_records = [te for te, rows in by_te.items()
                   if len({r["canonical_strand_probe"] for r in rows}) != len(rows)]
    safe_by_te, seen = defaultdict(list), set()
    for r in source_rows:
        if str(r.get("acquisition_safe", "0")) != "1" or str(r.get("source_exact_position_match", "1")) != "1":
            continue
        key = (r["te_id"], r["canonical_strand_probe"], r["terminal_side"], int(r["source_start"]),
               int(r["source_end"]), int(r["k"]))
        if key not in seen:
            seen.add(key)
            safe_by_te[r["te_id"]].append(dict(r))
    repaired, repair_log = [], []
    for te in sorted(all_te):
        rows = by_te.get(te, [])
        if te not in dup_records:
            for r in rows:
                r["library_role"] = "FINAL_SAFE_TERMINAL_ACQUISITION"
                repaired.append(r)
            continue
        new_rows = _choose_repair(safe_by_te[te], len(rows))
        if len(new_rows) != len(rows) or len({r["canonical_strand_probe"] for r in new_rows}) != len(rows):
            raise RuntimeError(f"Could not repair duplicate codebook for {te}")
        repaired.extend(new_rows)
        repair_log.append({
            "te_id": te, "family": all_te[te], "reason": "duplicate_canonical_probe_within_record",
            "old_probe_count": len(rows), "new_probe_count": len(new_rows),
            "old_probes": ";".join(r["canonical_strand_probe"] for r in rows),
            "new_probes": ";".join(r["canonical_strand_probe"] for r in new_rows),
        })
    repaired.sort(key=lambda r: (r["family"], r["te_id"], int(r["probe_rank_within_te_record"])))
    final_by_te = defaultdict(list)
    for r in repaired:
        final_by_te[r["te_id"]].append(r)
        k = int(r["k"])
        if k not in set(ALLOWED_K) or len(r["canonical_strand_probe"]) != k or set(r["canonical_strand_probe"]) - set("ACGT"):
            raise RuntimeError(f"Invalid canonical probe: {r}")
    for te, rows in final_by_te.items():
        if len({r["canonical_strand_probe"] for r in rows}) != len(rows):
            raise RuntimeError(f"Within-record duplicate remains for {te}")
    status_rows = []
    for te in sorted(all_te):
        n = len(final_by_te.get(te, []))
        status = "UNSUPPORTED_TERMINAL_SAFE" if n == 0 else "SINGLE_INDEPENDENT_BLOCK" if n == 1 else "MULTI_INDEPENDENT_BLOCK"
        status_rows.append({"te_id": te, "family": all_te[te], "selected_probe_count": n, "acquisition_status": status})
    return FinalTerminalResult(
        tsv_bytes(repaired, CODEBOOK_FIELDS, style="raw"),
        tsv_bytes(status_rows, ("te_id", "family", "selected_probe_count", "acquisition_status"), style="raw"),
        tsv_bytes(repair_log, ("te_id", "family", "reason", "old_probe_count", "new_probe_count", "old_probes", "new_probes"),
                  style="raw"),
        {"rows": len(repaired), "status_counts": dict(Counter(r["acquisition_status"] for r in status_rows)),
         "records_repaired_for_duplicate_sequence": dup_records},
    )
