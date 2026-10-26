"""Paired-end mate-anchored search, and strict reference-discriminative pair evidence.

Orientation-corrected mate-anchored search
    Candidate loci not within 50 bp of a same-family split-junction call are
    examined. For each candidate pair both TE-mate/host-mate directions are
    tested. The host mate needs one unique best local placement (<= 5 edits,
    midpoint within 250 bp of the candidate envelope). Side is geometric:
    LEFT = forward host left of the candidate midpoint, RIGHT = reverse host
    right of it. The TE mate needs a candidate-family-rich window whose unique
    best TE family matches the candidate (alias-compatible) and that is not
    explained by the local reference. One physical pair supports one
    (locus, side); an equal-strength conflict gives no winner. Loci with >= 2
    LEFT and >= 2 RIGHT pairs become range calls.

Reference-discriminative pair evidence
    One-sided non-passing loci are routed with bounded exact msBWT host seeds
    on the missing flank, locally prefiltered, globally host-mapped with
    Bowtie2, and rejected when an ordinary inward concordant pair within the
    fragment-size 95th percentile (the sample's own; 494 bp by default)
    explains the pair. The remaining pairs are TE-confirmed with
    strand-symmetric local-reference rejection. Only retained per-pair
    winners are emitted; they feed junction completion as strict side support.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from .alignment import hw_distance, hw_locations
from .candidates import CandidateLocus
from .events import IncumbentCall, IncumbentCallClass
from .external import HOST_GLOBAL_OPTIONS, cigar_ref_span, parse_sam_records, run_bowtie2, write_fasta
from .sequence import (
    distance_to_interval,
    interval_distance_point,
    interval_gap,
    interval_hits,
    kmer_hit_prefix,
    revcomp,
    shannon_entropy_symbols,
)
from .te_library import LocalTEIndex, family_best_classify, normalize_family

SCOPE_RADIUS = 800
CANDIDATE_SLOP = 250
CORE_GUARD_BP = 50
HOST_MATE_MAX_EDITS = 5
TE_K = 21
TE_WINDOW_BP = 45
TE_WINDOW_STEP_BP = 5
MIN_TE_WINDOW_HITS = 3
MIN_TE_WINDOW_FRACTION = 0.12
MAX_TE_WINDOWS_PER_MATE = 6
TE_LOCAL_REF_MAX_EDITS = 3
TE_FULL_READ_REF_MAX_EDITS = 8
MIN_SIDE_FRAGMENTS = 2
# Default 95th percentile of the fragment size (bp); a run passes the sample's own (``sample_geometry``).
FRAGMENT_Q95 = 494.0
EXISTING_GUARD_BP = 50
HOST_SEED_K = 25
HOST_SEED_STRIDE = 10
HOST_SEED_OUTER_PAD = 250
HOST_SEED_MIN_ENTROPY = 1.20
HOST_SEED_ROW_CAP = 64
MAX_RECRUITED_PAIRS_PER_LOCUS = 1024
GLOBAL_PLACEMENT_CAP = 20


def quantize(value: float) -> float:
    """Round to one decimal."""

    return float(f"{value:.1f}")


# --------------------------------------------------------------------------
# Shared primitives
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class HostAnchor:
    view: str
    start: int
    end: int
    midpoint: float
    edits: int
    candidate_distance: float


def best_local_host_anchor(
    read: str, scope: str, scope_start: int, locus_start: float, locus_end: float
) -> tuple[HostAnchor | None, str]:
    candidates: list[tuple[int, str, int, int]] = []
    for view, sequence in (("FWD", read.upper()), ("RC", revcomp(read))):
        distance, locations = hw_locations(sequence, scope, HOST_MATE_MAX_EDITS)
        if distance is None:
            continue
        for start, end in locations:
            candidates.append((distance, view, scope_start + start, scope_start + end + 1))
    if not candidates:
        return None, "HOST_MATE_NO_LOCAL_MATCH"
    best = min(row[0] for row in candidates)
    by_interval: dict[tuple[int, int], tuple[int, str, int, int]] = {}
    for row in candidates:
        if row[0] != best:
            continue
        key = (row[2], row[3])
        old = by_interval.get(key)
        if old is None or row[1] < old[1]:
            by_interval[key] = row
    if len(by_interval) != 1:
        return None, "HOST_MATE_AMBIGUOUS_BEST_PLACEMENT"
    distance, view, start, end = next(iter(by_interval.values()))
    midpoint = (start + end) / 2.0
    candidate_distance = distance_to_interval(midpoint, locus_start, locus_end)
    if candidate_distance > CANDIDATE_SLOP:
        return None, "HOST_MATE_OUTSIDE_CANDIDATE_ENVELOPE"
    return HostAnchor(view, start, end, midpoint, distance, candidate_distance), "HOST_MATE_UNIQUE_LOCAL_ANCHOR"


def classify_side(host_view: str, host_midpoint: float, candidate_midpoint: float) -> str | None:
    if host_view == "FWD" and host_midpoint < candidate_midpoint:
        return "LEFT"
    if host_view == "RC" and host_midpoint > candidate_midpoint:
        return "RIGHT"
    return None


def candidate_te_windows(read: str, family_kmers: set[str]) -> list[tuple[int, int, str, int, int, float]]:
    sequence = read.upper()
    if len(sequence) < TE_WINDOW_BP or not family_kmers:
        return []
    prefix = kmer_hit_prefix(sequence, family_kmers, TE_K)
    raw: list[tuple[int, int, str, int, int, float]] = []
    for start in range(0, len(sequence) - TE_WINDOW_BP + 1, TE_WINDOW_STEP_BP):
        end = start + TE_WINDOW_BP
        hits, total = interval_hits(prefix, start, end, TE_K)
        fraction = hits / total if total else 0.0
        if hits >= MIN_TE_WINDOW_HITS and fraction >= MIN_TE_WINDOW_FRACTION:
            raw.append((start, end, sequence[start:end], hits, total, fraction))
    final_start = len(sequence) - TE_WINDOW_BP
    if final_start >= 0 and all(row[0] != final_start for row in raw):
        hits, total = interval_hits(prefix, final_start, len(sequence), TE_K)
        fraction = hits / total if total else 0.0
        if hits >= MIN_TE_WINDOW_HITS and fraction >= MIN_TE_WINDOW_FRACTION:
            raw.append((final_start, len(sequence), sequence[final_start:], hits, total, fraction))
    raw.sort(key=lambda row: (-row[3], -row[5], row[0]))
    return raw[:MAX_TE_WINDOWS_PER_MATE]


@dataclass(frozen=True, slots=True)
class TEWindow:
    start: int
    end: int
    sequence: str
    seed_hits: int
    seed_total: int
    seed_fraction: float
    family: str
    fbte_id: str
    te_edits: int
    te_orientation: str
    local_ref_edits: int | None
    full_read_ref_edits: int | None


def _record_conflict(first: TEWindow, second: TEWindow) -> bool:
    return (
        first.te_edits == second.te_edits
        and first.seed_hits == second.seed_hits
        and round(first.seed_fraction, 6) == round(second.seed_fraction, 6)
        and (first.family, first.fbte_id) != (second.family, second.fbte_id)
    )


def best_te_mate_window(
    read: str, candidate_family: str, family_kmers: set[str], te_index: LocalTEIndex, scope: str
) -> tuple[TEWindow | None, str]:
    """Confirm the TE mate (any local-reference explanation rejects it)."""

    windows = candidate_te_windows(read, family_kmers)
    if not windows:
        return None, "TE_MATE_NO_RICH_WINDOW"
    full_reference = hw_distance(read.upper(), scope, TE_FULL_READ_REF_MAX_EDITS)
    candidates: list[TEWindow] = []
    saw_other = saw_local = saw_none = False
    for start, end, sequence, hits, total, fraction in windows:
        hit, _ = family_best_classify(te_index.record_hits(sequence))
        if hit is None:
            saw_none = True
            continue
        if normalize_family(hit.family) != normalize_family(candidate_family):
            saw_other = True
            continue
        local = hw_distance(sequence, scope, TE_LOCAL_REF_MAX_EDITS)
        if local is not None:
            saw_local = True
            continue
        if full_reference is not None:
            saw_local = True
            continue
        candidates.append(TEWindow(start, end, sequence, hits, total, fraction, hit.family, hit.record_id,
                                   hit.edit_distance, hit.orientation, local, full_reference))
    if not candidates:
        if saw_local:
            return None, "TE_MATE_REFERENCE_COMPATIBLE"
        if saw_other:
            return None, "TE_MATE_WRONG_UNIQUE_BEST_FAMILY"
        if saw_none:
            return None, "TE_MATE_NO_UNIQUE_TE_MATCH"
        return None, "TE_MATE_REJECT"
    candidates.sort(key=lambda row: (row.te_edits, -row.seed_hits, -row.seed_fraction, row.start, row.fbte_id))
    if len(candidates) > 1 and _record_conflict(candidates[0], candidates[1]):
        return None, "TE_MATE_EQUAL_STRENGTH_RECORD_CONFLICT"
    return candidates[0], "TE_MATE_CERTIFIED"


def _best_symmetric_local(sequence: str, scope: str, max_edits: int) -> int | None:
    values = [
        value
        for value in (hw_distance(sequence.upper(), scope, max_edits), hw_distance(revcomp(sequence.upper()), scope, max_edits))
        if value is not None
    ]
    return min(values) if values else None


def best_te_window_symmetric(
    read: str, candidate_family: str, family_kmers: set[str], te_index: LocalTEIndex, scope: str
) -> tuple[TEWindow | None, str]:
    """Confirm the TE mate with strand-symmetric local-reference rejection."""

    windows = candidate_te_windows(read, family_kmers)
    if not windows:
        return None, "TE_MATE_NO_RICH_WINDOW"
    full_local = _best_symmetric_local(read, scope, TE_FULL_READ_REF_MAX_EDITS)
    if full_local is not None:
        return None, "TE_MATE_FULL_READ_LOCAL_REFERENCE_COMPATIBLE"
    candidates: list[TEWindow] = []
    saw_other = saw_local = saw_none = False
    for start, end, sequence, hits, total, fraction in windows:
        hit, _ = family_best_classify(te_index.record_hits(sequence))
        if hit is None:
            saw_none = True
            continue
        if normalize_family(hit.family) != normalize_family(candidate_family):
            saw_other = True
            continue
        local = _best_symmetric_local(sequence, scope, TE_LOCAL_REF_MAX_EDITS)
        if local is not None:
            saw_local = True
            continue
        candidates.append(TEWindow(start, end, sequence, hits, total, fraction, hit.family, hit.record_id,
                                   hit.edit_distance, hit.orientation, local, full_local))
    if not candidates:
        if saw_local:
            return None, "TE_MATE_WINDOW_LOCAL_REFERENCE_COMPATIBLE"
        if saw_other:
            return None, "TE_MATE_WRONG_UNIQUE_BEST_FAMILY"
        if saw_none:
            return None, "TE_MATE_NO_UNIQUE_TE_MATCH"
        return None, "TE_MATE_REJECT"
    candidates.sort(key=lambda row: (row.te_edits, -row.seed_hits, -round(row.seed_fraction, 6), row.start, row.fbte_id))
    if len(candidates) > 1 and _record_conflict(candidates[0], candidates[1]):
        return None, "TE_MATE_EQUAL_STRENGTH_RECORD_CONFLICT"
    return candidates[0], "TE_MATE_CERTIFIED"


# --------------------------------------------------------------------------
# Orientation-corrected mate-anchored search
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MateEvidence:
    locus_id: str
    chrom: str
    candidate_family: str
    pair_id: int
    side: str
    te_mate: int
    host_mate: int
    host_view: str
    host_start: int
    host_end: int
    host_midpoint: float
    host_edits: int
    host_candidate_distance: float
    te_seed_hits: int
    te_seed_fraction: float
    te_edits: int
    extra_strength: tuple = ()

    @property
    def strength_key(self) -> tuple:
        return (
            self.host_edits + self.te_edits,
            self.host_edits,
            self.te_edits,
            -self.te_seed_hits,
            -round(self.te_seed_fraction, 6),
            self.host_candidate_distance,
            *self.extra_strength,
        )

    @property
    def decision_key(self) -> tuple:
        return (*self.strength_key, self.candidate_family, self.chrom, self.host_midpoint,
                self.locus_id, self.side, self.te_mate)

    @property
    def kept_key(self) -> tuple:
        return (self.locus_id, self.side, self.te_mate, round(self.host_midpoint, 3))


def resolve_pair_evidence(rows: Sequence[MateEvidence]) -> tuple[MateEvidence | None, str]:
    if not rows:
        return None, "NO_EVIDENCE"
    by_semantic: dict[tuple[str, str], MateEvidence] = {}
    for row in rows:
        key = (row.locus_id, row.side)
        old = by_semantic.get(key)
        if old is None or row.decision_key < old.decision_key:
            by_semantic[key] = row
    ordered = sorted(by_semantic.values(), key=lambda row: row.decision_key)
    best = ordered[0]
    if len(ordered) > 1:
        second = ordered[1]
        if second.strength_key == best.strength_key and (second.locus_id, second.side) != (best.locus_id, best.side):
            return None, "PAIR_EQUAL_STRENGTH_LOCUS_CONFLICT"
    return best, "PAIR_UNIQUE_WINNER"


@dataclass(frozen=True, slots=True)
class MateRescueLocus:
    """Row of the mate-anchored locus record (coordinates quantized to one decimal)."""

    locus_id: str
    chrom: str
    family: str
    candidate_start: float
    candidate_end: float
    candidate_midpoint: float
    left_pair_ids: tuple[int, ...]
    right_pair_ids: tuple[int, ...]
    passes: bool


@dataclass(frozen=True, slots=True)
class MateRescueResult:
    guarded_locus_ids: tuple[str, ...]
    ledger: tuple[MateRescueLocus, ...]
    calls: tuple[IncumbentCall, ...]
    counts: Mapping[str, int]


def locus_redundant_with_calls(locus: CandidateLocus, calls: Sequence[IncumbentCall], guard_bp: int) -> bool:
    family = normalize_family(locus.family)
    for call in calls:
        if call.chrom != locus.chrom or normalize_family(call.family) != family:
            continue
        if interval_gap(locus.start, locus.end, call.start, call.end) <= guard_bp:
            return True
    return False


def guard_candidate_loci(
    loci: Sequence[CandidateLocus], split_junction_calls: Sequence[IncumbentCall]
) -> tuple[list[CandidateLocus], list[str]]:
    """Return the loci not redundant with a split-junction call, and the guarded locus IDs."""

    eligible: list[CandidateLocus] = []
    guarded: list[str] = []
    for locus in loci:
        if locus_redundant_with_calls(locus, split_junction_calls, CORE_GUARD_BP):
            guarded.append(locus.locus_id)
        else:
            eligible.append(locus)
    return eligible, guarded


def mate_rescue_evidence(
    *,
    loci: Sequence[CandidateLocus],
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    te_index: LocalTEIndex,
    progress: Callable[[int, int], None] | None = None,
) -> tuple[list[MateEvidence], dict[str, int]]:
    """Return confirmed mate evidence of eligible ``loci`` in locus order (batchable)."""

    counts: Counter[str] = Counter()
    evidence: list[MateEvidence] = []
    for index, locus in enumerate(loci, 1):
        contig = reference[locus.chrom]
        center = int(round(locus.midpoint))
        scope_start = max(0, center - SCOPE_RADIUS)
        scope = contig[scope_start : min(len(contig), center + SCOPE_RADIUS + 1)]
        kmers = family_kmers.get(locus.family, set())
        if not kmers:
            counts["NO_CANDIDATE_FAMILY_KMER_SET"] += 1
            continue
        for pair_id in locus.pair_ids:
            if pair_id not in reads:
                counts["READ_MISSING"] += 1
                continue
            mate1, mate2 = reads[pair_id]
            for te_mate, te_read, host_mate, host_read in ((1, mate1, 2, mate2), (2, mate2, 1, mate1)):
                host, reason = best_local_host_anchor(host_read, scope, scope_start, locus.start, locus.end)
                if host is None:
                    counts[reason] += 1
                    continue
                side = classify_side(host.view, host.midpoint, locus.midpoint)
                if side is None:
                    counts["PAIR_GEOMETRY_INCONSISTENT"] += 1
                    continue
                window, reason = best_te_mate_window(te_read, locus.family, kmers, te_index, scope)
                if window is None:
                    counts[reason] += 1
                    continue
                evidence.append(MateEvidence(
                    locus.locus_id, locus.chrom, locus.family, pair_id, side, te_mate, host_mate,
                    host.view, host.start, host.end, host.midpoint, host.edits, host.candidate_distance,
                    window.seed_hits, window.seed_fraction, window.te_edits,
                ))
                counts["R4_PAIR_EVIDENCE_CERTIFIED"] += 1
        if progress is not None:
            progress(index, len(loci))
    return evidence, dict(counts)


def finalize_mate_rescue(
    sample_id: str,
    eligible: Sequence[CandidateLocus],
    guarded: Sequence[str],
    evidence: Sequence[MateEvidence],
    counts: Mapping[str, int] | None = None,
) -> MateRescueResult:
    by_pair: dict[int, list[MateEvidence]] = defaultdict(list)
    for row in evidence:
        by_pair[row.pair_id].append(row)
    winners: dict[str, list[MateEvidence]] = defaultdict(list)
    for pair_id in sorted(by_pair):
        winner, _ = resolve_pair_evidence(by_pair[pair_id])
        if winner is not None:
            winners[winner.locus_id].append(winner)
    ledger: list[MateRescueLocus] = []
    calls: list[IncumbentCall] = []
    for locus in sorted(eligible, key=lambda row: row.locus_id):
        rows = winners.get(locus.locus_id, [])
        left = tuple(sorted({row.pair_id for row in rows if row.side == "LEFT"}))
        right = tuple(sorted({row.pair_id for row in rows if row.side == "RIGHT"}))
        passes = len(left) >= MIN_SIDE_FRAGMENTS and len(right) >= MIN_SIDE_FRAGMENTS
        ledger.append(MateRescueLocus(
            locus.locus_id, locus.chrom, locus.family, quantize(locus.start), quantize(locus.end),
            quantize(locus.midpoint), left, right, passes,
        ))
        if not passes:
            continue
        start = int(math.floor(min(locus.start, locus.end)))
        end = int(math.ceil(max(locus.start, locus.end))) + 1
        if end <= start:
            end = start + 1
        union = tuple(sorted(set(left) | set(right)))
        calls.append(IncumbentCall(
            sample_id, f"V2R4_MATE_{len(calls) + 1:06d}", IncumbentCallClass.MATE_RESCUE,
            locus.chrom, start, end, locus.family, None, None, None, len(left), len(right),
            "two_sided_mate_cluster_rescue", float(len(union)), union,
        ))
    return MateRescueResult(tuple(guarded), tuple(ledger), tuple(calls), dict(counts or {}))


def run_mate_rescue(
    *,
    sample_id: str,
    loci: Sequence[CandidateLocus],
    split_junction_calls: Sequence[IncumbentCall],
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    te_index: LocalTEIndex,
    progress: Callable[[int, int], None] | None = None,
) -> MateRescueResult:
    eligible, guarded = guard_candidate_loci(loci, split_junction_calls)
    evidence, counts = mate_rescue_evidence(
        loci=eligible, reads=reads, reference=reference, family_kmers=family_kmers,
        te_index=te_index, progress=progress,
    )
    return finalize_mate_rescue(sample_id, eligible, guarded, evidence, counts)


# --------------------------------------------------------------------------
# Reference-discriminative strict pair evidence
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class OneSidedTarget:
    locus_id: str
    chrom: str
    family: str
    candidate_start: float
    candidate_end: float
    candidate_midpoint: float
    missing_side: str
    existing_left_ids: tuple[int, ...]
    existing_right_ids: tuple[int, ...]


def one_sided_targets(ledger: Sequence[MateRescueLocus]) -> tuple[OneSidedTarget, ...]:
    targets: list[OneSidedTarget] = []
    for row in ledger:
        if row.passes:
            continue
        left, right = len(row.left_pair_ids), len(row.right_pair_ids)
        if left >= MIN_SIDE_FRAGMENTS and right < MIN_SIDE_FRAGMENTS:
            missing = "RIGHT"
        elif right >= MIN_SIDE_FRAGMENTS and left < MIN_SIDE_FRAGMENTS:
            missing = "LEFT"
        else:
            continue
        targets.append(OneSidedTarget(
            row.locus_id, row.chrom, row.family, row.candidate_start, row.candidate_end,
            row.candidate_midpoint, missing, row.left_pair_ids, row.right_pair_ids,
        ))
    targets.sort(key=lambda row: (row.chrom, row.family, row.candidate_start, row.locus_id))
    return tuple(targets)


def generate_host_seeds(contig: str, start: float, end: float, missing_side: str) -> list[tuple[int, str]]:
    middle = int(round((start + end) / 2.0))
    if missing_side == "LEFT":
        lo = max(0, int(math.floor(start)) - HOST_SEED_OUTER_PAD)
        hi = min(len(contig), middle)
    elif missing_side == "RIGHT":
        lo = max(0, middle)
        hi = min(len(contig), int(math.ceil(end)) + HOST_SEED_OUTER_PAD)
    else:
        raise ValueError(missing_side)
    if hi - lo < HOST_SEED_K:
        return []
    starts = list(range(lo, hi - HOST_SEED_K + 1, HOST_SEED_STRIDE))
    final_start = hi - HOST_SEED_K
    if final_start >= lo and final_start not in starts:
        starts.append(final_start)
    out: list[tuple[int, str]] = []
    seen: set[str] = set()
    for position in starts:
        seed = contig[position : position + HOST_SEED_K].upper()
        if "N" in seed or shannon_entropy_symbols(seed) < HOST_SEED_MIN_ENTROPY:
            continue
        canonical = min(seed, revcomp(seed))
        if canonical in seen:
            continue
        seen.add(canonical)
        out.append((position, seed))
    return out


class SeedPairLocator:
    """Protocol-like base: bounded exact locate of a seed and its RC to pair IDs."""

    def seed_pairs(self, seed: str, row_cap: int) -> tuple[int, bool, set[int]]:  # pragma: no cover
        raise NotImplementedError


@dataclass(frozen=True, slots=True)
class RoutingResult:
    targets: tuple[OneSidedTarget, ...]
    pair_to_loci: Mapping[int, tuple[str, ...]]
    seed_count: Mapping[tuple[int, str], int]
    saturated_loci: frozenset[str]


def route_one_sided_targets(
    *,
    targets: Sequence[OneSidedTarget],
    reserved_pair_ids: set[int],
    reference: Mapping[str, str],
    locator: SeedPairLocator,
) -> RoutingResult:
    pair_to_loci: dict[int, set[str]] = defaultdict(set)
    seed_count: Counter[tuple[int, str]] = Counter()
    saturated: set[str] = set()
    for target in targets:
        seeds = generate_host_seeds(reference[target.chrom], target.candidate_start, target.candidate_end, target.missing_side)
        existing = set(target.existing_left_ids) | set(target.existing_right_ids)
        recruited: set[int] = set()
        for _, seed in seeds:
            _, complete, pair_ids = locator.seed_pairs(seed, HOST_SEED_ROW_CAP)
            if not complete:
                continue
            for pair_id in pair_ids - existing - reserved_pair_ids:
                recruited.add(pair_id)
                seed_count[(pair_id, target.locus_id)] += 1
            if len(recruited) > MAX_RECRUITED_PAIRS_PER_LOCUS:
                saturated.add(target.locus_id)
                break
        if target.locus_id in saturated:
            recruited = set()
        for pair_id in recruited:
            pair_to_loci[pair_id].add(target.locus_id)
    return RoutingResult(
        tuple(targets),
        {pair_id: tuple(sorted(values)) for pair_id, values in sorted(pair_to_loci.items())},
        dict(seed_count),
        frozenset(saturated),
    )


def ordinary_reference_pair_possible(
    te_read: str, scope: str, scope_start: int, host_start: int, host_end: int, host_view: str,
    fragment_q95: float = FRAGMENT_Q95,
) -> bool:
    for view, sequence in (("FWD", te_read.upper()), ("RC", revcomp(te_read.upper()))):
        distance, locations = hw_locations(sequence, scope, TE_FULL_READ_REF_MAX_EDITS)
        if distance is None:
            continue
        for start, end in locations:
            target_start = scope_start + start
            target_end = scope_start + end + 1
            # The 95th percentile of the fragment size: a pair whose span is above it is unlikely to be a
            # normal pair, so it cannot explain the read pair as an ordinary reference pair.
            if max(host_end, target_end) - min(host_start, target_start) > fragment_q95:
                continue
            if host_start <= target_start:
                inward = host_view == "FWD" and view == "RC"
            else:
                inward = view == "FWD" and host_view == "RC"
            if inward:
                return True
    return False


def placement_summary(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    valid = [row for row in rows if int(row["nm"]) <= HOST_MATE_MAX_EDITS and not row["supplementary"]]
    saturated = len(rows) >= GLOBAL_PLACEMENT_CAP
    if not valid:
        return {"mapped": False, "best_nm": None, "best": [], "best_count": 0, "saturated": saturated}
    best_nm = min(int(row["nm"]) for row in valid)
    dedup: dict[tuple, dict[str, object]] = {}
    for row in valid:
        if int(row["nm"]) == best_nm:
            dedup[(row["rname"], row["pos0"], row["rev"], row["cigar"])] = row
    return {"mapped": True, "best_nm": best_nm, "best": list(dedup.values()), "best_count": len(dedup), "saturated": saturated}


@dataclass(frozen=True, slots=True)
class StrictPairWinner:
    locus_id: str
    side: str
    pair_id: int


@dataclass(frozen=True, slots=True)
class StrictPairResult:
    winners: tuple[StrictPairWinner, ...]
    counts: Mapping[str, int]
    evidence: tuple[MateEvidence, ...] = ()


def certify_strict_pairs(
    *,
    routing: RoutingResult,
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    te_index: LocalTEIndex,
    bowtie2: str | Path,
    host_index_prefix: str | Path,
    work_directory: Path,
    fragment_q95: float = FRAGMENT_Q95,
) -> StrictPairResult:
    targets = {target.locus_id: target for target in routing.targets}
    counts: Counter[str] = Counter()
    prefiltered: list[tuple[int, str, int, int]] = []
    for pair_id, locus_ids in routing.pair_to_loci.items():
        mate1, mate2 = reads[pair_id]
        for locus_id in locus_ids:
            if locus_id in routing.saturated_loci:
                continue
            target = targets[locus_id]
            midpoint = (target.candidate_start + target.candidate_end) / 2.0
            contig = reference[target.chrom]
            center = int(round(midpoint))
            scope_start = max(0, center - SCOPE_RADIUS)
            scope = contig[scope_start : min(len(contig), center + SCOPE_RADIUS + 1)]
            for te_mate, host_mate, host_read in ((1, 2, mate2), (2, 1, mate1)):
                host, reason = best_local_host_anchor(
                    host_read, scope, scope_start, target.candidate_start, target.candidate_end
                )
                if host is None:
                    counts[reason] += 1
                    continue
                side = classify_side(host.view, host.midpoint, midpoint)
                if side is None:
                    counts["PAIR_GEOMETRY_INCONSISTENT_LOCAL"] += 1
                    continue
                if side != target.missing_side:
                    counts["LOCAL_HOST_NOT_MISSING_SIDE"] += 1
                    continue
                prefiltered.append((pair_id, locus_id, te_mate, host_mate))
    if not prefiltered:
        return StrictPairResult((), dict(counts), ())
    queries: dict[str, str] = {}
    for pair_id, _, _, host_mate in prefiltered:
        queries[f"P{pair_id}_M{host_mate}"] = reads[pair_id][host_mate - 1]
    work_directory.mkdir(parents=True, exist_ok=True)
    fasta = work_directory / "strict_pair_host_queries.fa"
    sam = work_directory / "strict_pair_host_queries.sam"
    write_fasta(fasta, sorted(queries.items()))
    run_bowtie2(bowtie2, [*HOST_GLOBAL_OPTIONS, "-k", str(GLOBAL_PLACEMENT_CAP), "-x", str(host_index_prefix),
                          "-U", str(fasta), "-S", str(sam)])
    alignments = parse_sam_records(sam)
    evidence: dict[int, list[MateEvidence]] = defaultdict(list)
    for pair_id, locus_id, te_mate, host_mate in prefiltered:
        target = targets[locus_id]
        midpoint = (target.candidate_start + target.candidate_end) / 2.0
        host_read = reads[pair_id][host_mate - 1]
        te_read = reads[pair_id][te_mate - 1]
        summary = placement_summary(alignments.get(f"P{pair_id}_M{host_mate}", []))
        if summary["saturated"]:
            counts["UNRESOLVED_HOST_GLOBAL_ENUMERATION_SATURATED"] += 1
            continue
        if not summary["mapped"]:
            counts["REJECT_HOST_NO_GLOBAL_MATCH"] += 1
            continue
        if summary["best_count"] != 1:
            counts["REJECT_HOST_GLOBAL_BEST_MULTIMAP"] += 1
            continue
        placement = summary["best"][0]
        host_start = int(placement["pos0"])
        host_end = host_start + cigar_ref_span(str(placement["cigar"]), len(host_read))
        host_midpoint = (host_start + host_end) / 2.0
        host_view = "RC" if placement["rev"] else "FWD"
        if not (placement["rname"] == target.chrom and interval_distance_point(
            host_midpoint, target.candidate_start, target.candidate_end) <= CANDIDATE_SLOP):
            counts["REJECT_HOST_GLOBAL_BEST_OFF_NOMINATED_LOCUS"] += 1
            continue
        side = classify_side(host_view, host_midpoint, midpoint)
        if side is None:
            counts["REJECT_GLOBAL_HOST_GEOMETRY_INCONSISTENT"] += 1
            continue
        if side != target.missing_side:
            counts["REJECT_GLOBAL_HOST_NOT_MISSING_SIDE"] += 1
            continue
        contig = reference[target.chrom]
        center = int(round(midpoint))
        scope_start = max(0, center - SCOPE_RADIUS)
        scope = contig[scope_start : min(len(contig), center + SCOPE_RADIUS + 1)]
        if ordinary_reference_pair_possible(te_read, scope, scope_start, host_start, host_end, host_view, fragment_q95):
            counts["REJECT_ORDINARY_CONCORDANT_REFERENCE_PAIR"] += 1
            continue
        kmers = family_kmers.get(target.family, set())
        if not kmers:
            counts["REJECT_NO_FAMILY_KMER_PANEL"] += 1
            continue
        window, reason = best_te_window_symmetric(te_read, target.family, kmers, te_index, scope)
        if window is None:
            counts[f"REJECT_{reason}"] += 1
            continue
        evidence[pair_id].append(MateEvidence(
            locus_id, target.chrom, target.family, pair_id, side, te_mate, host_mate, host_view,
            host_start, host_end, host_midpoint, int(summary["best_nm"]),
            interval_distance_point(host_midpoint, target.candidate_start, target.candidate_end),
            window.seed_hits, window.seed_fraction, window.te_edits,
            (-routing.seed_count.get((pair_id, locus_id), 0),),
        ))
        counts["R6_REFERENCE_DISCRIMINATIVE_EVIDENCE_CERTIFIED"] += 1
    winners: list[StrictPairWinner] = []
    for pair_id in sorted(evidence):
        winner, _ = resolve_pair_evidence(evidence[pair_id])
        if winner is not None:
            winners.append(StrictPairWinner(winner.locus_id, winner.side, winner.pair_id))
    return StrictPairResult(
        tuple(winners), dict(counts), tuple(row for pair_id in sorted(evidence) for row in evidence[pair_id])
    )
