"""Fixed-candidate split-junction reference competition (incumbent class 1).

Port of the frozen FCRC1 evidence layer.  For every candidate locus and both
mates of every candidate physical pair, terminal host|TE and TE|host cuts are
nominated from candidate-family 21-mer density, then certified only if:

1. the whole host-facing segment aligns to the local reference scope with at
   most three edits and yields exactly one best breakpoint, inside the
   candidate envelope plus 250 bp;
2. the TE-facing segment is not a reference continuation at that breakpoint
   (<= 3 edits, 3-bp pad);
3. the TE-facing segment has a unique best biological family in the full TE
   library (<= 3 edits);
4. the split model beats a whole-read local reference alignment by >= 4 edits.

One semantic winner per read context, fragment-level reconciliation of the two
mates, side clustering of physical-pair support, pairing of opposite sides
within 50 bp, and the frozen score >= 4 emission rule complete the layer.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

from .alignment import hw_distance, hw_locations
from .candidates import CandidateLocus
from .events import IncumbentCall, IncumbentCallClass
from .sequence import count_occurrences, distance_to_interval, interval_hits, kmer_hit_prefix, revcomp
from .te_library import LocalTEIndex, QualifiedObservation, SideCluster, cluster_observations, family_best_classify

SCOPE_RADIUS = 800
TE_K = 21
MIN_TE_KMER_ENTROPY = 1.10
MIN_TE_BP = 35
MIN_HOST_BP = 24
MIN_TE_SEED_HITS = 2
MIN_TE_SEED_FRACTION = 0.35
HOST_SEED_BP = 16
MAX_NOMINATIONS_PER_READ = 12
HOST_MAX_EDITS = 3
CANDIDATE_SLOP = 250
CONTINUATION_MAX_EDITS = 3
CONTINUATION_PAD = 3
REF_FULL_MAX_EDITS = 12
MIN_REF_ADVANTAGE = 4
PAIR_WINDOW_BP = 50
CENTRAL_SCORE_THRESHOLD = 4.0
HOST_SEED_OCCURRENCE_CAP = 100

__all__ = [
    "Nomination",
    "CertifiedJunction",
    "nominate_read",
    "certify_nomination",
    "resolve_context",
    "resolve_fragment",
    "pair_side_clusters",
    "SplitJunctionResult",
    "split_junction_hypotheses",
    "finalize_split_junctions",
    "run_split_junction_layer",
]


@dataclass(frozen=True, slots=True)
class Nomination:
    locus_id: str
    chrom: str
    candidate_family: str
    pair_id: int
    mate: int
    view: str
    side_kind: str
    split: int
    te_bp: int
    host_bp: int
    te_seed_hits: int
    te_seed_total: int
    te_seed_fraction: float
    host_seed_matches: int
    host_seed_first_pos: int | None
    nomination_score: float


@dataclass(frozen=True, slots=True)
class CertifiedJunction:
    locus_id: str
    chrom: str
    pair_id: int
    mate: int
    view: str
    event_side: str
    breakpoint: int
    candidate_family: str
    verified_family: str
    fbte_id: str
    te_orientation: str
    te_bp: int
    host_bp: int
    host_edits: int
    te_edits: int
    split_cost: int
    ref_full_edits: int | None
    ref_advantage: int
    continuation_edits: int | None
    candidate_distance: float
    nomination_score: float

    @property
    def context_key(self) -> tuple[int, int]:
        return (self.pair_id, self.mate)

    @property
    def semantic_key(self) -> tuple[str, str, int]:
        return (self.verified_family, self.event_side, self.breakpoint)

    @property
    def decision_key(self) -> tuple:
        return (
            self.split_cost,
            -self.ref_advantage,
            self.host_edits,
            self.te_edits,
            -self.host_bp,
            -self.te_bp,
            self.candidate_distance,
            self.verified_family,
            self.event_side,
            self.breakpoint,
            self.locus_id,
        )

    @property
    def strength_key(self) -> tuple:
        return (
            self.split_cost,
            self.ref_advantage,
            self.host_edits,
            self.te_edits,
            self.host_bp,
            self.te_bp,
            round(self.candidate_distance, 6),
        )


def nominate_read(
    *,
    locus: CandidateLocus,
    pair_id: int,
    mate: int,
    read: str,
    family_kmers: set[str],
    scope: str,
    scope_start: int,
) -> list[Nomination]:
    """Broad cut nomination; never evidence until certified."""

    if not family_kmers:
        return []
    raw: list[Nomination] = []
    for view, seq in (("FWD", read.upper()), ("RC", revcomp(read))):
        length = len(seq)
        prefix = kmer_hit_prefix(seq, family_kmers, TE_K)
        for split in range(MIN_HOST_BP, length - MIN_TE_BP + 1):
            hits, total = interval_hits(prefix, split, length, TE_K)
            fraction = hits / total if total else 0.0
            if hits >= MIN_TE_SEED_HITS and fraction >= MIN_TE_SEED_FRACTION and split >= HOST_SEED_BP:
                seed = seq[split - HOST_SEED_BP : split]
                if "N" not in seed:
                    matches, first = count_occurrences(scope, seed, stop_after=HOST_SEED_OCCURRENCE_CAP)
                    if matches:
                        score = fraction * 10.0 + min(hits, 10) + (5.0 if matches == 1 else 1.0 / matches)
                        raw.append(
                            Nomination(
                                locus.locus_id, locus.chrom, locus.family, pair_id, mate, view,
                                "HOST_LEFT_TE_RIGHT", split, length - split, split,
                                hits, total, fraction, matches,
                                None if first is None else scope_start + first, score,
                            )
                        )
            hits, total = interval_hits(prefix, 0, split, TE_K)
            fraction = hits / total if total else 0.0
            if hits >= MIN_TE_SEED_HITS and fraction >= MIN_TE_SEED_FRACTION and length - split >= HOST_SEED_BP:
                seed = seq[split : split + HOST_SEED_BP]
                if "N" not in seed:
                    matches, first = count_occurrences(scope, seed, stop_after=HOST_SEED_OCCURRENCE_CAP)
                    if matches:
                        score = fraction * 10.0 + min(hits, 10) + (5.0 if matches == 1 else 1.0 / matches)
                        raw.append(
                            Nomination(
                                locus.locus_id, locus.chrom, locus.family, pair_id, mate, view,
                                "TE_LEFT_HOST_RIGHT", split, split, length - split,
                                hits, total, fraction, matches,
                                None if first is None else scope_start + first, score,
                            )
                        )
    best: dict[tuple[str, str, int], Nomination] = {}
    for nomination in raw:
        key = (nomination.view, nomination.side_kind, nomination.split)
        old = best.get(key)
        if old is None or nomination.nomination_score > old.nomination_score:
            best[key] = nomination
    ordered = sorted(
        best.values(), key=lambda row: (-row.nomination_score, row.split, row.view, row.side_kind)
    )
    return ordered[:MAX_NOMINATIONS_PER_READ]


def _host_breakpoint(host: str, event_side: str, scope: str, scope_start: int) -> dict[str, object]:
    if not host or "N" in host:
        return {"ok": False, "reason": "HOST_HAS_N_OR_EMPTY"}
    distance, locations = hw_locations(host, scope, HOST_MAX_EDITS)
    if distance is None:
        return {"ok": False, "reason": "WHOLE_HOST_NO_MATCH_LE_MAX"}
    if not locations:
        return {"ok": False, "reason": "WHOLE_HOST_NO_LOCATION"}
    breakpoints = {
        scope_start + end + 1 if event_side == "L" else scope_start + start for start, end in locations
    }
    if len(breakpoints) != 1:
        return {
            "ok": False,
            "reason": "WHOLE_HOST_AMBIGUOUS_BEST_BREAKPOINT",
            "host_edits": distance,
            "placement_count": len(breakpoints),
        }
    return {"ok": True, "breakpoint": next(iter(breakpoints)), "host_edits": distance, "placement_count": 1}


def _continuation_distance(contig: str, breakpoint: int, event_side: str, clip: str) -> int | None:
    n = len(clip)
    if event_side == "L":
        start = max(0, breakpoint - CONTINUATION_PAD)
        end = min(len(contig), breakpoint + n + CONTINUATION_PAD)
    else:
        start = max(0, breakpoint - n - CONTINUATION_PAD)
        end = min(len(contig), breakpoint + CONTINUATION_PAD)
    return hw_distance(clip, contig[start:end], CONTINUATION_MAX_EDITS)


def certify_nomination(
    *,
    nomination: Nomination,
    locus: CandidateLocus,
    read: str,
    contig: str,
    scope: str,
    scope_start: int,
    te_index: LocalTEIndex,
    reference_cache: dict[tuple[str, str], int | None],
) -> tuple[CertifiedJunction | None, str, dict[str, object]]:
    oriented = read.upper() if nomination.view == "FWD" else revcomp(read)
    if nomination.side_kind == "HOST_LEFT_TE_RIGHT":
        host, clip, event_side = oriented[: nomination.split], oriented[nomination.split :], "L"
    elif nomination.side_kind == "TE_LEFT_HOST_RIGHT":
        clip, host, event_side = oriented[: nomination.split], oriented[nomination.split :], "R"
    else:
        return None, "UNKNOWN_SIDE_KIND", {}
    if not clip or "N" in clip:
        return None, "TE_CLIP_HAS_N_OR_EMPTY", {}
    placement = _host_breakpoint(host, event_side, scope, scope_start)
    if not placement.get("ok"):
        return None, str(placement.get("reason", "HOST_REJECT")), placement
    breakpoint = int(placement["breakpoint"])
    candidate_distance = distance_to_interval(breakpoint, locus.start, locus.end)
    if candidate_distance > CANDIDATE_SLOP:
        return None, "BREAKPOINT_OUTSIDE_CANDIDATE_ENVELOPE", {
            "breakpoint": breakpoint,
            "candidate_distance": candidate_distance,
            "host_edits": placement.get("host_edits"),
        }
    continuation = _continuation_distance(contig, breakpoint, event_side, clip)
    if continuation is not None:
        return None, "REFERENCE_CONTINUATION_COMPATIBLE", {
            "breakpoint": breakpoint,
            "host_edits": placement.get("host_edits"),
            "continuation_edits": continuation,
            "candidate_distance": candidate_distance,
        }
    hit, diagnostic = family_best_classify(te_index.record_hits(clip))
    if hit is None:
        return None, str(diagnostic.get("reason", "TE_VALIDATION_REJECT")), {
            "breakpoint": breakpoint,
            "host_edits": placement.get("host_edits"),
            "best_distance": diagnostic.get("best_distance", ""),
            "best_families": diagnostic.get("best_families", ""),
            "candidate_distance": candidate_distance,
        }
    split_cost = int(placement["host_edits"]) + int(hit.edit_distance)
    cache_key = (nomination.view, locus.locus_id)
    if cache_key not in reference_cache:
        reference_cache[cache_key] = hw_distance(oriented, scope, REF_FULL_MAX_EDITS)
    reference_distance = reference_cache[cache_key]
    advantage = (
        REF_FULL_MAX_EDITS + 1 - split_cost if reference_distance is None else reference_distance - split_cost
    )
    if advantage < MIN_REF_ADVANTAGE:
        return None, "REFERENCE_ONLY_COMPETITIVE", {
            "breakpoint": breakpoint,
            "host_edits": placement.get("host_edits"),
            "te_edits": hit.edit_distance,
            "split_cost": split_cost,
            "ref_full_edits": "" if reference_distance is None else reference_distance,
            "ref_advantage": advantage,
            "candidate_distance": candidate_distance,
            "verified_family": hit.family,
        }
    junction = CertifiedJunction(
        locus_id=locus.locus_id,
        chrom=locus.chrom,
        pair_id=nomination.pair_id,
        mate=nomination.mate,
        view=nomination.view,
        event_side=event_side,
        breakpoint=breakpoint,
        candidate_family=nomination.candidate_family,
        verified_family=hit.family,
        fbte_id=hit.record_id,
        te_orientation=hit.orientation,
        te_bp=len(clip),
        host_bp=len(host),
        host_edits=int(placement["host_edits"]),
        te_edits=int(hit.edit_distance),
        split_cost=split_cost,
        ref_full_edits=reference_distance,
        ref_advantage=advantage,
        continuation_edits=continuation,
        candidate_distance=float(candidate_distance),
        nomination_score=float(nomination.nomination_score),
    )
    return junction, "CERTIFIED", {}


def resolve_context(rows: Sequence[CertifiedJunction]) -> tuple[CertifiedJunction | None, str]:
    if not rows:
        return None, "NO_CERTIFIED_HYPOTHESIS"
    by_semantic: dict[tuple[str, str, int], CertifiedJunction] = {}
    for row in rows:
        old = by_semantic.get(row.semantic_key)
        if old is None or row.decision_key < old.decision_key:
            by_semantic[row.semantic_key] = row
    ordered = sorted(by_semantic.values(), key=lambda row: row.decision_key)
    best = ordered[0]
    if len(ordered) > 1 and ordered[1].strength_key == best.strength_key and ordered[1].semantic_key != best.semantic_key:
        return None, "READ_CONTEXT_EQUAL_STRENGTH_TIE"
    return best, "UNIQUE_CONTEXT_WINNER"


def resolve_fragment(rows: Sequence[CertifiedJunction]) -> tuple[list[CertifiedJunction], str]:
    if not rows:
        return [], "NO_CONTEXT_WINNER"
    ordered = sorted(rows, key=lambda row: row.decision_key)
    if len(ordered) == 1:
        return [ordered[0]], "ONE_MATE_WINNER"
    first, second = ordered[0], ordered[1]
    compatible = first.chrom == second.chrom and first.verified_family == second.verified_family
    gap = abs(first.breakpoint - second.breakpoint)
    if compatible and first.event_side != second.event_side and gap <= 50:
        return [first, second], "COMPATIBLE_TWO_SIDED_FRAGMENT"
    if compatible and first.event_side == second.event_side and gap <= 8:
        return [first, second], "COMPATIBLE_SAME_SIDE_FRAGMENT"
    if first.strength_key == second.strength_key:
        return [], "FRAGMENT_INCOMPATIBLE_EQUAL_STRENGTH"
    return [first], "FRAGMENT_INCOMPATIBLE_KEEP_STRONGER"


def _pair_ids(reads: set[str]) -> tuple[int, ...]:
    return tuple(sorted(int(value) for value in reads))


def pair_side_clusters(clusters: Sequence[SideCluster], sample_id: str) -> list[IncumbentCall]:
    """Pair opposite sides within 50 bp (nearest first) and number every event."""

    by: dict[tuple[str, str], dict[str, list[SideCluster]]] = defaultdict(lambda: {"L": [], "R": []})
    for cluster in clusters:
        by[(cluster.chrom, cluster.family)][cluster.side].append(cluster)
    rows: list[dict[str, object]] = []
    for (chrom, family), sides in sorted(by.items()):
        lefts = sorted(sides["L"], key=lambda item: item.coord)
        rights = sorted(sides["R"], key=lambda item: item.coord)
        candidates = []
        for left_index, left in enumerate(lefts):
            for right_index, right in enumerate(rights):
                gap = abs(right.coord - left.coord)
                if gap <= PAIR_WINDOW_BP:
                    candidates.append((gap, -len(left.reads | right.reads), left_index, right_index))
        used_left: set[int] = set()
        used_right: set[int] = set()
        for _, _, left_index, right_index in sorted(candidates):
            if left_index in used_left or right_index in used_right:
                continue
            used_left.add(left_index)
            used_right.add(right_index)
            left, right = lefts[left_index], rights[right_index]
            fragments = left.reads | right.reads
            rows.append({
                "chrom": chrom, "start": min(left.coord, right.coord), "end": max(left.coord, right.coord),
                "family": family, "left_coord": left.coord, "right_coord": right.coord,
                "delta": right.coord - left.coord, "left_support": len(left.reads),
                "right_support": len(right.reads), "evidence_type": "two_sided",
                "score": float(len(fragments) + 2), "pairs": _pair_ids(fragments),
            })
        for index, cluster in enumerate(lefts):
            if index not in used_left:
                rows.append({
                    "chrom": chrom, "start": cluster.coord, "end": cluster.coord, "family": family,
                    "left_coord": cluster.coord, "right_coord": -1, "delta": 0,
                    "left_support": len(cluster.reads), "right_support": 0,
                    "evidence_type": "one_sided_left", "score": float(len(cluster.reads)),
                    "pairs": _pair_ids(cluster.reads),
                })
        for index, cluster in enumerate(rights):
            if index not in used_right:
                rows.append({
                    "chrom": chrom, "start": cluster.coord, "end": cluster.coord, "family": family,
                    "left_coord": -1, "right_coord": cluster.coord, "delta": 0,
                    "left_support": 0, "right_support": len(cluster.reads),
                    "evidence_type": "one_sided_right", "score": float(len(cluster.reads)),
                    "pairs": _pair_ids(cluster.reads),
                })
    rows.sort(key=lambda row: (row["chrom"], row["start"], row["end"], row["family"], row["evidence_type"]))
    return [
        IncumbentCall(
            sample_id=sample_id,
            legacy_event_id=f"V2FCRC1_EVENT_{index:06d}",
            call_class=IncumbentCallClass.SPLIT_JUNCTION,
            chrom=str(row["chrom"]),
            start=int(row["start"]),
            end=int(row["end"]),
            family=str(row["family"]),
            left_coord=int(row["left_coord"]),
            right_coord=int(row["right_coord"]),
            delta=int(row["delta"]),
            left_support=int(row["left_support"]),
            right_support=int(row["right_support"]),
            evidence_type=str(row["evidence_type"]),
            score=float(row["score"]),
            fragment_pair_ids=tuple(row["pairs"]),
        )
        for index, row in enumerate(rows, 1)
    ]


@dataclass(frozen=True, slots=True)
class SplitJunctionResult:
    all_events: tuple[IncumbentCall, ...]
    calls: tuple[IncumbentCall, ...]
    certified_junctions: tuple[CertifiedJunction, ...]
    counts: Mapping[str, int]


def split_junction_hypotheses(
    *,
    loci: Sequence[CandidateLocus],
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    te_index: LocalTEIndex,
    progress: Callable[[int, int], None] | None = None,
    hypothesis_sink: Callable[[Nomination, str, CertifiedJunction | None, dict[str, object]], None] | None = None,
) -> tuple[list[CertifiedJunction], dict[str, int]]:
    """Certified split-junction hypotheses of ``loci`` in locus order.

    Loci are independent, so contiguous locus batches may be processed
    separately and concatenated in order before :func:`finalize_split_junctions`.
    """

    certified: list[CertifiedJunction] = []
    counts: dict[str, int] = defaultdict(int)
    for locus_index, locus in enumerate(loci, 1):
        contig = reference[locus.chrom]
        center = int(round(locus.midpoint))
        scope_start = max(0, center - SCOPE_RADIUS)
        scope_end = min(len(contig), center + SCOPE_RADIUS + 1)
        scope = contig[scope_start:scope_end]
        kmers = family_kmers.get(locus.family, set())
        if not kmers:
            counts["NO_CANDIDATE_FAMILY_KMER_SET"] += 1
            continue
        for pair_id in locus.pair_ids:
            mate1, mate2 = reads[pair_id]
            for mate, read in ((1, mate1), (2, mate2)):
                cache: dict[tuple[str, str], int | None] = {}
                for nomination in nominate_read(
                    locus=locus, pair_id=pair_id, mate=mate, read=read,
                    family_kmers=kmers, scope=scope, scope_start=scope_start,
                ):
                    junction, reason, diagnostic = certify_nomination(
                        nomination=nomination, locus=locus, read=read, contig=contig,
                        scope=scope, scope_start=scope_start, te_index=te_index, reference_cache=cache,
                    )
                    counts[reason] += 1
                    if hypothesis_sink is not None:
                        hypothesis_sink(nomination, reason, junction, diagnostic)
                    if junction is not None:
                        certified.append(junction)
        if progress is not None:
            progress(locus_index, len(loci))
    return certified, dict(counts)


def finalize_split_junctions(
    sample_id: str, certified: Sequence[CertifiedJunction], counts: Mapping[str, int] | None = None
) -> SplitJunctionResult:
    """Context and fragment competition, side clustering and pairing."""

    context_hypotheses: dict[tuple[int, int], list[CertifiedJunction]] = defaultdict(list)
    for junction in certified:
        context_hypotheses[junction.context_key].append(junction)
    context_winners: dict[int, list[CertifiedJunction]] = defaultdict(list)
    for rows in context_hypotheses.values():
        winner, _ = resolve_context(rows)
        if winner is not None:
            context_winners[winner.pair_id].append(winner)
    final: list[CertifiedJunction] = []
    for pair_id in sorted(context_winners):
        kept, _ = resolve_fragment(context_winners[pair_id])
        final.extend(kept)
    observations = [
        QualifiedObservation(
            read_id=str(row.pair_id), chrom=row.chrom, side=row.event_side, breakpoint=row.breakpoint,
            clip_len=row.te_bp, family=row.verified_family, fbte_id=row.fbte_id,
            edit_distance=row.te_edits, orientation=row.te_orientation,
        )
        for row in final
    ]
    events = pair_side_clusters(cluster_observations(observations), sample_id)
    calls = tuple(row for row in events if row.score >= CENTRAL_SCORE_THRESHOLD)
    return SplitJunctionResult(tuple(events), calls, tuple(final), dict(counts or {}))


def run_split_junction_layer(
    *,
    sample_id: str,
    loci: Sequence[CandidateLocus],
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    te_index: LocalTEIndex,
    progress: Callable[[int, int], None] | None = None,
    hypothesis_sink: Callable[[Nomination, str, CertifiedJunction | None, dict[str, object]], None] | None = None,
) -> SplitJunctionResult:
    """Run the frozen split-junction layer over every candidate locus."""

    certified, counts = split_junction_hypotheses(
        loci=loci, reads=reads, reference=reference, family_kmers=family_kmers, te_index=te_index,
        progress=progress, hypothesis_sink=hypothesis_sink,
    )
    return finalize_split_junctions(sample_id, certified, counts)
