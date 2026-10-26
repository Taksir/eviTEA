"""Insertion-specific junction groups on one-sided loci.

For every one-sided target locus, reads of msBWT-routed physical pairs are
scanned for contiguous host|TE or TE|host boundaries.  An accepted adjacency
needs: an exact globally unique host half near the nominated locus on the
missing genomic side, no contiguous reference explanation of the 40-bp
signature within two edits, and a complete bounded exact msBWT locate of the
signature.  Located pairs are collapsed into exact pair-sequence groups
(``sha256(mate1 + "|" + mate2)``), groups touching pairs used by earlier calls
or strict winners are excluded, one realization per (locus, side) competes
globally, and each surviving group adds at most one representative pair to one
locus side.

Output: the locus aggregation record and per-group junction source
for full TE-segment confirmation.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from .candidates import CandidateLocus
from .external import EXACTISH_OPTIONS, cigar_ref_span, parse_sam_records, run_bowtie2, write_fasta
from .mate_rescue import OneSidedTarget, RoutingResult, SeedPairLocator, StrictPairWinner
from .sequence import count_occurrences, interval_distance_point, interval_hits, kmer_hit_prefix, revcomp

SCOPE_RADIUS = 800
CANDIDATE_SLOP = 250
TE_K = 21
MIN_TE_BP = 35
MIN_HOST_BP = 24
MIN_TE_SEED_HITS = 2
MIN_TE_SEED_FRACTION = 0.35
JUNCTION_EACH_SIDE = 20
MAX_SPLITS_PER_READ = 4
MAX_DISCOVERY_SIGNATURES_PER_LOCUS = 16
REFERENCE_SIGNATURE_MAX_EDITS = 2
SIGNATURE_LOCATE_CAP = 256
HOST_PLACEMENT_REPORT = 20
MIN_SIGNATURE_HALF = 12
HOST_SEED_OCCURRENCE_CAP = 50


@dataclass(frozen=True, slots=True)
class ReadJunction:
    view: str
    side: str
    split: int
    te_bp: int
    host_bp: int
    te_seed_hits: int
    te_seed_total: int
    te_seed_fraction: float
    host_seed_matches_in_scope: int
    host_seed_first_pos: int | None
    signature: str


def infer_read_junctions(
    read: str, family_kmers: set[str], scope: str, scope_start: int
) -> list[ReadJunction]:
    """Return direct boundary candidates in one read (both views), strongest first."""

    if not family_kmers:
        return []
    candidates: list[tuple[float, ReadJunction]] = []
    for view, sequence in (("FWD", read.upper()), ("RC", revcomp(read))):
        length = len(sequence)
        prefix = kmer_hit_prefix(sequence, family_kmers, TE_K)
        for split in range(MIN_HOST_BP, length - MIN_TE_BP + 1):
            for side in ("HOST_LEFT_TE_RIGHT", "TE_LEFT_HOST_RIGHT"):
                if side == "HOST_LEFT_TE_RIGHT":
                    hits, total = interval_hits(prefix, split, length, TE_K)
                    host_ok = split >= JUNCTION_EACH_SIDE
                    host_seed = sequence[split - JUNCTION_EACH_SIDE : split]
                    te_bp, host_bp = length - split, split
                else:
                    hits, total = interval_hits(prefix, 0, split, TE_K)
                    host_ok = length - split >= JUNCTION_EACH_SIDE
                    host_seed = sequence[split : split + JUNCTION_EACH_SIDE]
                    te_bp, host_bp = split, length - split
                fraction = hits / total if total else 0.0
                if not (hits >= MIN_TE_SEED_HITS and fraction >= MIN_TE_SEED_FRACTION and host_ok):
                    continue
                if "N" in host_seed:
                    continue
                matches, first = count_occurrences(scope, host_seed, stop_after=HOST_SEED_OCCURRENCE_CAP)
                if not matches:
                    continue
                half = min(JUNCTION_EACH_SIDE, split, length - split)
                if half < MIN_SIGNATURE_HALF:
                    continue
                signature = sequence[split - half : split + half]
                junction = ReadJunction(
                    view, side, split, te_bp, host_bp, hits, total, fraction, matches,
                    None if first is None else scope_start + first, signature,
                )
                score = fraction * 10 + min(hits, 10) + (5 if matches == 1 else 1 / matches) + half / 100
                candidates.append((score, junction))
    best: dict[tuple, tuple[float, ReadJunction]] = {}
    for score, junction in candidates:
        key = (junction.view, junction.side, junction.signature, junction.host_seed_first_pos)
        if key not in best or score > best[key][0]:
            best[key] = (score, junction)
    ordered = sorted(best.values(), key=lambda item: (-item[0], item[1].split, item[1].signature))
    return [item[1] for item in ordered[:MAX_SPLITS_PER_READ]]


@dataclass(frozen=True, slots=True)
class JunctionProposal:
    qid: str
    locus_id: str
    chrom: str
    family: str
    missing_side: str
    candidate_start: float
    candidate_end: float
    candidate_midpoint: float
    source_pair_id: int
    source_mate: int
    read_view: str
    junction_semantic: str
    split: int
    te_bp: int
    te_seed_hits: int
    te_seed_fraction: float
    host_seed_matches_in_scope: int
    host_seed_first_pos: int
    signature: str
    canonical_signature: str
    host_half: str
    proposal_score: float

    @property
    def strength_key(self) -> tuple:
        return (
            -round(self.te_seed_fraction, 8),
            -self.te_seed_hits,
            self.host_seed_matches_in_scope,
            abs(self.host_seed_first_pos - self.candidate_midpoint),
            self.canonical_signature,
        )


@dataclass(frozen=True, slots=True)
class GroupHypothesis:
    group_hash: str
    representative_pair_id: int
    group_pair_ids: tuple[int, ...]
    locus_id: str
    side: str
    family: str
    canonical_signature: str
    source_pair_id: int
    proposal_score: float
    te_seed_hits: int
    te_seed_fraction: float
    raw_msbwt_rows: int
    host_candidate_distance: float

    @property
    def strength_key(self) -> tuple:
        return (
            -round(self.proposal_score, 8),
            -self.te_seed_hits,
            -round(self.te_seed_fraction, 8),
            self.raw_msbwt_rows,
            round(self.host_candidate_distance, 6),
        )

    @property
    def decision_key(self) -> tuple:
        return (*self.strength_key, self.family, self.locus_id, self.side,
                self.canonical_signature, self.representative_pair_id)


def resolve_group_hypotheses(rows: Sequence[GroupHypothesis]) -> tuple[GroupHypothesis | None, str]:
    if not rows:
        return None, "NO_JUNCTION_GROUP_EVIDENCE"
    by_semantic: dict[tuple[str, str], GroupHypothesis] = {}
    for row in rows:
        key = (row.locus_id, row.side)
        old = by_semantic.get(key)
        if old is None or row.decision_key < old.decision_key:
            by_semantic[key] = row
    ordered = sorted(by_semantic.values(), key=lambda row: row.decision_key)
    best = ordered[0]
    if len(ordered) > 1 and ordered[1].strength_key == best.strength_key and (
        ordered[1].locus_id, ordered[1].side) != (best.locus_id, best.side):
        return None, "JUNCTION_GROUP_EQUAL_STRENGTH_LOCUS_CONFLICT"
    return best, "JUNCTION_GROUP_UNIQUE_R7_WINNER"


def pair_group_hash(mate1: str, mate2: str) -> str:
    return hashlib.sha256((mate1 + "|" + mate2).encode("utf-8")).hexdigest()


def _unique_exact_host(rows: Sequence[dict[str, object]]) -> dict[str, object] | None:
    placements: dict[tuple, dict[str, object]] = {}
    for row in rows:
        if row["supplementary"] or int(row["nm"]) != 0:
            continue
        placements[(row["rname"], int(row["pos0"]), bool(row["rev"]), row["cigar"])] = row
    return next(iter(placements.values())) if len(placements) == 1 else None


@dataclass(frozen=True, slots=True)
class JunctionLocus:
    locus_id: str
    chrom: str
    family: str
    candidate_start: float
    candidate_end: float
    existing_left_ids: tuple[int, ...]
    existing_right_ids: tuple[int, ...]
    strict_left_ids: tuple[int, ...]
    strict_right_ids: tuple[int, ...]
    left_group_hashes: tuple[str, ...]
    right_group_hashes: tuple[str, ...]
    left_representatives: tuple[int, ...]
    right_representatives: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class JunctionGroupsResult:
    ledger: tuple[JunctionLocus, ...]
    winners: Mapping[str, GroupHypothesis]
    accepted_proposals: tuple[JunctionProposal, ...]
    proposals: tuple[JunctionProposal, ...]
    counts: Mapping[str, int]


def discover_proposals(
    targets: Sequence[OneSidedTarget],
    routed: Mapping[str, Sequence[int]],
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    counts: dict[str, int],
) -> list[JunctionProposal]:
    return discover_proposal_block(targets, routed, reads, reference, family_kmers, counts)[0]


def discover_proposal_block(
    targets: Sequence[OneSidedTarget],
    routed: Mapping[str, Sequence[int]],
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    counts: dict[str, int],
    counter_start: int = 0,
) -> tuple[list[JunctionProposal], int]:
    """Run discovery over a target block; return the next query counter.

    Query IDs ``J%08d`` number every local proposal in target order, so a
    block starting at ``counter_start`` continues it.
    """

    discovery: list[JunctionProposal] = []
    counter = counter_start
    for target in targets:
        start, end = target.candidate_start, target.candidate_end
        middle = target.candidate_midpoint
        contig = reference[target.chrom]
        scope_start = max(0, int(math.floor(start)) - SCOPE_RADIUS)
        scope = contig[scope_start : min(len(contig), int(math.ceil(end)) + SCOPE_RADIUS + 1)]
        kmers = family_kmers.get(target.family, set())
        if not kmers:
            counts["NO_FAMILY_KMER_PANEL"] = counts.get("NO_FAMILY_KMER_PANEL", 0) + 1
            continue
        local: list[JunctionProposal] = []
        for pair_id in sorted(routed.get(target.locus_id, ())):
            mate1, mate2 = reads[pair_id]
            for mate, read in ((1, mate1), (2, mate2)):
                for junction in infer_read_junctions(read, kmers, scope, scope_start):
                    if junction.host_seed_first_pos is None:
                        continue
                    if (junction.host_seed_first_pos < start - CANDIDATE_SLOP
                            or junction.host_seed_first_pos > end + CANDIDATE_SLOP):
                        continue
                    view_sequence = read.upper() if junction.view == "FWD" else revcomp(read)
                    split = junction.split
                    if split < JUNCTION_EACH_SIDE or split + JUNCTION_EACH_SIDE > len(view_sequence):
                        continue
                    if junction.side == "HOST_LEFT_TE_RIGHT":
                        host_half = view_sequence[split - JUNCTION_EACH_SIDE : split]
                    else:
                        host_half = view_sequence[split : split + JUNCTION_EACH_SIDE]
                    signature = junction.signature.upper()
                    if len(signature) != 2 * JUNCTION_EACH_SIDE:
                        continue
                    score = (
                        100.0 * float(junction.te_seed_fraction)
                        + 5.0 * min(int(junction.te_seed_hits), 10)
                        + (10.0 if junction.host_seed_matches_in_scope == 1
                           else 1.0 / max(1, junction.host_seed_matches_in_scope))
                    )
                    local.append(JunctionProposal(
                        f"J{counter:08d}", target.locus_id, target.chrom, target.family, target.missing_side,
                        start, end, middle, pair_id, mate, junction.view, junction.side, split,
                        junction.te_bp, junction.te_seed_hits, junction.te_seed_fraction,
                        junction.host_seed_matches_in_scope, junction.host_seed_first_pos,
                        signature, min(signature, revcomp(signature)), host_half, score,
                    ))
                    counter += 1
        best: dict[str, JunctionProposal] = {}
        for proposal in local:
            old = best.get(proposal.canonical_signature)
            if old is None or proposal.strength_key < old.strength_key:
                best[proposal.canonical_signature] = proposal
        chosen = sorted(best.values(), key=lambda row: (row.strength_key, row.source_pair_id, row.source_mate))
        discovery.extend(chosen[:MAX_DISCOVERY_SIGNATURES_PER_LOCUS])
    return discovery, counter


def screen_proposals(
    discovery: Sequence[JunctionProposal],
    *,
    bowtie2: str | Path,
    host_index_prefix: str | Path,
    work_directory: Path,
    counts: dict[str, int],
) -> tuple[list[tuple[JunctionProposal, float, str]], dict[str, str]]:
    """Screen by exact-unique host half and contiguous-reference competition."""

    work_directory.mkdir(parents=True, exist_ok=True)
    host_fasta = work_directory / "junction_host_halves.fa"
    signature_fasta = work_directory / "junction_signatures.fa"
    host_sam = work_directory / "junction_host_halves.sam"
    signature_sam = work_directory / "junction_signatures.sam"
    host_alignments: dict[str, list[dict[str, object]]] = {}
    signature_alignments: dict[str, list[dict[str, object]]] = {}
    if discovery:
        write_fasta(host_fasta, [(row.qid, row.host_half) for row in discovery])
        write_fasta(signature_fasta, [(row.qid, row.signature) for row in discovery])
        for fasta, sam in ((host_fasta, host_sam), (signature_fasta, signature_sam)):
            run_bowtie2(bowtie2, [*EXACTISH_OPTIONS, "-k", str(HOST_PLACEMENT_REPORT), "-x", str(host_index_prefix),
                                  "-U", str(fasta), "-S", str(sam)])
        host_alignments = parse_sam_records(host_sam)
        signature_alignments = parse_sam_records(signature_sam)
    accepted: list[tuple[JunctionProposal, float, str]] = []
    statuses: dict[str, str] = {}
    for proposal in discovery:
        host = _unique_exact_host(host_alignments.get(proposal.qid, []))
        reference_nms = [int(row["nm"]) for row in signature_alignments.get(proposal.qid, []) if not row["supplementary"]]
        reference_nm = min(reference_nms) if reference_nms else None
        if host is None:
            status = "REJECT_HOST_NOT_GLOBALLY_EXACT_UNIQUE"
        else:
            host_start = int(host["pos0"])
            host_midpoint = (host_start + host_start + cigar_ref_span(str(host["cigar"]), len(proposal.host_half))) / 2.0
            near = host["rname"] == proposal.chrom and interval_distance_point(
                host_midpoint, proposal.candidate_start, proposal.candidate_end) <= CANDIDATE_SLOP
            if not near:
                status = "REJECT_HOST_EXACT_UNIQUE_OFF_NOMINATED_LOCUS"
            else:
                side = ("LEFT" if host_midpoint < proposal.candidate_midpoint
                        else "RIGHT" if host_midpoint > proposal.candidate_midpoint else None)
                if side is None:
                    status = "REJECT_JUNCTION_SIDE_UNRESOLVED"
                elif side != proposal.missing_side:
                    status = "REJECT_JUNCTION_NOT_MISSING_SIDE"
                elif reference_nm is not None and reference_nm <= REFERENCE_SIGNATURE_MAX_EDITS:
                    status = "REJECT_CONTIGUOUS_REFERENCE_EXPLANATION"
                else:
                    status = "R7_INSERTION_SPECIFIC_JUNCTION_ACCEPTED_PREQUERY"
                    accepted.append((proposal, host_midpoint, side))
        statuses[proposal.qid] = status
        counts[status] = counts.get(status, 0) + 1
    return accepted, statuses


LocatedJunction = tuple[JunctionProposal, float, str, int, tuple[int, ...]]


def locate_accepted(
    accepted: Sequence[tuple[JunctionProposal, float, str]], locator: SeedPairLocator
) -> list[LocatedJunction]:
    located: list[LocatedJunction] = []
    for proposal, host_midpoint, side in accepted:
        rows, complete, pair_ids = locator.seed_pairs(proposal.signature, SIGNATURE_LOCATE_CAP)
        if complete:
            located.append((proposal, host_midpoint, side, rows, tuple(sorted(pair_ids))))
    return located


def build_group_hypotheses(
    located: Sequence[LocatedJunction],
    hit_reads: Mapping[int, tuple[str, str]],
    reserved: set[int],
) -> dict[str, list[GroupHypothesis]]:
    members: dict[str, set[int]] = defaultdict(set)
    group_of: dict[int, str] = {}
    for pair_id, (mate1, mate2) in hit_reads.items():
        group = pair_group_hash(mate1, mate2)
        group_of[pair_id] = group
        members[group].add(pair_id)
    hypotheses: dict[str, list[GroupHypothesis]] = defaultdict(list)
    for proposal, host_midpoint, side, rows, pair_ids in located:
        seen: set[str] = set()
        for pair_id in pair_ids:
            group = group_of[pair_id]
            if group in seen:
                continue
            seen.add(group)
            group_members = tuple(sorted(members[group]))
            if set(group_members) & reserved:
                continue
            hypotheses[group].append(GroupHypothesis(
                group, group_members[0], group_members, proposal.locus_id, side, proposal.family,
                proposal.canonical_signature, proposal.source_pair_id, proposal.proposal_score,
                proposal.te_seed_hits, proposal.te_seed_fraction, rows,
                interval_distance_point(host_midpoint, proposal.candidate_start, proposal.candidate_end),
            ))
    return hypotheses


def compete_and_link(
    hypotheses: Mapping[str, Sequence[GroupHypothesis]],
    targets: Sequence[OneSidedTarget],
    strict_winners: Sequence[StrictPairWinner],
    counts: dict[str, int],
) -> tuple[dict[str, GroupHypothesis], tuple[JunctionLocus, ...]]:
    strict_by_side: dict[tuple[str, str], set[int]] = defaultdict(set)
    for winner in strict_winners:
        strict_by_side[(winner.locus_id, winner.side)].add(winner.pair_id)
    winners: dict[str, GroupHypothesis] = {}
    by_side: dict[tuple[str, str], list[GroupHypothesis]] = defaultdict(list)
    for group in sorted(hypotheses):
        winner, reason = resolve_group_hypotheses(hypotheses[group])
        counts[reason] = counts.get(reason, 0) + 1
        if winner is not None:
            winners[group] = winner
            by_side[(winner.locus_id, winner.side)].append(winner)
    ledger: list[JunctionLocus] = []
    for target in targets:
        strict_left = strict_by_side.get((target.locus_id, "LEFT"), set())
        strict_right = strict_by_side.get((target.locus_id, "RIGHT"), set())
        already = set(target.existing_left_ids) | set(target.existing_right_ids) | strict_left | strict_right
        linked: dict[str, tuple[list[str], list[int]]] = {"LEFT": ([], []), "RIGHT": ([], [])}
        for side in ("LEFT", "RIGHT"):
            for hypothesis in sorted(by_side.get((target.locus_id, side), []), key=lambda row: row.decision_key):
                group_members = set(hypothesis.group_pair_ids)
                if group_members & already:
                    continue
                representative = min(group_members)
                linked[side][0].append(hypothesis.group_hash)
                linked[side][1].append(representative)
                already.add(representative)
        ledger.append(JunctionLocus(
            target.locus_id, target.chrom, target.family, target.candidate_start, target.candidate_end,
            target.existing_left_ids, target.existing_right_ids,
            tuple(sorted(strict_left)), tuple(sorted(strict_right)),
            tuple(linked["LEFT"][0]), tuple(linked["RIGHT"][0]),
            tuple(sorted(linked["LEFT"][1])), tuple(sorted(linked["RIGHT"][1])),
        ))
    return winners, tuple(ledger)


def run_junction_groups(
    *,
    targets: Sequence[OneSidedTarget],
    routing: RoutingResult,
    strict_winners: Sequence[StrictPairWinner],
    reserved_call_pair_ids: set[int],
    reads: Mapping[int, tuple[str, str]],
    reference: Mapping[str, str],
    family_kmers: Mapping[str, set[str]],
    locator: SeedPairLocator,
    bowtie2: str | Path,
    host_index_prefix: str | Path,
    work_directory: Path,
    pair_reader,
) -> JunctionGroupsResult:
    """``pair_reader(pair_ids)`` fetches exact-hit pairs as ``{pair_id: (mate1, mate2)}``."""

    counts: dict[str, int] = {}
    target_ids = {target.locus_id for target in targets}
    routed: dict[str, list[int]] = defaultdict(list)
    for pair_id, locus_ids in routing.pair_to_loci.items():
        for locus_id in locus_ids:
            if locus_id in target_ids:
                routed[locus_id].append(pair_id)
    discovery = discover_proposals(targets, routed, reads, reference, family_kmers, counts)
    accepted, _ = screen_proposals(
        discovery, bowtie2=bowtie2, host_index_prefix=host_index_prefix, work_directory=work_directory, counts=counts
    )
    located = locate_accepted(accepted, locator)
    hit_union: set[int] = set()
    for row in located:
        hit_union.update(row[4])
    reserved = set(reserved_call_pair_ids) | {winner.pair_id for winner in strict_winners}
    hypotheses = build_group_hypotheses(located, pair_reader(hit_union), reserved)
    winners, ledger = compete_and_link(hypotheses, targets, strict_winners, counts)
    return JunctionGroupsResult(ledger, winners, tuple(row[0] for row in accepted), tuple(discovery), counts)
