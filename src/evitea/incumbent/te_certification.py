"""Confirm junction group families from full TE-facing segments.

Only junction groups linked into a locus aggregation list ("used" groups) are
confirmed.  For each, the untrimmed TE-facing segment is taken from the
winning junction's source read, aligned end-to-end anywhere in the TE library
(either strand) with the built-in exhaustive Bowtie2 profile, collapsed to the
minimum NM per alias-normalized biological family, and accepted when the
candidate family is the unique best family and ``35 * NM <= 3 * length``.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from .external import TE_LIBRARY_OPTIONS, parse_sam_records, run_bowtie2, write_fasta
from .junction_groups import GroupHypothesis, JunctionLocus, JunctionProposal
from .sequence import revcomp
from .te_library import normalize_family

MIN_TE_BP = 35
RATE_NUM = 3
RATE_DEN = 35
T0_PREFIX = "T0__"


@dataclass(frozen=True, slots=True)
class UsedGroup:
    group_hash: str
    locus_id: str
    side: str
    candidate_family: str
    source_pair_id: int
    source_mate: int
    read_view: str
    junction_semantic: str
    split: int
    te_bp: int
    representative_pair_id: int


def used_groups(
    ledger: Sequence[JunctionLocus],
    winners: Mapping[str, GroupHypothesis],
    accepted: Sequence[JunctionProposal],
) -> tuple[UsedGroup, ...]:
    used: dict[str, tuple[str, str]] = {}
    for row in ledger:
        for side, hashes in (("LEFT", row.left_group_hashes), ("RIGHT", row.right_group_hashes)):
            for group in hashes:
                value = (row.locus_id, side)
                if group in used and used[group] != value:
                    raise ValueError(f"group {group} is linked to several locus sides")
                used[group] = value
    index: dict[tuple[str, str, int], list[JunctionProposal]] = defaultdict(list)
    for proposal in accepted:
        index[(proposal.locus_id, proposal.canonical_signature, proposal.source_pair_id)].append(proposal)
    out: list[UsedGroup] = []
    for group, (locus_id, side) in sorted(used.items()):
        winner = winners[group]
        rows = sorted(
            index[(locus_id, winner.canonical_signature, winner.source_pair_id)],
            key=lambda row: (-row.proposal_score, -row.te_seed_fraction, -row.te_seed_hits, row.source_mate, row.split),
        )
        if not rows:
            raise ValueError(f"used group {group} has no accepted source junction")
        source = rows[0]
        out.append(UsedGroup(
            group, locus_id, side, normalize_family(winner.family), source.source_pair_id,
            source.source_mate, source.read_view, source.junction_semantic, source.split, source.te_bp,
            winner.representative_pair_id,
        ))
    return tuple(out)


def full_te_segment(read: str, read_view: str, junction_semantic: str, split: int) -> str:
    sequence = read.upper() if read_view == "FWD" else revcomp(read.upper())
    if not 0 < split < len(sequence):
        raise ValueError(f"invalid split {split} for read length {len(sequence)}")
    if junction_semantic == "HOST_LEFT_TE_RIGHT":
        return sequence[split:]
    if junction_semantic == "TE_LEFT_HOST_RIGHT":
        return sequence[:split]
    raise ValueError(f"unknown junction semantic {junction_semantic}")


@dataclass(frozen=True, slots=True)
class TEAlignment:
    record_id: str
    reverse: bool
    pos1: int
    cigar: str
    sequence: str
    nm: int


def align_segments(
    segments: Mapping[str, str],
    *,
    bowtie2: str | Path,
    te_index_prefix: str | Path,
    work_directory: Path,
) -> dict[str, list[TEAlignment]]:
    """Align untrimmed segments (query ``T0__<group>``) and keep scored rows."""

    work_directory.mkdir(parents=True, exist_ok=True)
    fasta = work_directory / "te_segments.fa"
    sam = work_directory / "te_segments.sam"
    records = [(T0_PREFIX + group, sequence) for group, sequence in segments.items() if len(sequence) >= MIN_TE_BP]
    out: dict[str, list[TEAlignment]] = defaultdict(list)
    if not records:
        return out
    write_fasta(fasta, records)
    run_bowtie2(bowtie2, [*TE_LIBRARY_OPTIONS, "-x", str(te_index_prefix), "-U", str(fasta), "-S", str(sam)])
    for query, rows in parse_sam_records(sam).items():
        if not query.startswith(T0_PREFIX):
            continue
        for row in rows:
            if row["supplementary"] or int(row["nm"]) == 999999:
                continue
            out[query[len(T0_PREFIX):]].append(TEAlignment(
                str(row["rname"]), bool(row["rev"]), int(row["pos0"]) + 1, str(row["cigar"]),
                str(row["seq"]).upper(), int(row["nm"]),
            ))
    return out


def family_min_nm(alignments: Sequence[TEAlignment], record_families: Mapping[str, str]) -> dict[str, int]:
    best: dict[str, int] = {}
    for row in alignments:
        family = record_families.get(row.record_id)
        if family is None:
            continue
        old = best.get(family)
        if old is None or row.nm < old:
            best[family] = row.nm
    return best


def unique_best_rate(family_nm: Mapping[str, int], candidate: str, length: int) -> tuple[bool, str, int | None]:
    if length < MIN_TE_BP:
        return False, "REJECT_TE_FACING_SEGMENT_LT35", None
    if not family_nm:
        return False, "REJECT_NO_TE_LIBRARY_ALIGNMENT", None
    best = min(family_nm.values())
    winners = sorted(family for family, nm in family_nm.items() if nm == best)
    if len(winners) != 1:
        return False, "REJECT_EQUAL_BEST_BIOLOGICAL_FAMILY_TIE", best
    if winners[0] != normalize_family(candidate):
        return False, "REJECT_UNIQUE_BEST_DIFFERENT_FAMILY", best
    if RATE_DEN * best > RATE_NUM * length:
        return False, "REJECT_UNIQUE_BEST_BUT_RATE_ABOVE_3_PER_35", best
    return True, "PASS_UNIQUE_BEST_RATE_3_PER_35", best


@dataclass(frozen=True, slots=True)
class GroupCertification:
    group: UsedGroup
    segment_length: int
    passes: bool
    status: str
    best_nm: int | None
    alignments: tuple[TEAlignment, ...]


def certify_groups(
    groups: Sequence[UsedGroup],
    segments: Mapping[str, str],
    alignments: Mapping[str, Sequence[TEAlignment]],
    record_families: Mapping[str, str],
) -> tuple[GroupCertification, ...]:
    out: list[GroupCertification] = []
    for group in groups:
        length = len(segments[group.group_hash])
        rows = tuple(alignments.get(group.group_hash, ()))
        passes, status, best = unique_best_rate(family_min_nm(rows, record_families), group.candidate_family, length)
        out.append(GroupCertification(group, length, passes, status, best, rows))
    return tuple(out)
