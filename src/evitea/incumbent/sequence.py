"""Sequence, FASTA and interval helpers."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

from ..errors import InputContractError

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

__all__ = [
    "revcomp",
    "shannon_entropy_acgt",
    "shannon_entropy_symbols",
    "count_occurrences",
    "kmer_hit_prefix",
    "interval_hits",
    "fasta_records_full_header",
    "load_reference_contigs",
    "interval_distance_point",
    "distance_to_interval",
    "interval_gap",
]


def revcomp(sequence: str) -> str:
    """Return the upper-case reverse complement; ``N`` stays ``N``."""

    return sequence.translate(_COMPLEMENT)[::-1].upper()


def shannon_entropy_acgt(sequence: str) -> float:
    """Return entropy over A/C/G/T counts only (k-mer panel filter)."""

    if not sequence:
        return 0.0
    counts = [sequence.count(base) for base in "ACGT"]
    total = sum(counts)
    if total == 0:
        return 0.0
    value = 0.0
    for count in counts:
        if count:
            p = count / total
            value -= p * math.log2(p)
    return value


def shannon_entropy_symbols(sequence: str) -> float:
    """Return entropy over all symbols, normalized by sequence length (host seeds)."""

    if not sequence:
        return 0.0
    counts = Counter(sequence.upper())
    total = len(sequence)
    value = 0.0
    for count in counts.values():
        p = count / total
        value -= p * math.log2(p)
    return value


def count_occurrences(
    text: str, pattern: str, stop_after: int | None = None
) -> tuple[int, int | None]:
    """Return overlapping occurrence count (optionally capped) and first position."""

    count = 0
    first: int | None = None
    position = 0
    while True:
        index = text.find(pattern, position)
        if index < 0:
            break
        if first is None:
            first = index
        count += 1
        if stop_after is not None and count >= stop_after:
            break
        position = index + 1
    return count, first


def kmer_hit_prefix(sequence: str, kmers: set[str] | frozenset[str], k: int) -> list[int]:
    """Return prefix sums of k-mer panel membership over k-mer start positions."""

    n = max(0, len(sequence) - k + 1)
    prefix = [0] * (n + 1)
    for index in range(n):
        prefix[index + 1] = prefix[index] + (1 if sequence[index : index + k] in kmers else 0)
    return prefix


def interval_hits(prefix: list[int], start: int, end: int, k: int) -> tuple[int, int]:
    """Return hits and total k-mer starts inside ``sequence[start:end]``."""

    if end - start < k:
        return 0, 0
    last = end - k + 1
    return prefix[last] - prefix[start], last - start


def fasta_records_full_header(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(stripped header line, upper-case sequence)`` records."""

    name: str | None = None
    chunks: list[str] = []
    with Path(path).open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def load_reference_contigs(path: str | Path, wanted: set[str]) -> dict[str, str]:
    """Load selected contigs by first header token, tolerating a ``chr`` prefix."""

    out: dict[str, str] = {}
    if not wanted:
        return out
    for header, sequence in fasta_records_full_header(path):
        token = header.split()[0] if header.split() else header
        for name in list(wanted):
            if (
                token == name
                or header == name
                or token.removeprefix("chr") == name.removeprefix("chr")
            ):
                out[name] = sequence
        if len(out) == len(wanted):
            break
    missing = sorted(wanted - set(out))
    if missing:
        raise InputContractError(f"reference FASTA lacks selected contigs: {missing[:20]}")
    return out


def interval_distance_point(x: float, start: float, end: float) -> float:
    """Return distance from a point to a closed interval ``start <= end``."""

    if start <= x <= end:
        return 0.0
    return min(abs(x - start), abs(x - end))


def distance_to_interval(x: float, a: float, b: float) -> float:
    """Return distance from a point to the closed interval spanned by ``a`` and ``b``."""

    lo, hi = sorted((float(a), float(b)))
    if x < lo:
        return lo - x
    if x > hi:
        return x - hi
    return 0.0


def interval_gap(a_start: float, a_end: float, b_start: float, b_end: float) -> float:
    """Return gap between closed intervals; zero when they touch or overlap."""

    a0, a1 = sorted((float(a_start), float(a_end)))
    b0, b1 = sorted((float(b_start), float(b_end)))
    if not (a1 < b0 or b1 < a0):
        return 0.0
    if a1 < b0:
        return b0 - a1
    return a0 - b1
