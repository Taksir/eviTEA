"""Read length and fragment size of the sample, resolved once per run.

Read length
-----------

The library layer measures TE coverage in read-length windows.  The length comes, in this order, from

1. ``--read-length`` (``runtime.read_length``), when given;
2. the prepared sample (``backend_profile.uniform_length``, recorded by ``prepare-sample``);
3. the msBWT index itself: a uniform-length index of ``n`` reads holds ``n`` end markers and ``n * (L + 1)``
   symbols (``totalCounts.npy``), so ``L = total / n - 1``.

``DEFAULT_READ_LENGTH`` (125 bp) is the length of the development samples; it is used only when none of the above
is available (for example by direct callers of the library functions).

Fragment size
-------------
Two percentiles of the sample's fragment-size distribution are used by the library layer:

* the 95th percentile (``fragment_q95``) decides whether an ordinary pair on the reference explains a read pair: a
  pair whose span is above it is unlikely to be a normal pair, so it is not taken as one;
* the 99th percentile (``fragment_q99``) bounds the distance between a mate and the breakpoint it supports, and is
  Bowtie2's maximum fragment length (``-X``) in the paired reference mappings.

They are estimated from the sample's own properly paired reads (``estimate_fragment_sizes``): up to
``FRAGMENT_SAMPLE_PAIRS`` pairs, taken at evenly spaced pair ordinals, are mapped to the reference with Bowtie2
(end-to-end, ``--no-mixed --no-discordant``, ``-X FRAGMENT_SEARCH_MAX``).  Spans of pairs with both mates at
MAPQ >= ``FRAGMENT_MIN_MAPQ`` are kept; spans above median + ``FRAGMENT_OUTLIER_MADS`` x MAD (the median absolute
deviation) are dropped as outliers, as Picard's insert-size metrics do.  The percentiles are nearest-rank percentiles
of the remaining spans, so the same reads give the same numbers.  With fewer than ``MIN_FRAGMENT_PAIRS`` usable pairs
the development values (``DEFAULT_FRAGMENT_Q95`` 494 bp, ``DEFAULT_FRAGMENT_Q99`` 805 bp) are used and a warning is
logged.  ``--fragment-q95`` and ``--fragment-q99`` override the estimate.
"""

from __future__ import annotations

import logging
import math
import statistics
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .errors import InputContractError

__all__ = [
    "DEFAULT_FRAGMENT_Q95",
    "DEFAULT_FRAGMENT_Q99",
    "DEFAULT_READ_LENGTH",
    "FragmentSizes",
    "estimate_fragment_sizes",
    "fragment_sizes_from_spans",
    "nearest_rank_percentile",
    "properly_paired_spans",
    "resolve_fragment_sizes",
    "MIN_READ_LENGTH",
    "ReadLength",
    "read_length_from_index",
    "resolve_read_length",
    "validate_read_length",
]

DEFAULT_READ_LENGTH = 125
# Shortest read the window geometry accepts: the library layer's longest exact word is 21 bp.
MIN_READ_LENGTH = 21

DEFAULT_FRAGMENT_Q95 = 494.0
DEFAULT_FRAGMENT_Q99 = 805
FRAGMENT_SAMPLE_PAIRS = 200_000
MIN_FRAGMENT_PAIRS = 1_000
FRAGMENT_SEARCH_MAX = 5_000
FRAGMENT_MIN_MAPQ = 20
FRAGMENT_OUTLIER_MADS = 10
FRAGMENT_MAPPING_OPTIONS = (
    "--end-to-end", "--no-mixed", "--no-discordant", "-I", "0", "-X", str(FRAGMENT_SEARCH_MAX),
    "--no-unal", "--no-hd", "--reorder",
)

_log = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class ReadLength:
    length: int
    source: str  # "option", "prepared-sample", "msbwt-index" or "default"


def validate_read_length(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < MIN_READ_LENGTH:
        raise InputContractError(f"read length must be an integer of at least {MIN_READ_LENGTH} bp, got {value!r}")
    return value


def read_length_from_index(index_dir: str | Path) -> int:
    """Read length of a uniform-length msBWT index, from its symbol counts."""

    import numpy as np

    path = Path(index_dir) / "totalCounts.npy"
    try:
        counts = np.load(path, allow_pickle=False)
    except (OSError, ValueError) as exc:
        raise InputContractError(f"cannot read the msBWT symbol counts: {path}") from exc
    reads, total = int(counts[0]), int(counts.sum())
    if reads <= 0 or total % reads:
        raise InputContractError(
            f"the msBWT index {index_dir} does not hold reads of one length; pass --read-length"
        )
    return validate_read_length(total // reads - 1)


def resolve_read_length(
    option: int | None,
    backend_profile: Mapping[str, Any] | None,
    msbwt_index: str | Path | None,
) -> ReadLength:
    if option is not None:
        return ReadLength(validate_read_length(option), "option")
    if backend_profile is not None and backend_profile.get("uniform_length") is not None:
        return ReadLength(validate_read_length(int(backend_profile["uniform_length"])), "prepared-sample")
    if msbwt_index is not None and (Path(msbwt_index) / "totalCounts.npy").is_file():
        return ReadLength(read_length_from_index(msbwt_index), "msbwt-index")
    return ReadLength(DEFAULT_READ_LENGTH, "default")


@dataclass(frozen=True, slots=True)
class FragmentSizes:
    q95: float
    q99: int
    pairs: int  # usable properly paired reads the estimate came from (0 when not estimated)
    source: str  # "estimate", "option", "estimate+option" or "default" (too few pairs, or not estimated)

    def record(self) -> dict[str, object]:
        return {"fragment_q95": self.q95, "fragment_q99": self.q99, "fragment_pairs": self.pairs,
                "fragment_source": self.source}


def validate_fragment_size(value: Any, name: str) -> Any:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0 or not math.isfinite(value):
        raise InputContractError(f"{name} must be a positive number, got {value!r}")
    return value


def nearest_rank_percentile(sorted_values: list[int], percent: float) -> int:
    """The smallest value with at least ``percent`` % of the values at or below it."""

    if not sorted_values:
        raise ValueError("no values")
    rank = max(1, math.ceil(percent / 100 * len(sorted_values)))
    return sorted_values[rank - 1]


def properly_paired_spans(sam_lines: Iterable[str], min_mapq: int = FRAGMENT_MIN_MAPQ) -> list[int]:
    """Fragment spans (|TLEN|) of properly paired, primary pairs with both mates at ``min_mapq`` or above."""

    first: dict[str, tuple[int, int]] = {}
    spans: list[int] = []
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t", 9)
        flag, mapq, tlen = int(fields[1]), int(fields[4]), abs(int(fields[8]))
        if flag & 0x904 or not flag & 0x2 or not flag & 0x1 or tlen == 0:
            continue
        seen = first.pop(fields[0], None)
        if seen is None:
            first[fields[0]] = (mapq, tlen)
        elif min(seen[0], mapq) >= min_mapq:
            spans.append(tlen)
    return spans


def fragment_sizes_from_spans(spans: Iterable[int], *, min_pairs: int = MIN_FRAGMENT_PAIRS) -> FragmentSizes:
    values = sorted(int(x) for x in spans)
    if len(values) < min_pairs:
        return FragmentSizes(DEFAULT_FRAGMENT_Q95, DEFAULT_FRAGMENT_Q99, len(values), "default")
    median = statistics.median(values)
    mad = statistics.median(abs(x - median) for x in values)
    kept = [x for x in values if x <= median + FRAGMENT_OUTLIER_MADS * mad] if mad > 0 else values
    if len(kept) < min_pairs:
        return FragmentSizes(DEFAULT_FRAGMENT_Q95, DEFAULT_FRAGMENT_Q99, len(kept), "default")
    return FragmentSizes(float(nearest_rank_percentile(kept, 95)), nearest_rank_percentile(kept, 99), len(kept), "estimate")


def estimate_fragment_sizes(
    sample,
    *,
    bowtie2: str | Path,
    host_index_prefix: str | Path,
    work_directory: Path,
    threads: int = 1,
    sample_pairs: int = FRAGMENT_SAMPLE_PAIRS,
    min_pairs: int = MIN_FRAGMENT_PAIRS,
) -> FragmentSizes:
    """Estimate from ``sample`` (an ``incumbent.MsbwtSample``): see the module notes."""

    from .incumbent.external import run_bowtie2

    total = int(sample.pair_count)
    count = min(total, sample_pairs)
    pair_ids = sorted({i * total // count for i in range(count)}) if count else []
    work = Path(work_directory)
    work.mkdir(parents=True, exist_ok=True)
    read1, read2, sam = work / "fragment_R1.fa", work / "fragment_R2.fa", work / "fragment.sam"
    reads = sample.pair_reads(pair_ids)
    with read1.open("w", encoding="ascii", newline="\n") as a, read2.open("w", encoding="ascii", newline="\n") as b:
        for pid in pair_ids:
            first, second = reads[pid]
            a.write(f">P{pid}\n{first}\n")
            b.write(f">P{pid}\n{second}\n")
    run_bowtie2(bowtie2, ["-x", str(host_index_prefix), "-f", "-1", str(read1), "-2", str(read2),
                          *FRAGMENT_MAPPING_OPTIONS, "-p", str(max(1, threads)), "-S", str(sam)])
    with sam.open("r", encoding="ascii") as handle:
        sizes = fragment_sizes_from_spans(properly_paired_spans(handle), min_pairs=min_pairs)
    for path in (read1, read2, sam):
        path.unlink()
    return sizes


def resolve_fragment_sizes(
    q95_option: float | None,
    q99_option: int | None,
    estimate: Callable[[], FragmentSizes],
    *,
    warn: Callable[[str], None] | None = None,
) -> FragmentSizes:
    """Options win; the estimate (run only when an option is missing) fills the rest."""

    if q95_option is not None and q99_option is not None:
        return FragmentSizes(float(q95_option), int(q99_option), 0, "option")
    estimated = estimate()
    if estimated.source == "default":
        message = (
            f"fragment size: only {estimated.pairs} usable properly paired reads (at least {MIN_FRAGMENT_PAIRS} "
            f"needed); using the fallback 95th / 99th percentiles {DEFAULT_FRAGMENT_Q95:g} / {DEFAULT_FRAGMENT_Q99} bp"
        )
        _log.warning(message)
        if warn is not None:
            warn(message)
    if q95_option is None and q99_option is None:
        return estimated
    return FragmentSizes(
        float(q95_option) if q95_option is not None else estimated.q95,
        int(q99_option) if q99_option is not None else estimated.q99,
        estimated.pairs,
        f"{estimated.source}+option",
    )
