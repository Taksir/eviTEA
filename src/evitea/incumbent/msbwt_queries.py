"""Bounded exact msBWT queries with physical-pair provenance (incumbent path).

The sample index is the multi-string BWT built from the synchronized mate
FASTQs: ``about[dollar_id] = (file_id, record_ordinal)`` with file 0 = mate 1
and file 1 = mate 2, so ``record_ordinal`` is the physical ``pair_id``.  Every
quantity the incumbent derives from the index (interval widths, the set of
reads and pairs carrying a query, recovered read sequences) is invariant to the
index construction order, so any msBWT over the same reads reproduces it.
"""

from __future__ import annotations

from collections.abc import Iterable, MutableMapping
from pathlib import Path
from typing import Any

import numpy as np

from ..errors import AssetCompatibilityError, BackendCompatibilityError
from .mate_rescue import SeedPairLocator
from .sequence import revcomp

__all__ = ["MsbwtSample", "MsbwtPairLocator"]

_SCAN_CHUNK = 2_000_000


class MsbwtSample(SeedPairLocator):
    """Exact interval, dollar-ID, provenance and read-recovery access."""

    def __init__(
        self,
        index_directory: str | Path,
        *,
        bwt: Any = None,
        about: Any = None,
        interval_cache: MutableMapping[str, tuple[int, int]] | None = None,
        dollar_cache: MutableMapping[str, np.ndarray] | None = None,
    ):
        directory = Path(index_directory)
        if bwt is None:
            try:
                from MUSCython import MultiStringBWTCython
            except ImportError as exc:  # pragma: no cover - environment dependent
                raise BackendCompatibilityError("msbwt-modern3 is not installed") from exc
            bwt = MultiStringBWTCython.loadBWT(str(directory), useMemmap=True)
        if about is None:
            about = np.load(directory / "about.npy", mmap_mode="r", allow_pickle=False)
        names = getattr(about.dtype, "names", None)
        if about.ndim != 1 or not names or len(names) < 2:
            raise AssetCompatibilityError("msBWT about.npy must be a (file_id, record_ordinal) table")
        self.bwt = bwt
        self.about = about
        self.file_field, self.seq_field = names[0], names[1]
        self.record_count = int(about.shape[0])
        if self.record_count % 2:
            raise AssetCompatibilityError("msBWT about.npy must index both mates of every pair")
        self.pair_count = self.record_count // 2
        self._intervals: MutableMapping[str, tuple[int, int]] = {} if interval_cache is None else interval_cache
        self._dollars: MutableMapping[str, np.ndarray] = {} if dollar_cache is None else dollar_cache

    # -- exact queries -------------------------------------------------
    def interval(self, pattern: str, *, cache: bool = True) -> tuple[int, int]:
        hit = self._intervals.get(pattern)
        if hit is None:
            lo, hi = self.bwt.findIndicesOfStr(pattern)
            hit = (int(lo), int(hi))
            if hit[0] < 0 or hit[1] < hit[0]:
                raise BackendCompatibilityError(f"msBWT returned an invalid interval for {pattern}")
            if cache:
                self._intervals[pattern] = hit
        return hit

    def width(self, pattern: str) -> int:
        lo, hi = self.interval(pattern)
        return hi - lo

    def dollar_ids(self, pattern: str, *, cache: bool = True) -> np.ndarray:
        """Dollar IDs of every interval row of ``pattern`` (row order)."""

        cached = self._dollars.get(pattern)
        if cached is not None:
            return cached
        lo, hi = self.interval(pattern, cache=cache)
        values = np.fromiter(
            (int(self.bwt.getSequenceDollarID(row)) for row in range(lo, hi)), dtype=np.uint64, count=hi - lo
        )
        if cache:
            self._dollars[pattern] = values
        return values

    # -- provenance ----------------------------------------------------
    def file_and_pair(self, dollars: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if len(dollars) == 0:
            return np.zeros(0, dtype=np.uint8), np.zeros(0, dtype=np.uint64)
        rows = self.about[np.asarray(dollars, dtype=np.uint64)]
        files = np.asarray(rows[self.file_field], dtype=np.uint8)
        pairs = np.asarray(rows[self.seq_field], dtype=np.uint64)
        if np.any(files > 1):
            raise AssetCompatibilityError("msBWT about.npy contains a non-paired file id")
        return files, pairs

    def verify_paired_layout(self) -> None:
        """Fail closed unless files 0/1 each index record ordinals ``0..pairs-1``."""

        counts = np.zeros(2, dtype=np.int64)
        lo = [None, None]
        hi = [None, None]
        for start in range(0, self.record_count, _SCAN_CHUNK):
            chunk = self.about[start : start + _SCAN_CHUNK]
            files = np.asarray(chunk[self.file_field], dtype=np.uint8)
            seqs = np.asarray(chunk[self.seq_field], dtype=np.uint64)
            if np.any(files > 1):
                raise AssetCompatibilityError("msBWT about.npy contains a non-paired file id")
            for fid in (0, 1):
                values = seqs[files == fid]
                if len(values):
                    counts[fid] += len(values)
                    mn, mx = int(values.min()), int(values.max())
                    lo[fid] = mn if lo[fid] is None else min(lo[fid], mn)
                    hi[fid] = mx if hi[fid] is None else max(hi[fid], mx)
        if counts[0] != self.pair_count or counts[1] != self.pair_count:
            raise AssetCompatibilityError("msBWT about.npy mate counts differ")
        if lo != [0, 0] or hi != [self.pair_count - 1, self.pair_count - 1]:
            raise AssetCompatibilityError("msBWT about.npy record ordinals are not 0..pairs-1")

    def mate_dollars(self, wanted: dict[int, int]) -> dict[int, int]:
        """Map ``pair_id -> file_id`` requests to the dollar ID of that mate."""

        if not wanted:
            return {}
        ids = np.asarray(sorted(wanted), dtype=np.uint64)
        want_file = np.asarray([wanted[int(x)] for x in ids.tolist()], dtype=np.uint8)
        found: dict[int, int] = {}
        for start in range(0, self.record_count, _SCAN_CHUNK):
            chunk = self.about[start : start + _SCAN_CHUNK]
            seqs = np.asarray(chunk[self.seq_field], dtype=np.uint64)
            files = np.asarray(chunk[self.file_field], dtype=np.uint8)
            idx = np.searchsorted(ids, seqs)
            safe = np.minimum(idx, len(ids) - 1)
            hit = (idx < len(ids)) & (ids[safe] == seqs) & (want_file[safe] == files)
            for j in np.flatnonzero(hit).tolist():
                pair_id = int(seqs[j])
                if pair_id in found:
                    raise AssetCompatibilityError(f"duplicate mate record for pair {pair_id}")
                found[pair_id] = start + int(j)
        missing = set(wanted) - set(found)
        if missing:
            raise AssetCompatibilityError(f"msBWT lacks {len(missing)} requested mates, e.g. {sorted(missing)[:5]}")
        return found

    def recover(self, dollar_id: int) -> str:
        raw = self.bwt.recoverString(int(dollar_id))
        if isinstance(raw, (bytes, bytearray)):
            raw = bytes(raw).decode("ascii")
        text = str(raw).strip()
        if text.startswith("$"):
            text = text[1:]
        text = text.upper()
        if not text or any(base not in "ACGNT" for base in text):
            raise BackendCompatibilityError(f"msBWT recovered an invalid read for dollar {dollar_id}")
        return text

    def pair_reads(self, pair_ids: Iterable[int]) -> dict[int, tuple[str, str]]:
        """``{pair_id: (mate1, mate2)}`` recovered from the index."""

        pairs = sorted({int(x) for x in pair_ids})
        if not pairs:
            return {}
        mate1 = self.mate_dollars({pid: 0 for pid in pairs})
        mate2 = self.mate_dollars({pid: 1 for pid in pairs})
        return {pid: (self.recover(mate1[pid]), self.recover(mate2[pid])) for pid in pairs}

    # -- bounded seed location (R6 routing) ------------------------------
    def seed_pairs(self, seed: str, row_cap: int) -> tuple[int, bool, set[int]]:
        """Rows over seed+RC, whether located, and physical pair IDs if located."""

        patterns = [seed.upper()]
        reverse = revcomp(seed.upper())
        if reverse != patterns[0]:
            patterns.append(reverse)
        intervals = [self.interval(pattern, cache=False) for pattern in patterns]
        total = sum(hi - lo for lo, hi in intervals)
        if total > row_cap:
            return total, False, set()
        dollars = np.unique(np.concatenate([self.dollar_ids(p, cache=False) for p in patterns]))
        files, pairs = self.file_and_pair(dollars)
        return total, True, {int(x) for x in pairs[files <= 1].tolist()}


# Backwards-compatible name used by the R6 routing tests.
MsbwtPairLocator = MsbwtSample
