"""The reference word tables of the detection, kept on disk so that a job loads them instead of building them again.

``detector.reference_multi`` (the multi-copy 21- and 31-letter words) and ``walks.RefWords`` (every 31-letter word) depend only on the
reference and on the detection code. A cache folder holds them under a key of those two inputs: ``<cache_dir>/<key>/`` with
``m21.npy``, ``m31.npy``, ``words31.npy`` (when the walk-out needs it) and ``about.json``. The arrays are only searched with
``np.searchsorted``, so a memory-mapped array answers exactly as the built one, and worker processes share its pages.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
from collections.abc import Callable, Mapping
from pathlib import Path

import numpy as np

from . import _p4_lib as P
from .walks import RefWords

__all__ = ["SCHEMA", "cache_key", "load_or_build"]

SCHEMA = 1
M21, M31, WORDS, ABOUT = "m21.npy", "m31.npy", "words31.npy", "about.json"


def cache_key(reference_sha256: str, code_sha256: str) -> str:
    """sha256 of the inputs the tables depend on."""

    text = json.dumps({"schema": SCHEMA, "code_sha256": code_sha256, "reference_sha256": reference_sha256}, sort_keys=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _save(path: Path, array: np.ndarray) -> None:
    with path.open("wb") as handle:
        np.save(handle, array)


def _load(path: Path) -> np.ndarray:
    """Memory-mapped, except an empty array (it cannot be mapped)."""

    array = np.load(path, mmap_mode="r")
    if array.size == 0:
        array = np.load(path)
    if array.dtype != np.uint64:
        raise ValueError(f"{path}: expected uint64 codes, found {array.dtype}")
    return array


def _words_object(codes: np.ndarray) -> RefWords:
    words = RefWords.__new__(RefWords)
    words.codes = codes
    return words


def _write_words(folder: Path, genome: Mapping[str, str]) -> tuple[np.ndarray, float]:
    started = time.time()
    codes = RefWords(genome).codes
    tmp = folder / f"{WORDS}.tmp.{os.getpid()}"
    _save(tmp, codes)
    tmp.replace(folder / WORDS)
    return codes, round(time.time() - started, 1)


def load_or_build(genome: Mapping[str, str], *, reference_sha256: str, code_sha256: str, cache_dir: Path, need_words: bool,
                  progress: Callable[[str], None]) -> tuple[P.RefMulti, RefWords | None]:
    """(``RefMulti``, ``RefWords`` or None) of the reference. Loaded from ``cache_dir`` when the folder of the key exists, else built
    and saved (a temporary folder renamed at the end, so a half-written cache is never read). ``need_words`` False builds or loads
    no word table; a folder written without it gets ``words31.npy`` when a later run needs it."""

    cache_dir = Path(cache_dir)
    key = cache_key(reference_sha256, code_sha256)
    folder = cache_dir / key
    if not (folder / M21).is_file() or not (folder / M31).is_file():
        started = time.time()
        cache_dir.mkdir(parents=True, exist_ok=True)
        tmp = cache_dir / f"{key}.tmp.{os.getpid()}"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir()
        progress("library-free: reference multi-copy words")
        m21 = P.multi_codes(genome, P.REP_K)
        m31 = P.multi_codes(genome, P.K)
        about = {"schema": SCHEMA, "code_sha256": code_sha256, "reference_sha256": reference_sha256,
                 "m21_words": int(len(m21)), "m31_words": int(len(m31)), "multi_seconds": round(time.time() - started, 1)}
        _save(tmp / M21, m21)
        _save(tmp / M31, m31)
        (tmp / ABOUT).write_text(json.dumps(about, indent=2), encoding="utf-8")
        try:
            tmp.rename(folder)
            progress(f"library-free: reference word tables saved in {folder} ({about['multi_seconds']} s)")
        except OSError:  # another job renamed its folder first: its tables are the same
            shutil.rmtree(tmp, ignore_errors=True)
            if not (folder / M21).is_file():
                raise
    else:
        progress(f"library-free: reference multi-copy words loaded from {folder}")
    rm = P.RefMulti(_load(folder / M21), _load(folder / M31))
    if not need_words:
        return rm, None
    if not (folder / WORDS).is_file():
        progress("library-free: reference words for the walk-out")
        codes, seconds = _write_words(folder, genome)
        about_path = folder / ABOUT
        try:
            about = json.loads(about_path.read_text(encoding="utf-8"))
            about.update(words31=int(len(codes)), words_seconds=seconds)
            tmp_about = folder / f"{ABOUT}.tmp.{os.getpid()}"
            tmp_about.write_text(json.dumps(about, indent=2), encoding="utf-8")
            tmp_about.replace(about_path)
        except (OSError, ValueError):
            pass
    else:
        progress("library-free: reference words for the walk-out loaded from the cache")
    return rm, _words_object(_load(folder / WORDS))
