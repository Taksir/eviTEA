"""One contig of the detection scanned by several processes, with the places of ``scan.scan_contig_passes``.

The scan of a contig has three parts. The first look counts the 31-letter words (and, with the haplotype scan, the 41-letter words)
at stride 10: the counts of one word do not depend on any other. The candidate runs follow from the whole contig's counts. Each run
is then walked on its own: a place depends on its run and on values of the whole contig (local depth, depth floor, the 41/31 ratio,
the repeat cutoff), never on another place. The only state shared between places is a count cache, which changes how often a count is
computed and never its value.

So the contig is cut in two ways, and the pieces are joined in order:

* the first look, in chunks of word positions: each worker counts the words of its chunk, the main process joins the arrays and runs
  the rest of ``scan_contig_passes`` that is not a walk (local depth, runs, the haplotype-scan masks, the depth floor);
* the runs, in tiles of ``tile`` letters: a tile holds the first-pass runs and the haplotype-scan runs that begin in it. A worker
  walks the first-pass runs of its tile and then the haplotype-scan runs of the same tile with one counter, so the second pass finds
  the counts of the first. A walk may reach past the tile edge: every worker holds the whole contig.

Workers are forked, so they share the msBWT pages, the reference tables, the minimap2 index and the contig.
"""

from __future__ import annotations

import gc
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import _p4_lib as P
from . import scan as S
from .options import Options

__all__ = ["S1_CHUNK", "TILE", "resolve_workers", "scan_contig_tiled"]

TILE = 250_000        # letters of a contig in one tile
S1_CHUNK = 20_000     # word positions in one chunk of the first look
_SHARED: dict = {}    # the per-contig inputs; set before the pool forks


def resolve_workers(workers: int) -> int:
    """``--workers``: 0 = every CPU this process may use."""

    if workers < 0:
        raise ValueError("--workers must be 0 (all CPUs) or a positive number")
    return workers or max(1, os.process_cpu_count() or 1)


def _join_stats(a: dict, b: dict) -> dict:
    """The stats of two tiles added: numbers summed, nested dicts joined the same way, ``rep_max`` and ``mode`` copied (they must
    agree). Keys keep the order of first appearance."""

    if list(a) != list(b):
        raise RuntimeError(f"tile stats differ in their keys: {list(a)} and {list(b)}")
    out: dict = {}
    for key, value in a.items():
        other = b[key]
        if isinstance(value, dict):
            out[key] = _join_stats(value, other)
        elif key in ("rep_max", "mode"):
            if value != other:
                raise RuntimeError(f"tile stats disagree on {key}: {value} and {other}")
            out[key] = value
        else:
            if not (type(value) is int and type(other) is int):
                raise RuntimeError(f"tile stat {key} is not a whole number: {value!r}, {other!r}")
            out[key] = value + other
    return out


def _s1_task(index: int):
    """Counts of word positions [ja, jb) of the first look (and the 41-letter counts with the haplotype scan)."""

    sh = _SHARED
    ja, jb = sh["s1_bounds"][index]
    counter = sh["make_counter"](False)
    ref = sh["ref"]
    _, c = P.scan_counts(counter, ref[ja * P.STRIDE:(jb - 1) * P.STRIDE + P.K])
    x = S.count_x(counter, ref, sh["pos"][ja:jb]) if sh["het"] else None
    return c, x


def _tile_task(index: int):
    """The first-pass places and, with the haplotype scan, its places for the runs that begin in tile ``index``."""

    sh = _SHARED
    counter = sh["make_counter"](True)
    (f0, f1), (h0, h1) = sh["first_range"][index], sh["het_range"][index]
    info: dict = {}
    first, _, _, xfirst = S._first_pass(counter, sh["ref"], sh["chrom"], sh["rm"], sh["firewall"], sh["opts"], sh["words"],
                                        sh["cands"][f0:f1], info, sh["record"], sh["d_min"], sh["cross"])
    second, het_info = None, {}
    if sh["het"]:
        second = S._het_pass(counter, sh["ref"], sh["chrom"], sh["rm"], sh["firewall"], sh["opts"], sh["hcands"][h0:h1], sh["pos"],
                             sh["m"], het_info, sh["record"], sh["d_min"])
    return first, xfirst, info, second, het_info


def _run_pool(workers: int, task, count: int):
    """``task(i)`` for i in range(count) on forked workers, results in order. A worker error, or a worker that dies (for example
    killed for memory), ends the pool and is raised; no result of a failed pool is used."""

    context = multiprocessing.get_context("fork")
    gc.freeze()
    pool = ProcessPoolExecutor(max_workers=min(workers, count), mp_context=context)
    try:
        results = list(pool.map(task, range(count), chunksize=1))
        pool.shutdown(wait=True)
        return results
    except BaseException:
        for process in list((getattr(pool, "_processes", None) or {}).values()):
            process.terminate()
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    finally:
        gc.unfreeze()


def _ranges(firsts: np.ndarray, ntiles: int, tile: int) -> list[tuple[int, int]]:
    """Index range of the candidates whose ``first`` lies in each tile."""

    if len(firsts) > 1 and (np.diff(firsts) < 0).any():
        raise RuntimeError("the candidate runs are not in order of their first position; the tiles cannot be cut")
    edges = np.searchsorted(firsts, np.arange(ntiles + 1, dtype=np.int64) * tile, side="left")
    edges[-1] = len(firsts)
    return [(int(edges[t]), int(edges[t + 1])) for t in range(ntiles)]


def scan_contig_tiled(make_counter, sequence, chrom, rm, firewall, opts: Options, words=None, rep_max=None, *, workers: int,
                      tile: int = TILE, s1_chunk: int = S1_CHUNK, timing: dict | None = None):
    """(first-pass places, haplotype-scan places or None, info) equal by value to ``scan.scan_contig_passes``, except that ``info``
    has no ``walk_log`` and ``cross_log``. ``make_counter(use)`` returns a counter (``use`` False for the first look, True for the
    walks). With one worker the call is ``scan_contig_passes`` itself."""

    if workers < 1:
        raise ValueError("workers must be at least 1")
    if workers == 1:
        return S.scan_contig_passes(make_counter(False), sequence, chrom, rm, firewall, opts, words, rep_max)
    ref = sequence.upper()
    het = opts.het_scan != "off"
    info: dict = {}
    n = max(0, (len(ref) - P.K) // P.STRIDE + 1)
    pos = np.arange(n, dtype=np.int64) * P.STRIDE
    bounds = [(ja, min(n, ja + s1_chunk)) for ja in range(0, n, s1_chunk)]
    _SHARED.clear()
    _SHARED.update(make_counter=make_counter, ref=ref, pos=pos, het=het, s1_bounds=bounds)
    started = time.time()
    try:
        if bounds:
            parts = _run_pool(workers, _s1_task, len(bounds))
            c = np.concatenate([p[0] for p in parts])
            x = np.concatenate([p[1] for p in parts]) if het else None
        else:
            c = np.full(0, -1, dtype=np.int64)
            x = np.full(0, -1, dtype=np.int64) if het else None
        depth, nfrac = P.local_depth(c)
        cands = P.candidates(pos, c, depth, nfrac)
        hcands, m = [], None
        if het:
            m = S.low_masks(c, nfrac, depth, x)
            hcands = S.candidates_from_low(pos, m["low"], depth)
            info.update(r0=m["r0"], n_old=int(m["old"].sum()), n_fwd=int(m["fwd"].sum()), n_bwd=int(m["bwd"].sum()),
                        n_low=int(m["low"].sum()))
        record = opts.check_scope != "none"
        d_min = S.depth_floor(depth, opts)
        if d_min is not None:
            info["depth_floor"] = round(d_min, 2)
        cross = dict(pos=pos, c=c, rep_max=P.REP_MAX if rep_max is None else rep_max) if opts.any_cross else None
        ntiles = len(ref) // tile + 1
        first_range = _ranges(np.array([cd["first"] for cd in cands], dtype=np.int64), ntiles, tile)
        het_range = _ranges(np.array([cd["first"] for cd in hcands], dtype=np.int64), ntiles, tile)
        busy = [t for t in range(ntiles) if first_range[t][1] > first_range[t][0] or het_range[t][1] > het_range[t][0]]
        if not busy:  # nothing to walk: one empty tile gives the stats of an empty scan
            busy, first_range[0], het_range[0] = [0], (0, 0), (0, 0)
        if timing is not None:
            timing["first_look"] = timing.get("first_look", 0.0) + time.time() - started
        started = time.time()
        _SHARED.update(first_range=[first_range[t] for t in busy], het_range=[het_range[t] for t in busy], cands=cands,
                       hcands=hcands, m=m, rm=rm, firewall=firewall, opts=opts, words=words, chrom=chrom, record=record,
                       d_min=d_min, cross=cross)
        results = _run_pool(workers, _tile_task, len(busy))
    finally:
        _SHARED.clear()
    if timing is not None:
        timing["walks"] = timing.get("walks", 0.0) + time.time() - started
    first: list = []
    xfirst: list = []
    second: list = []
    i_first = i_cross = i_het = None
    for f, xf, it, sec, ih in results:
        first += f
        if xf is not None:
            xfirst += xf
        i_first = it["first_pass"] if i_first is None else _join_stats(i_first, it["first_pass"])
        if "cross_rules" in it:
            i_cross = it["cross_rules"] if i_cross is None else _join_stats(i_cross, it["cross_rules"])
        if het:
            second += sec
            i_het = ih["het_pass"] if i_het is None else _join_stats(i_het, ih["het_pass"])
    info["first_pass"] = i_first
    if i_cross is not None:
        info["cross_rules"] = i_cross
    if cross is not None:
        info["first_merge"] = first
        first = xfirst
    if het:
        info["het_pass"] = i_het
    return first, second if het else None, info
