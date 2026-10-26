"""Run per-locus steps in deterministic locus batches.

Split-junction and mate-anchored search evidence is computed independently per
candidate locus.  Contiguous locus batches are evaluated (optionally in worker
processes) and their results concatenated in locus order, so the downstream
competitions see exactly the sequential input order.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from typing import Any

from .mate_rescue import mate_rescue_evidence
from .split_junction import split_junction_hypotheses

__all__ = ["map_locus_batches"]

_STATE: dict[str, Any] = {}
_LAYERS: dict[str, Callable[..., tuple[list, dict[str, int]]]] = {
    "split_junction": split_junction_hypotheses,
    "mate_rescue": mate_rescue_evidence,
}


def _initialize(state: Mapping[str, Any]) -> None:
    _STATE.clear()
    _STATE.update(state)


def _run_batch(layer: str, loci: Sequence) -> tuple[list, dict[str, int]]:
    return _LAYERS[layer](loci=loci, **_STATE)


def map_locus_batches(
    layer: str,
    loci: Sequence,
    state: Mapping[str, Any],
    *,
    workers: int = 1,
    batch_size: int = 25,
    progress: Callable[[int, int], None] | None = None,
) -> tuple[list, dict[str, int]]:
    """Run ``layer`` over ``loci`` and return (ordered rows, merged counts)."""

    batches = [list(loci[i : i + batch_size]) for i in range(0, len(loci), batch_size)]
    rows: list = []
    counts: Counter[str] = Counter()
    if workers <= 1 or len(batches) <= 1:
        _initialize(state)
        results = (_run_batch(layer, batch) for batch in batches)
        for done, (batch_rows, batch_counts) in enumerate(results, 1):
            rows.extend(batch_rows)
            counts.update(batch_counts)
            if progress is not None:
                progress(done, len(batches))
        _STATE.clear()
        return rows, dict(counts)
    with ProcessPoolExecutor(max_workers=workers, initializer=_initialize, initargs=(dict(state),)) as pool:
        futures = [pool.submit(_run_batch, layer, batch) for batch in batches]
        for done, future in enumerate(futures, 1):
            batch_rows, batch_counts = future.result()
            rows.extend(batch_rows)
            counts.update(batch_counts)
            if progress is not None:
                progress(done, len(batches))
    return rows, dict(counts)
