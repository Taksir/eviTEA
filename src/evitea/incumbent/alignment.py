"""Pure-Python infix edit-distance kernel with edlib ``HW`` semantics.

It reproduces two results of ``edlib.align(query, target, mode="HW",
task=..., k=max_edits)`` without a compiled dependency:

* ``hw_distance``: the minimum Levenshtein distance between the complete query
  and any substring of the target, or ``None`` when it exceeds ``max_edits``;
* ``hw_locations``: that distance plus every ``(start, end)`` pair, where
  ``end`` ranges over all target positions attaining the minimum (ascending)
  and ``start`` is the *smallest* start whose substring attains the same
  distance for that end. This is edlib's rule: it reruns a prefix-anchored
  (``SHW``) search on the reversed query and target and takes the last
  reported position.

Characters compare by exact identity (``N`` is an ordinary symbol), as in
edlib's default equality. The implementation is Myers' bit-parallel algorithm
on Python integers. A pigeonhole exact-piece filter only proves that no
alignment within ``max_edits`` exists, so it cannot change any result.
"""

from __future__ import annotations

from collections.abc import Sequence

__all__ = ["hw_distance", "hw_locations", "hw_column_scores"]


def _peq(query: str) -> dict[str, int]:
    table: dict[str, int] = {}
    for index, symbol in enumerate(query):
        table[symbol] = table.get(symbol, 0) | (1 << index)
    return table


def _pigeonhole_rejects(query: str, target: str, max_edits: int) -> bool:
    """Return True only when no substring can be within ``max_edits``."""

    pieces = max_edits + 1
    length = len(query)
    if length < pieces:
        return False
    base, extra = divmod(length, pieces)
    cursor = 0
    for index in range(pieces):
        size = base + (1 if index < extra else 0)
        if query[cursor : cursor + size] in target:
            return False
        cursor += size
    return True


def hw_column_scores(query: str, target: str, *, prefix_anchored: bool = False) -> list[int]:
    """Return the bottom-row score after every target character.

    ``prefix_anchored=False`` is infix (HW) mode.  ``True`` is SHW mode: the
    alignment must start at target position 0 but may end anywhere.
    """

    m = len(query)
    peq = _peq(query)
    mask = (1 << m) - 1
    high = 1 << (m - 1)
    pv = mask
    mv = 0
    score = m
    hin = 1 if prefix_anchored else 0
    scores: list[int] = []
    append = scores.append
    get = peq.get
    for symbol in target:
        eq = get(symbol, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & mask)
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph = ((ph << 1) | hin) & mask
        mh = (mh << 1) & mask
        pv = (mh | ~(xv | ph)) & mask
        mv = ph & xv
        append(score)
    return scores


def _valid(query: str, target: str, max_edits: int) -> bool:
    return bool(query) and bool(target) and max_edits >= 0


def _candidate_windows(query: str, target: str, max_edits: int) -> list[tuple[int, int]] | None:
    """Return merged target windows that contain every alignment within ``max_edits``.

    Any alignment with at most ``max_edits`` edits leaves one of the
    ``max_edits + 1`` disjoint query pieces unedited, so it contains an exact
    piece occurrence. Return ``None`` when the query is too short to split.
    """

    pieces = max_edits + 1
    m = len(query)
    if m < pieces:
        return None
    base, extra = divmod(m, pieces)
    n = len(target)
    windows: list[tuple[int, int]] = []
    offset = 0
    for index in range(pieces):
        size = base + (1 if index < extra else 0)
        piece = query[offset : offset + size]
        position = target.find(piece)
        while position >= 0:
            windows.append(
                (max(0, position - offset - max_edits), min(n, position - offset + m + max_edits))
            )
            position = target.find(piece, position + 1)
        offset += size
    if not windows:
        return []
    windows.sort()
    merged = [windows[0]]
    for start, end in windows[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            if end > last_end:
                merged[-1] = (last_start, end)
        else:
            merged.append((start, end))
    return merged


def hw_distance(query: str, target: str, max_edits: int) -> int | None:
    """Return the infix edit distance of the complete query, or ``None`` if it exceeds ``max_edits``."""

    if not _valid(query, target, max_edits):
        return None
    m = len(query)
    windows = _candidate_windows(query, target, max_edits)
    if windows is None:
        best = min(m, min(hw_column_scores(query, target)))
    elif not windows:
        return None
    else:
        best = m
        for start, end in windows:
            value = min(hw_column_scores(query, target[start:end]))
            if value < best:
                best = value
    return best if best <= max_edits else None


def hw_locations(
    query: str, target: str, max_edits: int
) -> tuple[int | None, tuple[tuple[int, int], ...]]:
    """Return the distance and all edlib-equivalent ``(start, end)`` locations (0-based, inclusive)."""

    if not _valid(query, target, max_edits):
        return None, ()
    if _pigeonhole_rejects(query, target, max_edits):
        return None, ()
    scores = hw_column_scores(query, target)
    m = len(query)
    best = min(m, min(scores))
    if best > max_edits:
        return None, ()
    reversed_query = query[::-1]
    locations: list[tuple[int, int]] = []
    if best == m:
        # edlib reports the empty alignment before the target as end = -1.
        locations.append((0, -1))
    for end, value in enumerate(scores):
        if value != best:
            continue
        # Prefix scores beyond m + best columns necessarily exceed best.
        span_limit = min(end + 1, m + best + 1)
        reversed_target = target[end::-1][:span_limit]
        prefix_scores = hw_column_scores(
            reversed_query, reversed_target, prefix_anchored=True
        )
        last = max(index for index, score in enumerate(prefix_scores) if score == best)
        locations.append((end - last, end))
    return best, tuple(locations)


def best_symmetric(scores: Sequence[int | None]) -> int | None:
    """Return the minimum over non-missing distances, or ``None`` when all are missing."""

    values = [value for value in scores if value is not None]
    return min(values) if values else None
