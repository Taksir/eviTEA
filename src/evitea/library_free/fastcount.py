"""An exact counter with the interface of ``_p4_lib.Counter`` that searches the msBWT less.

``Counter.cnt(w)`` runs two complete backward searches (the word and its reverse complement) and keeps its answer under the word as
walked. This counter gives the same numbers, by three changes that alter how often an answer is computed and never what it is:

* The interval of a word is the interval of the word without its first letter, extended by that letter (``findIndicesOfStr`` with
  ``givenRange``; a backward search is that same sequence of steps whether or not it is cut in two). The walk asks ``cnt(ctx + b)``
  for the four letters after one context and then ``cnt(ctx)`` for the next one, so most words cost one step.
* A word that is absent is found early: its last letters are searched first and the search stops when they are absent
  (``countOccurrencesOfSeq`` stops at an empty interval, ``findIndicesOfStr`` does not).
* The count cache is keyed on the strand-free word (the smaller of the word and its reverse complement), because the count is the
  sum over both strands.

During the first look (``use`` False) nothing is kept, as in ``Counter``. The caches empty themselves above a size cap; that changes
speed only. An index without ``findIndicesOfStr`` (the test fakes) gets plain searches with the strand-free cache.
"""

from __future__ import annotations

import random

from . import _p4_lib as P

__all__ = ["COUNT_CAP", "INTERVAL_CAP", "SPLIT", "FastCounter", "make_counter", "plain_count", "self_check"]

SPLIT = 14                 # letters searched first for an absent word
COUNT_CAP = 4_000_000      # counts kept before the count cache empties
INTERVAL_CAP = 8_000_000   # intervals kept before the interval cache empties
_COMP = P._C


def plain_count(bwt, w: str) -> int:
    """The count of ``_p4_lib.Counter`` for an uncached word: both strands, 0 for a word with a letter other than A, C, G, T."""

    if not P.is_acgt(w):
        return 0
    return int(bwt.countOccurrencesOfSeq(w)) + int(bwt.countOccurrencesOfSeq(P.rc(w)))


def _has_intervals(bwt) -> bool:
    """True when ``bwt`` continues a search from an interval (pymsbwt's ``findIndicesOfStr(seq, givenRange)``)."""

    find = getattr(bwt, "findIndicesOfStr", None)
    if find is None:
        return False
    try:
        low, high = find("A", (0, 0))
        return int(low) == int(high)
    except (TypeError, ValueError):
        return False


class FastCounter:
    """``Counter``-compatible attributes ``b``, ``n`` (index searches made), ``use``; methods ``cnt`` and ``clear``. It has no
    attribute ``cache``: the scan empties ``counter.cache`` by length, and this counter keeps its own limits."""

    def __init__(self, bwt, use: bool = False, *, count_cap: int = COUNT_CAP, interval_cap: int = INTERVAL_CAP, split: int = SPLIT):
        self.b = bwt
        self.n = 0
        self.use = use
        self.count_cap, self.interval_cap, self.split = count_cap, interval_cap, split
        self.intervals = _has_intervals(bwt)
        self._counts: dict[str, int] = {}
        self._iv: dict[str, tuple[int, int]] = {}

    def clear(self) -> None:
        self._counts.clear()
        self._iv.clear()

    def _search(self, s: str) -> tuple[int, int]:
        """Interval of ``s`` (a non-empty word), stopping early when its last letters are absent (the interval is then empty)."""

        find = self.b.findIndicesOfStr
        n = len(s)
        if n > self.split + 4:
            r = find(s[n - self.split:])
            self.n += 1
            if r[0] == r[1]:
                return r
            self.n += 1
            return find(s[:n - self.split], (r[0], r[1]))
        self.n += 1
        return find(s)

    def _interval(self, s: str) -> tuple[int, int]:
        iv = self._iv
        r = iv.get(s)
        if r is not None:
            return r
        if len(s) == 1:
            self.n += 1
            r = self.b.findIndicesOfStr(s)
        else:
            tail = s[1:]
            r1 = iv.get(tail)
            if r1 is None:
                r1 = self._search(tail)
                iv[tail] = r1
            if r1[0] == r1[1]:
                r = r1
            else:
                self.n += 1
                r = self.b.findIndicesOfStr(s[0], (r1[0], r1[1]))
        iv[s] = r
        return r

    def cnt(self, w: str) -> int:
        if not self.use:
            self.n += 2
            return plain_count(self.b, w)
        rw = w.translate(_COMP)[::-1]
        key = w if w <= rw else rw
        v = self._counts.get(key)
        if v is not None:
            return v
        if not P.is_acgt(w):
            return 0
        if self.intervals:
            a = self._interval(w)
            b = self._interval(rw)
            v = int(a[1]) - int(a[0]) + int(b[1]) - int(b[0])
            if len(self._iv) > self.interval_cap:
                self._iv.clear()
        else:
            self.n += 2
            v = plain_count(self.b, w)
        self._counts[key] = v
        if len(self._counts) > self.count_cap:
            self._counts.clear()
        return v


def make_counter(bwt, use: bool = False, **limits) -> FastCounter:
    return FastCounter(bwt, use, **limits)


def self_check(bwt, genome, contigs, n: int = 2000, seed: int = 0) -> dict:
    """Compare the fast counts with the plain counts on ``n`` words of the contigs to scan: 31-letter and 41-letter words at seeded
    positions, the same with one letter changed, and short walks (a word, the next word one letter on, and the four letters after a
    context). Any difference raises ``RuntimeError`` naming the word. Returns the number of words checked and whether the index
    continues searches from an interval."""

    fast = FastCounter(bwt, use=True)
    rng = random.Random(seed)
    names = [c for c in contigs if len(genome[c]) >= 2 * P.K + 20]
    checked = 0

    def check(word: str) -> None:
        nonlocal checked
        want = plain_count(bwt, word)
        got = fast.cnt(word)
        checked += 1
        if got != want:
            raise RuntimeError(f"fast counter self-check failed: {word} counts {got}, the plain search counts {want}")

    def change(word: str) -> str:
        i = rng.randrange(len(word))
        return word[:i] + rng.choice([x for x in "ACGT" if x != word[i]]) + word[i + 1:]

    while names and checked < n:
        seq = genome[rng.choice(names)].upper()
        p = rng.randrange(0, len(seq) - 2 * P.K - 10)
        w31, w41 = seq[p:p + P.K], seq[p:p + P.K + 10]
        for word in (w31, w41, change(w31), change(w41)):
            check(word)
        ctx = seq[p:p + P.K - 1]
        for i in range(p + P.K - 1, p + P.K + 3):
            rb = seq[i]
            check(ctx + rb)
            check(ctx)
            for b in "ACGT":
                check(ctx + b)
            ctx = ctx[1:] + rb
    return {"words": checked, "intervals": fast.intervals}
