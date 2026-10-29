"""Local heterozygous guard: the library-free calls that a region of homozygous insertions cannot hold.

An inbred or haploid sample carries its insertions in one allele, so a call whose reads still show the reference path (e above
HOM_MAX) among neighbours that all read as homozygous is most often noise: a diverged copy of a repeat, a read-through. A pool or
an outbred sample carries heterozygous insertions among the homozygous ones, and there the same call is common and real. The guard
asks the neighbours, so one rule fits both.

Confirmed places (``confirmed``): first-pass family calls with status ADDED or NEAR_EVITEA (judged after the resident-copy test and
before the single-copy rule and this guard, so neither drop changes the set), both junctions placed (bL and bR), a target-site
duplication of 2 to 20 bp and a measured e. A confirmed place is heterozygous-looking when e is above HOM_MAX.

A call that would be ADDED, in either pass, with e above HOM_MAX is DROPPED_HOM_REGION when fewer than BAR of the K confirmed places
nearest to it on the same contig (by distance of bL; ties to the left; the call itself counts when it is a confirmed place) are
heterozygous-looking. A contig with fewer than M confirmed places uses the share of the whole sample. A sample with fewer than NMIN
confirmed places has no guard (nothing is dropped). The guard acts after the single-copy rule, on added calls only.

The merge view takes the main view's confirmed places (``Guard.fixed``): both views judge a call against the same neighbours.

``--het-guard off`` switches it off. A recent private insertion in an inbred line is heterozygous, and the guard drops it when its
neighbours are homozygous: turn the guard off for a screen for new insertions in an inbred line.

The guard reads the family calls only (outside ``runtime.CODE_FILES``), so saved scan pieces are reused with it on or off.
"""

from __future__ import annotations

import bisect
import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from . import _p4_lib as P

__all__ = ["DROPPED", "MODES", "Guard", "Profile", "apply", "confirmed", "dropped_candidate"]

DROPPED = "DROPPED_HOM_REGION"
HOM_MAX = P.HOM_MAX        # the detector's bound of the homozygous genotype
K = 10                     # confirmed places nearest to a call
M = 5                      # a contig with fewer confirmed places uses the share of the whole sample
BAR = 0.25                 # drop when fewer than this share of the K nearest are heterozygous-looking
NMIN = 10                  # a sample with fewer confirmed places has no guard
TSD_MIN, TSD_MAX = 2, 20   # target-site duplication of a confirmed place
NEAR_BAR = 0.05            # report: dropped calls with e within this of HOM_MAX
LIST_MAX = 200             # dropped calls near the bar listed in the summary
MODES = ("on", "off")


def confirmed(scan: str, status: str, place: Mapping[str, str], e: float) -> bool:
    """Is this family call a confirmed place (see the module text)? ``status`` is the status before the single-copy rule."""

    if scan != "first" or status not in ("ADDED", "NEAR_EVITEA") or math.isnan(e):
        return False
    if place.get("bL", "") in ("", None) or place.get("bR", "") in ("", None):
        return False
    try:
        tsd = int(float(place["tsd"]))
    except (KeyError, TypeError, ValueError):
        return False
    return TSD_MIN <= tsd <= TSD_MAX


def dropped_candidate(e: float) -> bool:
    """Can the guard drop a call with this e? Only a heterozygous-looking one; no measured e is kept."""

    return not math.isnan(e) and e > HOM_MAX


@dataclass(frozen=True)
class Profile:
    """The confirmed places of a sample: per contig the (bL, heterozygous-looking) pairs in order of bL."""

    places: Mapping[str, tuple[tuple[int, bool], ...]]
    n: int
    het: int

    @classmethod
    def build(cls, items: Iterable[tuple[str, int, bool]]) -> "Profile":
        """From (contig, bL, heterozygous-looking) of every confirmed place."""

        by: dict[str, list[tuple[int, bool]]] = defaultdict(list)
        n = het = 0
        for chrom, pos, is_het in items:
            by[chrom].append((int(pos), bool(is_het)))
            n += 1
            het += bool(is_het)
        return cls({c: tuple(sorted(v)) for c, v in by.items()}, n, het)

    @property
    def sample_share(self) -> float:
        return self.het / self.n if self.n else 0.0

    @property
    def active(self) -> bool:
        return self.n >= NMIN

    def share(self, chrom: str, pos: int) -> float:
        """Share of heterozygous-looking places among the K confirmed places nearest to ``pos`` on ``chrom``; the share of the
        sample on a contig with fewer than M; 1.0 (never dropped) in a sample with fewer than NMIN."""

        if not self.active:
            return 1.0
        v = self.places.get(chrom, ())
        if len(v) < M:
            return self.sample_share
        keys = [p for p, _ in v]
        hi = bisect.bisect_left(keys, pos)
        lo = hi - 1
        het = taken = 0
        while taken < K and (lo >= 0 or hi < len(v)):
            if hi >= len(v) or (lo >= 0 and pos - v[lo][0] <= v[hi][0] - pos):
                het += v[lo][1]
                lo -= 1
            else:
                het += v[hi][1]
                hi += 1
            taken += 1
        return het / taken

    def drop(self, chrom: str, pos: int, e: float) -> bool:
        return dropped_candidate(e) and self.share(chrom, pos) < BAR

    def summary(self) -> dict:
        return dict(confirmed=self.n, heterozygous_looking=self.het, sample_share=round(self.sample_share, 4),
                    contigs=len(self.places))


class Guard:
    """``--het-guard``: ``on``/``off`` and, for the merge view, the main view's profile (``fixed``). After ``apply`` the profile that
    was built from the view's own calls is in ``built``."""

    def __init__(self, on: bool = True, fixed: Profile | None = None, source: str = ""):
        self.on = on
        self.fixed = fixed
        self.source = source
        self.built: Profile | None = None

    @classmethod
    def parse(cls, text: str) -> "Guard":
        if text not in MODES:
            raise ValueError(f"--het-guard: give on or off, got {text!r}")
        return cls(text == "on")

    def for_merge_view(self) -> "Guard":
        """The guard of the merge view: the main view's confirmed places."""

        return Guard(self.on, self.built if self.built is not None else self.fixed, "main view")


def apply(guard: Guard, items: Iterable[tuple[str, int, bool]], passes: Sequence[tuple]) -> dict:
    """Drop the added calls the guard rejects and return the summary. ``items``: (contig, bL, heterozygous-looking) of the confirmed
    places of this view, taken before the single-copy rule; ``passes``: (pass, status dict) with ``pass.scan``, ``pass.places``,
    ``pass.calls`` and ``pass.e``; the status of a call that would be ADDED becomes DROPPED_HOM_REGION."""

    if not guard.on:
        return dict(mode="off")
    profile = guard.fixed if guard.fixed is not None else Profile.build(items)
    if guard.fixed is None:
        guard.built = profile
    out: dict = dict(mode="on", source=guard.source or "own calls", K=K, M=M, bar=BAR, nmin=NMIN, **profile.summary())
    if not profile.active:
        out.update(active=False, reason=f"fewer than {NMIN} confirmed places: no guard", dropped={"first": 0, "second": 0})
        return out
    drops = {"first": 0, "second": 0}
    per_contig: dict[str, int] = defaultdict(int)
    near = []
    for pa, st in passes:
        for i in pa.calls:
            if st[i] != "ADDED":
                continue
            p = pa.places[i]
            e = pa.e[i]
            if profile.drop(p["chrom"], int(p["bL"]), e):
                st[i] = DROPPED
                drops[pa.scan] += 1
                per_contig[p["chrom"]] += 1
                if e <= HOM_MAX + NEAR_BAR:
                    near.append(dict(scan=pa.scan, chrom=p["chrom"], bL=int(p["bL"]), e=round(e, 4)))
    out.update(active=True, dropped=drops, dropped_per_contig=dict(sorted(per_contig.items())),
               near_bar=dict(band=NEAR_BAR, dropped=len(near), calls=near[:LIST_MAX]))
    return out
