"""Single-copy contigs: the library-free calls a contig carried in one copy cannot hold.

A contig that the sample carries in one copy (the X and Y of a male, a monosomic chromosome, every contig of a haploid genome)
has one haplotype, so a real insertion there is carried by every read at the break (reference-path fraction e near 0, the
homozygous genotype). The haplotype scan looks for junctions that only part of the reads carry, and the first pass keeps calls
with the heterozygous genotype; on such a contig both give only noise. After every other status, a call that would be ADDED on
a single-copy contig is DROPPED_SINGLE_COPY when

S1  it is a haplotype-scan (second-scan) call, or
S2  it is a first-pass call with e above HOM_MAX (the detector's bound of the homozygous genotype).

Grouping, families, the tolerance and the other statuses do not change, and no call outside the added calls changes.

Which contigs (``--single-copy-contigs``):

auto   measured. A contig is single copy when its depth Dc is at most DEPTH_RATIO_MAX of the genome depth Dg (one copy gives
       0.5, two give 1, a pool of equal numbers of males and females 0.75 on the X) and its library-confirmed first-pass calls
       (NEAR_EVITEA) are mostly homozygous: at least GUARD_MIN_CALLS such calls with a measured e, at most GUARD_HET_MAX of them
       above HOM_MAX. Depth cannot tell one male from a pool of males (both 0.5 on the X), and a pool carries heterozygous
       insertions, which the guard sees. With fewer library-confirmed calls on the contig, auto leaves it alone.
       Depth: about N_WORDS words of K letters at an even stride on the measured contigs (the default contigs of the
       reference and the run's contigs), each occurring once in the reference (both strands), counted in the sample msBWT
       (both strands). Dc = median count of the contig's words that the sample carries (count above zero), Dg = the same over
       all measured contigs; a contig with fewer than MIN_WORDS such words is not judged.
all    every contig (a haploid genome: depth cannot show haploidy, since every contig then has the same depth).
none   off.
LIST   the named contigs, comma-separated (for example chrX,chrY for a male).

The measure is outside the detection code (``runtime.CODE_FILES``), so saved scan pieces are reused with any choice; it is kept
in the work directory (``copy_profile.json``) with the identity of the reference and the msBWT.
"""

from __future__ import annotations

import math
import statistics
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from . import _p4_lib as P

__all__ = ["DROPPED", "Choice", "Decision", "dropped", "measure", "report"]

DROPPED = "DROPPED_SINGLE_COPY"
K = P.K                    # word length of the depth measure (as the scan's words)
N_WORDS = 200_000          # sampled words over the measured contigs
MIN_WORDS = 1_000          # a contig with fewer counted words is not judged
DEPTH_RATIO_MAX = 0.625    # midway between one copy (0.5) and a pool of equal numbers of males and females (0.75)
GUARD_MIN_CALLS = 10       # library-confirmed first-pass calls on the contig
GUARD_HET_MAX = 0.25       # at most this share of them with e above HOM_MAX
NEAR_BAR = 0.05            # report: calls with e within this of HOM_MAX
CHUNK = 8_000_000          # bases per step of the reference pass
PREFIX = 12                # letters of the prefix filter of the reference pass
MODES = ("auto", "all", "none")


@dataclass(frozen=True)
class Choice:
    """``--single-copy-contigs``: mode auto, all, none or list (with the names)."""

    mode: str = "auto"
    names: tuple[str, ...] = ()

    @classmethod
    def parse(cls, text: str) -> "Choice":
        text = (text or "").strip()
        if text in MODES:
            return cls(text)
        names = tuple(x.strip() for x in text.split(",") if x.strip())
        if not names:
            raise ValueError("--single-copy-contigs: give auto, all, none or comma-separated contig names")
        if len(set(names)) != len(names):
            raise ValueError(f"--single-copy-contigs: a contig is named twice: {text}")
        return cls("list", names)

    def as_text(self) -> str:
        return ",".join(self.names) if self.mode == "list" else self.mode


def dropped(scan: str, e: float, hom_max: float = P.HOM_MAX) -> bool:
    """S1 / S2 for a call that would be added on a single-copy contig. A first-pass call with no measured e is kept."""

    return scan == "second" or (not math.isnan(e) and e > hom_max)


@dataclass(frozen=True)
class Decision:
    """The choice and, for auto, the depth profile (``measure``) or the reason there is none. ``fixed``: the contigs another view
    decided on (the merge view takes the main view's decision, since it describes the sample, not a rule set)."""

    choice: Choice
    profile: dict | None = None
    reason: str = ""
    fixed: tuple[str, ...] | None = None
    source: str = ""

    @classmethod
    def from_view(cls, choice: Choice, flagged: Iterable[str], source: str) -> "Decision":
        return cls(choice, fixed=tuple(sorted(flagged)), source=source)

    def decide(self, near_calls: Iterable[tuple[str, float]], contigs: Iterable[str],
               hom_max: float = P.HOM_MAX) -> tuple[frozenset, dict]:
        """(flagged contigs, summary). ``near_calls``: (contig, e) of the first-pass NEAR_EVITEA calls; ``contigs``: the contigs
        of the family calls."""

        out: dict = {"mode": self.choice.as_text()}
        if self.fixed is not None:
            out.update(flagged=list(self.fixed), source=self.source)
            return frozenset(self.fixed), out
        if self.choice.mode == "none":
            return frozenset(), out
        if self.choice.mode == "all":
            flagged = frozenset(contigs)
            out["flagged"] = sorted(flagged)
            return flagged, out
        if self.choice.mode == "list":
            out["flagged"] = list(self.choice.names)
            return frozenset(self.choice.names), out
        guard = defaultdict(lambda: [0, 0])
        for chrom, e in near_calls:
            if not math.isnan(e):
                guard[chrom][0] += 1
                guard[chrom][1] += e > hom_max
        if self.profile is None:
            out.update(flagged=[], reason=self.reason or "no depth profile")
            return frozenset(), out
        flagged, per = set(), {}
        for chrom, m in self.profile["contigs"].items():
            n, het = guard.get(chrom, (0, 0))
            low = m["ratio"] is not None and m["words"] >= MIN_WORDS and m["ratio"] <= DEPTH_RATIO_MAX
            hom = n >= GUARD_MIN_CALLS and het <= GUARD_HET_MAX * n
            per[chrom] = dict(depth=m["depth"], words=m["words"], ratio=m["ratio"], near_evitea=n, near_evitea_het=het,
                              low_depth=low, guard=hom, single_copy=low and hom)
            if low and hom:
                flagged.add(chrom)
        out.update(genome_depth=self.profile["genome_depth"], words=self.profile["words"], contigs=per, flagged=sorted(flagged))
        return frozenset(flagged), out


def _codes(words: Sequence[str]) -> np.ndarray:
    """2-bit codes (A0 C1 G2 T3) of ACGT words of K letters."""

    if not words:
        return np.zeros(0, dtype=np.uint64)
    a = P._ENC[np.frombuffer("".join(words).encode(), dtype=np.uint8)].reshape(-1, K).astype(np.uint64)
    out = np.zeros(len(words), dtype=np.uint64)
    for j in range(K):
        out = (out << np.uint64(2)) | a[:, j]
    return out


def _codes16(v: np.ndarray) -> np.ndarray:
    """Codes of the 16-letter words of the letter codes ``v`` (0 to 3), by doubling: ``len(v) - 15`` of them."""

    f = v.astype(np.uint32)
    for h in (1, 2, 4, 8):
        m = len(f) - h
        f = (f[:m] << np.uint32(2 * h)) | f[h:h + m]
    return f


def _sample(genome: Mapping[str, str], contigs: Sequence[str], n_words: int):
    """(contig of each word, words) at an even stride over ``contigs``; words with a non-ACGT letter are left out."""

    total = sum(len(genome[c]) for c in contigs)
    stride = max(1, total // max(1, n_words))
    owner, words = [], []
    for chrom in contigs:
        seq = genome[chrom]
        for x in range(stride // 2, len(seq) - K + 1, stride):
            w = seq[x:x + K].upper()
            if P.is_acgt(w):
                owner.append(chrom)
                words.append(w)
    return owner, words


def _reference_copies(genome: Mapping[str, str], words: Sequence[str]) -> np.ndarray:
    """Copies of each word in the reference, both strands (forward occurrences of the word and of its reverse complement), in
    one pass over every contig. A prefix table of the wanted words picks the positions whose full word is looked up."""

    if K != 31:
        raise AssertionError("the reference pass joins two 16-letter codes into a 31-letter one")
    fwd, rev = _codes(words), _codes([P.rc(w) for w in words])
    q = np.unique(np.concatenate([fwd, rev]))
    occ = np.zeros(len(q), dtype=np.int64)
    if not len(q):
        return np.zeros(0, dtype=np.int64)
    table = np.zeros(1 << (2 * PREFIX), dtype=bool)
    table[(q >> np.uint64(2 * (K - PREFIX))).astype(np.int64)] = True
    low = np.uint32(2 * (16 - PREFIX))
    for seq in genome.values():
        for s in range(0, max(0, len(seq) - K + 1), CHUNK):
            a = P._ENC[np.frombuffer(seq[s:s + CHUNK + K - 1].upper().encode(), dtype=np.uint8)]
            n = len(a) - K + 1
            bad = np.concatenate([[0], np.cumsum(a == 4)])
            valid = (bad[K:K + n] - bad[:n]) == 0
            f16 = _codes16(np.where(a == 4, 0, a))
            at = np.flatnonzero(table[f16[:n] >> low] & valid)
            if not len(at):
                continue
            c = ((f16[at] >> np.uint32(2)).astype(np.uint64) << np.uint64(32)) | f16[at + 15].astype(np.uint64)
            j = np.searchsorted(q, c)
            j[j == len(q)] = 0
            hit = q[j] == c
            occ += np.bincount(j[hit], minlength=len(q))
    return occ[np.searchsorted(q, fwd)] + occ[np.searchsorted(q, rev)]


def measure(genome: Mapping[str, str], count: Callable[[str], int], contigs: Sequence[str], *, n_words: int = N_WORDS,
            progress: Callable[[str], None] | None = None) -> dict:
    """Depth profile of ``contigs`` (``count``: sample count of a word, both strands; ``_p4_lib.Counter.cnt``)."""

    say = progress or (lambda message: None)
    owner, words = _sample(genome, contigs, n_words)
    say(f"library-free: single-copy contigs: {len(words)} sampled words; copies in the reference")
    copies = _reference_copies(genome, words)
    say(f"library-free: single-copy contigs: {int((copies == 1).sum())} words once in the reference; counting them in the sample")
    by = defaultdict(list)
    counted = 0
    for chrom, w, n in zip(owner, words, copies):
        if n != 1:
            continue
        counted += 1
        c = int(count(w))
        if c > 0:
            by[chrom].append(c)
    kept = [c for v in by.values() for c in v]
    genome_depth = float(statistics.median(kept)) if kept else None
    out = {}
    for chrom in contigs:
        v = by.get(chrom, [])
        depth = float(statistics.median(v)) if v else None
        ratio = None if depth is None or not genome_depth else round(depth / genome_depth, 4)
        out[chrom] = dict(depth=depth, words=len(v), ratio=ratio)
    return dict(genome_depth=genome_depth, words=dict(sampled=len(words), unique=counted, counted=len(kept)), contigs=out)


def report(calls: Iterable[tuple[str, str, str, str, float]], flagged: frozenset, hom_max: float = P.HOM_MAX) -> dict:
    """Counts after the rule, from (scan, cluster, contig, status, e) of every family call: calls dropped per scan; families
    (scan, cluster) with more than half of their calls dropped, and their added calls on other contigs (a family the rule
    takes apart elsewhere); first-pass calls on single-copy contigs with e within NEAR_BAR of HOM_MAX, added and dropped."""

    calls = list(calls)
    drops = {"first": 0, "second": 0}
    fam = defaultdict(lambda: [0, 0])
    for scan, cl, _chrom, st, _e in calls:
        fam[(scan, cl)][0] += 1
        if st == DROPPED:
            fam[(scan, cl)][1] += 1
            drops[scan] += 1
    mostly = {k for k, (n, d) in fam.items() if d * 2 > n}
    elsewhere = sum(1 for scan, cl, chrom, st, _e in calls if st == "ADDED" and (scan, cl) in mostly and chrom not in flagged)
    kept = sum(1 for scan, _cl, chrom, st, e in calls if scan == "first" and chrom in flagged and st == "ADDED"
               and not math.isnan(e) and hom_max - NEAR_BAR < e <= hom_max)
    cut = sum(1 for scan, _cl, chrom, st, e in calls if scan == "first" and chrom in flagged and st == DROPPED
              and e <= hom_max + NEAR_BAR)
    return dict(dropped=drops, families_mostly_dropped=len(mostly), added_elsewhere_in_them=elsewhere,
                near_bar=dict(added=kept, dropped=cut, band=NEAR_BAR))
