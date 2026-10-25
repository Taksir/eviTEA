"""Record building: consensus of the member segments of a family, reference copies from exact words, end-word codebooks.

Record      the member segments of a family are oriented against a seed, piled on the seed per end, and voted per column.
            The seed is the member with the most letters whose two segments both hold at least ORIENT_K letters. Up to
            SEED_CANDIDATES seeds are tried, and the one that places most members wins. A letter wins a column when it has more
            than half of the covering members and at least FAMILY_MIN members, otherwise the column is N. The two end consensuses
            join where a suffix of the first is exactly a prefix of the second (JOIN_K letters at least, because the walks met
            inside a short element). Otherwise they are separated by SPACER N.
Copies      a reference copy of a record is a chain of exact WORD_K-letter words of the record on one contig and strand. Each
            word lies at most COPY_GAP letters after the chain's last word, and its offset (reference position minus record
            position) is within DIAG_SLACK of the chain's. A chain is a copy when it holds two non-overlapping words.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

FAMILY_MIN, ORIENT_K, JOIN_K, SPACER, WORD_K, SEED_CANDIDATES = 2, 15, 31, 100, 21, 20
COPY_GAP, DIAG_SLACK = 100, 30
PAD_TO = 100

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass(frozen=True)
class Record:
    seq: str
    members: int          # usable members (both segments together hold >= ORIENT_K letters)
    placed: int           # members placed against the seed
    a_len: int            # letters of the first-end consensus
    b_len: int            # letters of the second-end consensus
    joined: str           # "overlap" | "spacer" | "one end"


def _mappy():
    import mappy

    return mappy


def _best(al, q):
    if al is None or len(q) < ORIENT_K:
        return None
    hits = [h for h in al.map(q) if h.strand == 1]
    return max(hits, key=lambda h: h.mlen) if hits else None


def _pile(h, q, counts, depth, member):
    qi, ri = h.q_st, h.r_st
    for n, op in h.cigar:
        if op in (0, 7, 8):
            for j in range(n):
                counts[ri + j][q[qi + j]] += 1
                depth[ri + j].add(member)
            qi += n
            ri += n
        elif op == 1:
            qi += n
        elif op == 2:
            for j in range(n):
                counts[ri + j]["-"] += 1
                depth[ri + j].add(member)
            ri += n


def _vote(counts, depth):
    letters = []
    for c, d in zip(counts, depth):
        cover = sum(c.values())
        if len(d) < FAMILY_MIN or not cover:
            letters.append("")
            continue
        letter, n = c.most_common(1)[0]
        letters.append(letter if 2 * n > cover else "N")
    ok = [i for i, x in enumerate(letters) if x and x not in "N-"]
    if not ok:
        return ""
    return "".join(x if x else "N" for x in letters[ok[0]:ok[-1] + 1]).replace("-", "")


def _place(alA, alB, lseg, rseg):
    L, R = lseg.upper(), rseg.upper()
    opts = []
    for qa, qb in ((L, R), (rc(R), rc(L))):
        ha, hb = _best(alA, qa), _best(alB, qb)
        opts.append(((ha.mlen if ha else 0) + (hb.mlen if hb else 0), qa, ha, qb, hb))
    opts.sort(key=lambda x: -x[0])
    return None if opts[0][0] == 0 or opts[0][0] == opts[1][0] else opts[0]


def join_ends(a: str, b: str) -> tuple[str, str]:
    """Merge the two ends only when a suffix of the first is exactly a prefix of the second. A shared word alone is not enough,
    because terminal repeats put the same words at both ends."""
    for p in range(max(0, len(a) - len(b)), len(a) - JOIN_K + 1):
        if b.startswith(a[p:]):
            return a[:p] + b, "overlap"
    return a + "N" * SPACER + b, "spacer"


def build_record(members) -> Record | None:
    """Build a record from ``members``, the (left segment, right segment) pairs of one family. Return None when fewer than
    FAMILY_MIN members are usable or no column reaches the vote."""
    mp = _mappy()
    usable = [m for m in members if len(m[0]) + len(m[1]) >= ORIENT_K]
    if len(usable) < FAMILY_MIN:
        return None
    both = [m for m in usable if len(m[0]) >= ORIENT_K and len(m[1]) >= ORIENT_K]
    cands = sorted(both or usable, key=lambda m: -(len(m[0]) + len(m[1])))[:SEED_CANDIDATES]
    best = None
    for seed in cands:
        A, B = seed[0].upper(), seed[1].upper()
        alA = mp.Aligner(seq=A, preset="map-ont") if len(A) >= ORIENT_K else None
        alB = mp.Aligner(seq=B, preset="map-ont") if len(B) >= ORIENT_K else None
        placed = [(i, x) for i, x in ((i, _place(alA, alB, *m)) for i, m in enumerate(usable)) if x]
        if best is None or len(placed) > len(best[2]):
            best = (A, B, placed)
    A, B, placed = best
    cA, dA = [Counter() for _ in A], [set() for _ in A]
    cB, dB = [Counter() for _ in B], [set() for _ in B]
    for i, (_, qa, ha, qb, hb) in placed:
        if ha:
            _pile(ha, qa, cA, dA, i)
        if hb:
            _pile(hb, qb, cB, dB, i)
    a, b = _vote(cA, dA), _vote(cB, dB)
    if not a and not b:
        return None
    if a and b:
        seq, how = join_ends(a, b)
    else:
        seq, how = a or b, "one end"
    return Record(seq, len(usable), len(placed), len(a), len(b), how)


def _reference_chunks(path):
    """Yield (contig name, uppercase sequence) of a FASTA."""
    from ..languages.common import parse_fasta

    for header, seq in parse_fasta(path):
        yield header.split()[0], seq


def word_copies(records: dict[str, str], reference_fasta) -> list[tuple[str, int, int, str]]:
    """Return the reference copies of each record from its own WORD_K-letter words, as rows (contig, begin, end, record id)."""
    index: dict[str, list] = {}
    for te, seq in records.items():
        for i in range(len(seq) - WORD_K + 1):
            w = seq[i:i + WORD_K]
            if "N" not in w:
                index.setdefault(w, []).append((te, i, 1))
                index.setdefault(rc(w), []).append((te, i, -1))
    rows = []
    for name, seq in _reference_chunks(reference_fasta):
        hits: dict = {}
        for p in range(len(seq) - WORD_K + 1):
            for te, i, s in index.get(seq[p:p + WORD_K], ()):
                hits.setdefault((te, s), []).append((p, p - i if s == 1 else p + i, i))
        for (te, s), hs in hits.items():
            hs.sort()
            active: dict = {}
            chains = []
            for n, (p, d, i) in enumerate(hs):
                b = d // DIAG_SLACK
                ch = next((c for bb in (b - 1, b, b + 1) for c in active.get(bb, ())
                           if p - c[0] <= COPY_GAP and abs(d - c[1]) <= DIAG_SLACK), None)
                if ch is None:
                    ch = [p, d, i, i, p]
                    active.setdefault(b, []).append(ch)
                    chains.append(ch)
                elif ch[1] // DIAG_SLACK != b:
                    active[ch[1] // DIAG_SLACK].remove(ch)
                    active.setdefault(b, []).append(ch)
                ch[0], ch[1], ch[2], ch[3] = p, d, min(ch[2], i), max(ch[3], i)
                if n % 256 == 0:
                    for bb in list(active):
                        active[bb] = [c for c in active[bb] if p - c[0] <= COPY_GAP]
            rows += [(name, c[4], c[0] + WORD_K, te) for c in chains if c[3] - c[2] >= WORD_K]
    return rows


def reference_absent_share(seq: str, ref_count: dict[str, int]) -> float:
    """Return the fraction of the record's WORD_K-letter words (both strands as one) that occur nowhere in the reference."""
    ws = [w for w in (seq[i:i + WORD_K] for i in range(len(seq) - WORD_K + 1)) if "N" not in w]
    if not ws:
        return 0.0
    canon = [min(w, rc(w)) for w in ws]
    return round(sum(ref_count.get(c, 0) == 0 for c in canon) / len(canon), 4)


def count_record_words(records: dict[str, str], reference_fasta) -> dict[str, int]:
    """Map each word of the records (both strands as one) to its occurrences in the reference."""
    look = {}
    for seq in records.values():
        for i in range(len(seq) - WORD_K + 1):
            w = seq[i:i + WORD_K]
            if "N" not in w:
                c = min(w, rc(w))
                look[w] = c
                look[rc(w)] = c
    cnt: Counter = Counter()
    for _, seq in _reference_chunks(reference_fasta):
        for i in range(len(seq) - WORD_K + 1):
            c = look.get(seq[i:i + WORD_K])
            if c is not None:
                cnt[c] += 1
    return dict(cnt)


def languages(te_fa, fmap, bed, ref):
    """Return the end-word codebooks (terminal, internal-only) of the library builder. When its gap fill stops because no record
    is short of end words, the rebuilt words pass through."""
    from ..languages import internal as I
    from ..languages import terminal as T
    from ..languages.builder import build_languages, g1a_tables

    try:
        b = build_languages(te_fasta=te_fa, family_map=fmap, reference_fasta=ref, repeatmasker_bed=bed)
        return b.terminal_codebook, b.internal_only_codebook
    except RuntimeError as e:
        if "No underpowered records" not in str(e):
            raise
    cat, prov = g1a_tables(te_fa, fmap)
    assets = dict(te_fasta=te_fa, family_map=fmap, reference_fasta=ref, repeatmasker_bed=bed)
    rb = T.rebuild_terminal21(catalogue=cat, provenance=prov, **assets)
    hyb = []
    for r in T.read_tsv_bytes(rb.codebook):
        r = dict(r)
        r["k"] = 21
        r["library_role"] = r.get("library_role", "REBUILT_SAFE_TERMINAL21_ACQUISITION")
        hyb.append(r)
    fin = T.finalize_terminal(hybrid_codebook=T.tsv_bytes(hyb, T.CODEBOOK_FIELDS), rebuilt_source_audit=rb.source_audit,
                              gapfill_source_audit=b"", rebuilt_te_summary=rb.te_record_summary)
    it = I.build_internal21(catalogue=cat, provenance=prov, terminal_codebook=fin.codebook, **assets)
    return fin.codebook, it.internal_only_codebook


def sha12(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:12]
