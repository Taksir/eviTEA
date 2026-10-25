"""P4 library-free detector, stages S1-S7 (p4_library_free/P4_DESIGN_MEMO.md section 4). Step 2b-i.
No TE library, Dfam sequence, long-read insert, truth or hostile coordinate is read anywhere in this file.

Counter interface (as g4b oracle_r2_walk.Counter): C.cnt(w) = occurrences of w + occurrences of rc(w) in the reads,
0 when w contains a letter other than A/C/G/T.

Coordinates: 0-based, half-open. A sample allele ref[:bL] + INSERT + ref[bR:] has left break bL (first reference letter the
left walk cannot follow) and right break bR (first reference letter the sample shares again, seen from the right).
delta = bL - bR; a target-site duplication of t letters gives delta = t (TSD = ref[bR:bL]).
"""
import warnings
import numpy as np

from .. import depth_window as _dw

# ---- parameters (memo section 4; "fixed" unless marked synthetic) ----
K = 31                  # word length (S1, walk)
STRIDE = 10             # S1 sampling stride (fixed; stride <= 15 hits a homozygous insertion >= 2 times)
LOWF = 0.65             # S1 drop threshold, x local depth (fixed; oracle LOW constant)
DFLANK = 2000           # S1 local depth = median of sampled counts within +-2 kb (fixed)
NFRAC = 0.20            # S1 skip windows with > 20% non-ACGT words (fixed)
MIN_LOW = 1             # S1 a candidate needs >= 1 low sampled word: an insertion with a t-letter TSD removes only 30 - t
                        # reference words, so stride 10 sees >= 1 of them for t <= 20 and >= 2 only for t <= 10 (F-031)
RUN_GAP = 2 * STRIDE    # S1 low words <= 20 bp apart join one run (added 2b-i, see prereg section 3)
REF_OK = 0.15           # walk branch support, x D (oracle_r2_walk)
D_MIN = _dw.D_MIN        # no place below this depth (fixed)
D_MAX = _dw.upper_bound()  # no place above this depth: none unless --depth-cap / EVITEA_DEPTH_CAP sets one
INDEL_MAX = 50          # S2 small indel limit (fixed)
EXT_MAX = 600           # S3 extension length per side (fixed)
STOP_FRAC, STOP_MIN = 0.10, 3   # S3 extension stops when best letter < max(3, 0.10 x support at break) (G4c rule)
FORK = 0.60             # S3 extension stops when the top letter carries < 60% of the support (synthetic: 0.5-0.8)
PAIR_LO, PAIR_HI = -20, 50      # S3 two breaks pair when delta in [-20, +50] (fixed)
HOM_MAX, HET_MAX = 0.15, 0.85   # S3 genotype from the reference-path fraction at the break (G4b)
SEG_MIN = 60            # S4 minimum segment length (fixed)
ANCHOR_RC = 2.0         # S4 abstain if an anchor's read count > 2 x D (fixed)
REP_K, REP_FLANK, REP_MAX = 21, 500, 0.50   # S4 REP definition (class_strata cs_lib.rep_fracs)
LOCAL_FLANK, LOCAL_ID, LOCAL_COV = 2000, 0.90, 0.80   # S4 locally explained (fixed)
NUMT_ID, NUMT_COV = 0.80, 0.80                        # S4 NUMT (fixed)
SD_ID, SD_COV = 0.95, 0.80                            # S4 segmental duplication, 1 or 2 loci (fixed)
SAT_PMAX, SAT_COPIES, SAT_COV, SAT_WIN_ID = 50, 3, 0.80, 0.85   # S4 satellite (fixed; window identity added 2b-i)
LC_PMAX, LC_RUN = 3, 10  # S5 low-complexity mask: runs with s[i]==s[i+p], p <= 3, length >= 10 (added 2b-i)
GK = 15                 # S5 canonical word length (synthetic)
G_SHARE = 0.20          # S5 shared words / smaller word set (synthetic)
G_ALN, G_ID = 50, 0.80  # S5 alignment confirmation: >= 50 columns at >= 80% identity (synthetic)
BAND = 16               # S5 banded alignment half-width around a shared-word diagonal (added 2b-i)
SW_MATCH, SW_MISM, SW_GAP = 2, -3, -5   # S5 local alignment scores (added 2b-i)
TSD_LO, TSD_HI = 2, 20  # S6 V-TSD exact TSD length (fixed)
MERGE_D = 500           # S7 one call per place; family calls added where EviTEA has no call within 500 bp (fixed)

_C = str.maketrans("ACGT", "TGCA")
_DEL = str.maketrans("", "", "ACGT")


def rc(s):
    return s.translate(_C)[::-1]


def is_acgt(w):
    return len(w) > 0 and w.translate(_DEL) == ""


# ---------------------------------------------------------------- counters
class Counter:
    """msBWT counter, both strands; cache may be cleared between places to bound memory."""
    def __init__(self, bwt, cache=True):
        self.b, self.cache, self.n, self.use = bwt, {}, 0, cache

    def cnt(self, w):
        v = self.cache.get(w)
        if v is None:
            if not is_acgt(w):
                v = 0
            else:
                self.n += 2
                v = int(self.b.countOccurrencesOfSeq(w)) + int(self.b.countOccurrencesOfSeq(rc(w)))
            if self.use:
                self.cache[w] = v
        return v

    def clear(self):
        self.cache.clear()


class HashCounter:
    """Exact counter from a read list for words of the given lengths (tests and synthetic checks only)."""
    def __init__(self, reads, ks=(30, 31)):
        self.occ, self.n = {}, 0
        for r in reads:
            r = r.upper()
            for k in ks:
                for i in range(len(r) - k + 1):
                    w = r[i:i + k]
                    if is_acgt(w):
                        self.occ[w] = self.occ.get(w, 0) + 1

    def cnt(self, w):
        if not is_acgt(w):
            return 0
        self.n += 2
        return self.occ.get(w, 0) + self.occ.get(rc(w), 0)

    def clear(self):
        pass


# ---------------------------------------------------------------- S1 scan
def scan_counts(C, ref, k=K, stride=STRIDE):
    """Counts of the reference word starting every `stride` bases; -1 for words with a non-ACGT letter."""
    n = max(0, (len(ref) - k) // stride + 1)
    pos = np.arange(n, dtype=np.int64) * stride
    c = np.full(n, -1, dtype=np.int64)
    for j in range(n):
        p = j * stride
        w = ref[p:p + k]
        if is_acgt(w):
            c[j] = C.cnt(w)
    return pos, c


def local_depth(c, flank=DFLANK, stride=STRIDE, chunk=20000):
    """D[j] = median of valid counts among samples j-h..j+h (h = flank // stride, clipped to the arm);
    nfrac[j] = invalid samples / in-range samples in that window. D is NaN where no sample is valid."""
    h = flank // stride
    n = len(c)
    x = c.astype(float)
    x[c < 0] = np.nan
    xp = np.concatenate([np.full(h, np.nan), x, np.full(h, np.nan)])
    inr = np.concatenate([np.zeros(h), np.ones(n), np.zeros(h)])
    D = np.full(n, np.nan)
    nf = np.zeros(n)
    sw = np.lib.stride_tricks.sliding_window_view
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        win = sw(xp[s:e + 2 * h], 2 * h + 1)
        rin = sw(inr[s:e + 2 * h], 2 * h + 1).sum(axis=1)
        valid = np.sum(~np.isnan(win), axis=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            D[s:e] = np.nanmedian(win, axis=1)
        nf[s:e] = (rin - valid) / rin
    return D, nf


def candidates(pos, c, D, nf):
    """Runs of low sampled words: c < LOWF*D, valid word, window N fraction <= NFRAC. Low words <= RUN_GAP apart join one run;
    a run needs >= MIN_LOW low words. Returns dicts with first/last low word start, n_low and D (median of D over the run)."""
    with np.errstate(invalid="ignore"):
        low = (c >= 0) & (nf <= NFRAC) & ~np.isnan(D) & (c < LOWF * D)
    idx = np.flatnonzero(low)
    out = []
    if len(idx) == 0:
        return out
    g0 = 0
    for t in range(1, len(idx) + 1):
        if t == len(idx) or pos[idx[t]] - pos[idx[t - 1]] > RUN_GAP:
            g = idx[g0:t]
            if len(g) >= MIN_LOW:
                out.append(dict(first=int(pos[g[0]]), last=int(pos[g[-1]]), n_low=int(len(g)),
                                D=float(np.median(D[g])), j0=int(g[0]), j1=int(g[-1])))
            g0 = t
    return out


# ---------------------------------------------------------------- S2/S3 walk
def _occ(s, w, lo, hi):
    """Start positions j of w in s with lo <= j <= hi. The search is bounded at hi + len(w): an unbounded str.find scanned the
    rest of a 28-Mb arm whenever w was absent from the window (F-035; same result, ~100x slower)."""
    out, lo = [], max(0, lo)
    if hi < lo:
        return out
    end = hi + len(w)                       # > 0 here (a negative end would count from the string's end)
    j = s.find(w, lo, end)
    while j != -1:
        out.append(j)
        j = s.find(w, j + 1, end)
    return out


def _occ_unbounded(s, w, lo, hi):
    """The first version of _occ (reference for the equivalence test only)."""
    out, j = [], s.find(w, max(0, lo))
    while j != -1 and j <= hi:
        out.append(j)
        j = s.find(w, j + 1)
    return out


def indel_rejoin(C, s, i, ctx, alt, D):
    """Greedy most-supported extension from ctx+alt, up to INDEL_MAX+K letters; after each letter test whether the last K letters of
    the sample path occur in s at j (i-INDEL_MAX-K <= j <= i+INDEL_MAX+m) with |indel| = |m - (j+K-i)| <= INDEL_MAX.
    Returns (new_ctx, resume_index, indel) or None. Ties: smallest |indel|, then smallest j."""
    thr = REF_OK * D
    path, m = ctx + alt, 1
    while True:
        w = path[-K:]
        best = None
        for j in _occ(s, w, i - INDEL_MAX - K, i + INDEL_MAX + m):
            ind = m - (j + K - i)
            if abs(ind) <= INDEL_MAX:
                key = (abs(ind), j)
                if best is None or key < best[0]:
                    best = (key, j, ind)
        if best is not None:
            j = best[1]
            return s[j + 1:j + K], j + K, best[2]
        if m >= INDEL_MAX + K:
            return None
        cur = path[-(K - 1):]
        n = [C.cnt(cur + b) for b in "ACGT"]
        bi = max(range(4), key=lambda q: n[q])
        if n[bi] < thr:
            return None
        path += "ACGT"[bi]
        m += 1


def walk(C, s, i, ctx, end, D, shortcut=True):
    """Reference-guided walk on strand string s from index i, ctx = the 30 sample letters before i, until index end.
    Rules as g4b oracle_r2_walk.classify (branch support REF_OK*D; a supported alternative that rejoins the reference within 3
    letters is a substitution) plus small-indel rejoin (indel_rejoin). Exact shortcut: when cnt(ctx) - cnt(ctx+ref) < REF_OK*D no
    alternative branch can reach REF_OK*D (sum_b cnt(ctx+b) <= cnt(ctx)), so the alternatives are not queried.
    status: through (reached end) | break (unexplained divergence; i = break index) | dead (no supported branch) | refN."""
    thr = REF_OK * D
    n_snp = n_ind = 0
    guard = 4 * (end - i) + 1000
    while i < end:
        guard -= 1
        if guard < 0:
            return dict(status="dead", i=i, n_snp=n_snp, n_ind=n_ind, why="guard")
        rb = s[i]
        if not is_acgt(s[i:i + 4]) or i + 4 > len(s):
            return dict(status="refN", i=i, n_snp=n_snp, n_ind=n_ind)
        nref = C.cnt(ctx + rb)
        if shortcut and nref >= thr and C.cnt(ctx) - nref < thr:
            ctx = ctx[1:] + rb
            i += 1
            continue
        n = {b: (nref if b == rb else C.cnt(ctx + b)) for b in "ACGT"}
        T = sum(n.values())
        alt = max((b for b in "ACGT" if b != rb), key=lambda b: n[b])
        if n[alt] >= thr:
            if all(C.cnt((ctx + alt + s[i + 1:i + 1 + j])[-K:]) >= thr for j in (1, 2, 3)):
                n_snp += 1
                ctx = ctx[1:] + (rb if nref >= thr else alt)
                i += 1
                continue
            r = indel_rejoin(C, s, i, ctx, alt, D)
            if r is None:
                return dict(status="break", i=i, ctx=ctx, alt=alt, n_alt=n[alt], nref=nref, T=T,
                            frac=nref / T if T else 0.0, n_snp=n_snp, n_ind=n_ind)
            n_ind += 1
            if nref >= thr:
                ctx = ctx[1:] + rb
                i += 1
            else:
                ctx, i = r[0], r[1]
            continue
        if nref >= thr:
            ctx = ctx[1:] + rb
            i += 1
            continue
        return dict(status="dead", i=i, n_snp=n_snp, n_ind=n_ind, why="no_branch")
    return dict(status="through", i=i, n_snp=n_snp, n_ind=n_ind)


def extend(C, ctx, alt, n_alt, maxlen=EXT_MAX, fork=FORK):
    """Most-supported extension from ctx+alt (G4c junction.extend) with the fork rule. Returns (segment, stop reason)."""
    seq, cur = alt, (ctx + alt)[-(K - 1):]
    stop = max(STOP_MIN, STOP_FRAC * n_alt)
    while len(seq) < maxlen:
        n = [C.cnt(cur + b) for b in "ACGT"]
        bi = max(range(4), key=lambda q: n[q])
        T = sum(n)
        if n[bi] < stop:
            return seq, "support"
        if n[bi] < fork * T:
            return seq, "fork"
        seq += "ACGT"[bi]
        cur = cur[1:] + "ACGT"[bi]
    return seq, "max"


def left_walk(C, ref, cand):
    """S2: walk from the last sampled word before the run to last + K + INDEL_MAX."""
    a = cand["first"] - STRIDE
    if a < 0 or not is_acgt(ref[a:a + K]):
        return dict(status="edge", i=a), a
    end = min(len(ref) - 4, cand["last"] + K + INDEL_MAX)
    return walk(C, ref, a + K, ref[a + 1:a + K], end, cand["D"]), a


def right_walk(C, ref, cand):
    """S3: mirror walk on the reverse strand from the first sampled word after the run to first - INDEL_MAX - 1.
    Returns (walk dict with forward break bR when status == break, b, window lo, window hi)."""
    b = cand["last"] + STRIDE
    if b + K > len(ref) or not is_acgt(ref[b:b + K]):
        return dict(status="edge"), b, 0, 0
    lo = max(0, cand["first"] - 2 * INDEL_MAX - 2 * K - 10)
    hi = min(len(ref), b + K + 2 * INDEL_MAX + K + 10)
    R = rc(ref[lo:hi])
    r0 = hi - b                                    # rc index right after the anchor word
    end_r = min(len(R) - 4, hi - (cand["first"] - INDEL_MAX))
    w = walk(C, R, r0, R[r0 - K + 1:r0], end_r, cand["D"])
    if w["status"] == "break":
        w["bR"] = (hi - 1 - w["i"]) + 1
    return w, b, lo, hi


def genotype(frac):
    return "hom" if frac <= HOM_MAX else ("het" if frac <= HET_MAX else "low")


# ---------------------------------------------------------------- S4 helpers
def satellite(seg):
    """Tandem-periodic: for some period p <= SAT_PMAX, windows of length max(2p, 20) whose positions match at lag p in
    >= SAT_WIN_ID of cases, each extended by p, cover >= SAT_COV of the segment, and the periodic stretch spans >= SAT_COPIES copies.
    Also low complexity: the most common letter covers >= 80% of the segment."""
    L = len(seg)
    if L == 0:
        return False
    if max(seg.count(b) for b in "ACGT") >= 0.80 * L:
        return True
    a = np.frombuffer(seg.encode(), dtype=np.uint8)
    ok = np.isin(a, np.frombuffer(b"ACGT", dtype=np.uint8))
    for p in range(1, SAT_PMAX + 1):
        w = max(2 * p, 20)
        if L < SAT_COPIES * p or L - p < w:
            continue
        m = ((a[:L - p] == a[p:]) & ok[:L - p]).astype(np.int64)
        cs = np.concatenate([[0], np.cumsum(m)])
        good = np.flatnonzero(cs[w:] - cs[:-w] >= SAT_WIN_ID * w)      # window starts i, 0 <= i <= L-p-w
        if len(good) == 0:
            continue
        diff = np.zeros(L + 1, dtype=np.int64)
        np.add.at(diff, good, 1)
        np.add.at(diff, np.minimum(L, good + w + p), -1)
        cov = np.cumsum(diff[:L]) > 0
        if cov.sum() < SAT_COV * L:
            continue
        # longest covered stretch must span >= SAT_COPIES copies
        edges = np.flatnonzero(np.diff(np.concatenate([[0], cov.astype(np.int8), [0]])))
        longest = int(np.max(edges[1::2] - edges[0::2]))
        if longest >= SAT_COPIES * p:
            return True
    return False


def union_cov(intervals, L):
    """Fraction of [0, L) covered by the union of [a, b) intervals."""
    if L <= 0:
        return 0.0
    iv = sorted(t for t in ((max(0, a), min(L, b)) for a, b in intervals) if t[1] > t[0])
    tot, ca, cb = 0, None, None
    for a, b in iv:
        if ca is None:
            ca, cb = a, b
        elif a <= cb:
            cb = max(cb, b)
        else:
            tot += cb - ca
            ca, cb = a, b
    if ca is not None:
        tot += cb - ca
    return tot / L


def local_rule(hits, L):
    """hits: (q_st, q_en, mlen, blen) against ref[x-2000:x+2000]."""
    return union_cov([(a, b) for a, b, ml, bl in hits if bl and ml / bl >= LOCAL_ID], L) >= LOCAL_COV


def numt_rule(hits, L):
    return union_cov([(a, b) for a, b, ml, bl in hits if bl and ml / bl >= NUMT_ID], L) >= NUMT_COV


def sd_rule(hits, L):
    """hits: (q_st, q_en, mlen, blen, ctg, r_st) against the genome without the mitochondrial contig. SD if 1 or 2 distinct loci align at >= 95% identity
    over >= 80% of the segment."""
    loci = {(h[4], h[5] // 1000) for h in hits if h[3] and h[2] / h[3] >= SD_ID and (h[1] - h[0]) >= SD_COV * L}
    return 1 <= len(loci) <= 2


class MappyFirewall:
    """S4 alignment-based rules with mappy (minimap2). genome: dict name -> sequence (upper case) without the mitochondrial contig;
    chrm: the mitochondrial sequence (None: no mitochondrial test)."""
    def __init__(self, genome_fa, chrm, preset="map-ont"):
        """genome_fa: FASTA of the whole reference (the mitochondrial contig included; NUMTs are tested first); chrm: the sequence of
        the mitochondrial contig, or None when the reference has none."""
        import mappy
        self.mp, self.preset = mappy, preset
        self.g = mappy.Aligner(genome_fa, preset=preset, best_n=5)
        if not self.g:
            raise RuntimeError("mappy index failed: " + genome_fa)
        self.m = mappy.Aligner(seq=chrm, preset=preset) if chrm else None

    def _hits(self, al, seg):
        if al is None:
            return []
        return [(h.q_st, h.q_en, h.mlen, h.blen, h.ctg, h.r_st) for h in al.map(seg)]

    def local(self, seg, ref, x):
        w = ref[max(0, x - LOCAL_FLANK):x + LOCAL_FLANK]
        if len(w) < 50:
            return False
        al = self.mp.Aligner(seq=w, preset=self.preset)
        return local_rule([h[:4] for h in self._hits(al, seg)], len(seg))

    def numt(self, seg):
        return numt_rule([h[:4] for h in self._hits(self.m, seg)], len(seg))

    def sd(self, seg):
        return sd_rule(self._hits(self.g, seg), len(seg))


# ---------------------------------------------------------------- reference word multiplicity (S4)
_ENC = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _ENC[_b] = _i


def kmer_codes(seq, k):
    """(canonical codes uint64, valid mask) for every k-word start of seq (upper case); k <= 32."""
    a = _ENC[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(a) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=bool)
    bad = np.concatenate([[0], np.cumsum(a == 4)])
    valid = (bad[k:k + n] - bad[:n]) == 0
    v = np.where(a == 4, 0, a).astype(np.uint64)
    f = np.zeros(n, dtype=np.uint64)
    r = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        f = (f << np.uint64(2)) | v[j:j + n]
        r = r | ((np.uint64(3) - v[j:j + n]) << np.uint64(2 * j))
    return np.minimum(f, r), valid


def multi_codes(genome, k):
    """Sorted unique canonical codes of k-words that occur >= 2 times in the genome counting both strands."""
    parts = []
    for s in genome.values():
        c, v = kmer_codes(s.upper(), k)
        parts.append(c[v])
    allc = np.concatenate(parts) if parts else np.zeros(0, dtype=np.uint64)
    allc.sort()
    if len(allc) < 2:
        return np.zeros(0, dtype=np.uint64)
    dup = allc[1:] == allc[:-1]
    return np.unique(allc[1:][dup])


class RefMulti:
    def __init__(self, m21, m31):
        self.m = {REP_K: m21, K: m31}

    def _isin(self, codes, k):
        arr = self.m[k]
        if len(arr) == 0 or len(codes) == 0:
            return np.zeros(len(codes), dtype=bool)
        ix = np.searchsorted(arr, codes)
        ix[ix >= len(arr)] = len(arr) - 1
        return arr[ix] == codes

    def word_multi(self, w):
        c, v = kmer_codes(w.upper(), len(w))
        return bool(v[0] and self._isin(c[:1], len(w))[0])

    def rep_frac(self, ref, x, flank=REP_FLANK):
        seg = ref[max(0, x - flank):x + flank]
        c, v = kmer_codes(seg.upper(), REP_K)
        c = c[v]
        if len(c) == 0:
            return None
        return float(np.mean(self._isin(c, REP_K)))


# ---------------------------------------------------------------- per-place pipeline (S2-S4)
def process(C, ref, cand, rm, fw, chrom="", fork=FORK):
    """S2-S4 for one candidate. Returns a place dict with 'stage' = where it stopped:
    S2_through | S2_dead | S2_refN | S2_edge | S4_abstain | S4_reject | S4_pass."""
    D = cand["D"]
    pl = dict(chrom=chrom, first=cand["first"], last=cand["last"], n_low=cand["n_low"], D=round(D, 2))
    lw, a = left_walk(C, ref, cand)
    if lw["status"] != "break":
        pl["stage"] = "S2_" + lw["status"]
        return pl
    pl["bL"] = lw["i"]
    # S4 abstention checks that need no segment (done first to save walking)
    why = []
    if not (D_MIN <= D <= D_MAX):
        why.append("depth")
    b = cand["last"] + STRIDE
    for nm, p in (("anchorL", a), ("anchorR", b)):
        w = ref[p:p + K]
        if len(w) < K or not is_acgt(w) or rm.word_multi(w) or C.cnt(w) > ANCHOR_RC * D:
            why.append(nm)
    rf = rm.rep_frac(ref, lw["i"])
    pl["rep"] = None if rf is None else round(rf, 3)
    if rf is not None and rf >= REP_MAX:
        why.append("rep")
    if why:
        pl.update(stage="S4_abstain", why=",".join(why))
        return pl
    # S3 two-sided walk and extraction
    lseg, lstop = extend(C, lw["ctx"], lw["alt"], lw["n_alt"], fork=fork)
    rw, _, _, _ = right_walk(C, ref, cand)
    pl.update(fracL=round(lw["frac"], 3), lseg=lseg, lstop=lstop, rstat=rw["status"])
    rseg = ""
    if rw["status"] == "break":
        rs, rstop = extend(C, rw["ctx"], rw["alt"], rw["n_alt"], fork=fork)
        rseg = rc(rs)                                # forward orientation, ends at the insert's right end
        pl.update(bR=rw["bR"], fracR=round(rw["frac"], 3), rstop=rstop)
    pl["rseg"] = rseg
    two = rw["status"] == "break" and PAIR_LO <= pl["bL"] - rw["bR"] <= PAIR_HI
    pl["two_sided"] = int(two)
    pl["delta"] = pl["bL"] - rw["bR"] if rw["status"] == "break" else None
    pl["tsd"] = max(0, pl["delta"]) if two else None
    pl["geno"] = genotype(min(lw["frac"], rw["frac"]) if rw["status"] == "break" else lw["frac"])
    pl["pos"] = pl["bL"]
    # S4 firewall per segment
    keep, rej = [], []
    for side, sg in (("L", lseg), ("R", rseg)):
        if len(sg) < SEG_MIN:
            if sg:
                rej.append(side + ":short")
            continue
        r = None
        if fw.local(sg, ref, pl["bL"]):
            r = "local"
        elif fw.numt(sg):
            r = "numt"
        elif fw.sd(sg):
            r = "sd"
        elif satellite(sg):
            r = "sat"
        if r:
            rej.append(side + ":" + r)
        else:
            keep.append(side)
    pl["seg_keep"] = ",".join(keep)
    pl["seg_reject"] = ",".join(rej)
    pl["stage"] = "S4_pass" if keep else "S4_reject"
    return pl


def target_candidate(pos, D, half=50):
    """Targeted mode: a pseudo-run centred on pos (used for hostile sites off the arm and for synthetic checks)."""
    f = (pos - half) // STRIDE * STRIDE
    return dict(first=f, last=f + 2 * half // STRIDE * STRIDE, n_low=0, D=D)


def local_D(C, ref, pos, flank=DFLANK, stride=STRIDE):
    """Median of valid sampled-word counts within +-flank of pos (targeted mode)."""
    v = []
    for p in range((max(0, pos - flank) // stride) * stride, min(len(ref) - K, pos + flank) + 1, stride):
        w = ref[p:p + K]
        if is_acgt(w):
            v.append(C.cnt(w))
    return float(np.median(v)) if v else float("nan")


# ---------------------------------------------------------------- S5 grouping
def mask_lc(s):
    """Replace low-complexity runs (s[i] == s[i+p] for p <= LC_PMAX over >= LC_RUN positions) by N."""
    L = len(s)
    out = list(s)
    for p in range(1, LC_PMAX + 1):
        i = 0
        while i < L - p:
            if s[i] == s[i + p]:
                j = i
                while j < L - p and s[j] == s[j + p]:
                    j += 1
                if j - i >= LC_RUN:
                    for q in range(i, min(L, j + p)):
                        out[q] = "N"
                i = j + 1
            else:
                i += 1
    return "".join(out)


def kset(s, k=GK):
    c, v = kmer_codes(s, k)
    return set(c[v].tolist())


def _sw(a, b, band=None, diag=0):
    """Local alignment, linear gap (SW_MATCH, SW_MISM, SW_GAP). Each cell carries (score, matches, columns) of its best path;
    ties prefer diagonal, then up, then left. band: only cells with |(j - i) - diag| <= band. Returns best (score, matches, cols);
    ties: first cell in row-major order. N never matches."""
    n, m = len(a), len(b)
    best = (0, 0, 0)
    prev = [(0, 0, 0)] * (m + 1)
    for i in range(1, n + 1):
        cur = [(0, 0, 0)] * (m + 1)
        ai = a[i - 1]
        if band is None:
            jlo, jhi = 1, m
        else:
            jlo, jhi = max(1, i + diag - band), min(m, i + diag + band)
        for j in range(jlo, jhi + 1):
            mt = ai == b[j - 1] and ai != "N"
            d = prev[j - 1]
            ds = d[0] + (SW_MATCH if mt else SW_MISM)
            u = prev[j]
            us = u[0] + SW_GAP
            lf = cur[j - 1]
            ls = lf[0] + SW_GAP
            if ds >= us and ds >= ls and ds > 0:
                cell = (ds, d[1] + (1 if mt else 0), d[2] + 1)
            elif us >= ls and us > 0:
                cell = (us, u[1], u[2] + 1)
            elif ls > 0:
                cell = (ls, lf[1], lf[2] + 1)
            else:
                cell = (0, 0, 0)
            cur[j] = cell
            if cell[0] > best[0]:
                best = cell
        prev = cur
    return best


def aln_ok(res):
    sc, mt, cols = res
    return cols >= G_ALN and mt >= G_ID * cols


def confirm(a, b, k=GK, top=2):
    """Banded alignment around the top shared-word diagonals, both strands. True if the best alignment on some tried diagonal has
    >= G_ALN columns at >= G_ID identity."""
    for bb in (b, rc(b)):
        pos_a = {}
        for i in range(len(a) - k + 1):
            w = a[i:i + k]
            if is_acgt(w):
                pos_a.setdefault(w, []).append(i)
        dc = {}
        for j in range(len(bb) - k + 1):
            for i in pos_a.get(bb[j:j + k], ()):
                dc[j - i] = dc.get(j - i, 0) + 1
        for d, _ in sorted(dc.items(), key=lambda t: (-t[1], t[0]))[:top]:
            if aln_ok(_sw(a, bb, band=BAND, diag=d)):
                return True
    return False


def confirm_brute(a, b):
    """Full local alignment, both strands (reference for confirm)."""
    return aln_ok(_sw(a, b)) or aln_ok(_sw(a, rc(b)))


class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def f(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def u(self, a, b):
        a, b = self.f(a), self.f(b)
        if a != b:
            self.p[max(a, b)] = min(a, b)
            return True
        return False


def group(places, confirm_fn=confirm):
    """S5. places: list of dicts with 'segs' (list of kept segments). Returns (cluster id per place, number of alignments run).
    Join two places if some pair of their segments shares >= G_SHARE x the smaller canonical-word set and confirm_fn agrees.
    Pairs are tried in order of decreasing share (ties: place indices); pairs already connected are skipped (single linkage)."""
    segs = []
    for pi, p in enumerate(places):
        for s in p["segs"]:
            m = mask_lc(s)
            segs.append((pi, m, kset(m)))
    idx = {}
    for si, (_, _, ks) in enumerate(segs):
        for w in ks:
            idx.setdefault(w, []).append(si)
    shared = {}
    for si, (pi, _, ks) in enumerate(segs):
        for w in ks:
            for sj in idx[w]:
                if sj > si and segs[sj][0] != pi:
                    shared[(si, sj)] = shared.get((si, sj), 0) + 1
    pairs = []
    for (si, sj), n in shared.items():
        mn = min(len(segs[si][2]), len(segs[sj][2]))
        if mn and n >= G_SHARE * mn:
            pairs.append((-n / mn, segs[si][0], segs[sj][0], si, sj))
    pairs.sort()
    uf = UF(len(places))
    nal = 0
    for _, pa, pb, si, sj in pairs:
        if uf.f(pa) == uf.f(pb):
            continue
        nal += 1
        if confirm_fn(segs[si][1], segs[sj][1]):
            uf.u(pa, pb)
    return [uf.f(i) for i in range(len(places))], nal


def dedupe(places):
    """One place per insertion (S7 'one call per place'): places with a left break (stage S4_*) on one chromosome whose breaks
    chain within MERGE_D are one place. Representative: S4_pass two-sided > S4_pass > S4_reject > S4_abstain; ties: leftmost
    break, then input order. The others get dup = 1 and are left out of grouping. Two candidates of one insertion often break at
    the same letter (seen in SYN-FAM); without this they could form a two-place 'family' from one insertion."""
    rank = {"S4_pass": 1, "S4_reject": 2, "S4_abstain": 3}
    idx = [i for i, p in enumerate(places) if p.get("stage") in rank and p.get("bL") is not None]
    idx.sort(key=lambda i: (places[i]["chrom"], places[i]["bL"], i))
    groups, cur = [], []
    for i in idx:
        if cur and (places[i]["chrom"] != places[cur[-1]]["chrom"] or places[i]["bL"] - places[cur[-1]]["bL"] > MERGE_D):
            groups.append(cur); cur = []
        cur.append(i)
    if cur:
        groups.append(cur)
    for g in groups:
        best = min(g, key=lambda i: (rank[places[i]["stage"]] - (0.5 if places[i]["stage"] == "S4_pass" and places[i].get("two_sided") else 0),
                                     places[i]["bL"], i))
        for i in g:
            places[i]["dup"] = int(i != best)
            places[i]["_rep"] = places[best]
    return places


# ---------------------------------------------------------------- S6 family rule, S7 output
def families(places, cl):
    """S6. Returns {cluster: {'n': places, 'two': two-sided places, 'tsd': places with TSD in [TSD_LO, TSD_HI],
    'pair': V-pair family, 'vtsd': V-TSD family}}."""
    out = {}
    for p, c in zip(places, cl):
        f = out.setdefault(c, dict(n=0, two=0, tsd=0))
        f["n"] += 1
        f["two"] += p.get("two_sided", 0)
        t = p.get("tsd")
        f["tsd"] += int(t is not None and TSD_LO <= t <= TSD_HI)
    for f in out.values():
        f["pair"] = int(f["n"] >= 2 and f["two"] >= 2)
        f["vtsd"] = int(f["pair"] and f["tsd"] >= 2)
    return out


def classify_place(p, fam):
    """S7 track of one grouped place: family | orphan (two-sided, alone) | singleton (one-sided, alone) | member (cluster that is not
    a family)."""
    if fam["pair"]:
        return "family"
    if fam["n"] == 1:
        return "orphan" if p.get("two_sided") else "singleton"
    return "member"
