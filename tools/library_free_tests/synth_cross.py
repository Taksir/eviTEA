"""Synthetic sample for the cross-species rules (``evitea.library_free.cross``): one scenario per rule and a control for each.

DIV       an insertion inside a block where the sample differs from the reference at one letter in 42 (about half the 31-letter
          reference words are broken); 1% of the reads come from a copy of the block without the substitutions (shared-word depth)
LOWCOV    an insertion inside a block read at a fifth of the depth (shared-word depth must leave it)
REP       an insertion between two copies of a reference repeat: the repeat content around the place is above REP_MAX, below
          that of the most repetitive windows (a block of copies elsewhere), and the local depth is that of single-copy words
          (repeat cutoff)
POLYA     a poly-A element on each strand whose tail length varies between reads (34-54 letters) (poly-A tails)
FAM       one element inserted at three places with target-site duplications (the family the recurring sequence needs)
RECUR     the same element inserted 1 kb from an old reference copy of it (93.5% identity), which the local test rejects
LOCALDUP  a 400-letter copy of the flank 1 kb upstream inserted with a duplication (a local copy that must stay rejected)
"""
import random

from evitea.library_free import _p4_lib as P

ACGT = "ACGT"


def rnd(r, n):
    return "".join(r.choice(ACGT) for _ in range(n))


def mutate(r, s, rate):
    """s with each letter substituted with probability ``rate``."""
    out = list(s)
    for i in range(len(s)):
        if r.random() < rate:
            out[i] = r.choice([x for x in ACGT if x != out[i]])
    return "".join(out)


def substitute(r, s, n):
    """s with ``n`` letters at distinct random places substituted."""
    out = list(s)
    for i in r.sample(range(len(s)), n):
        out[i] = r.choice([x for x in ACGT if x != out[i]])
    return "".join(out)


def build(seed=11, depth=30, rlen=125, n_var=20, undiverged=0.01):
    r = random.Random(seed)
    ref_len = 150_000
    ref = list(rnd(r, ref_len))
    for _ in range(150):                    # A runs, so that 30 A in a row are multi-copy as in a real genome
        x = r.randrange(0, ref_len - 50)
        L = r.randint(32, 40)
        ref[x:x + L] = list("A" * L if r.random() < 0.5 else "T" * L)
    # repeat block 30-34 kb: copies of R between 100-letter unique stretches (the most repetitive windows); at 40 kb two copies
    # of the first 300 letters of R on either side of a 200-letter unique stretch, where the REP place is
    R = rnd(r, 700)
    blk = "".join(R + rnd(r, 100) for _ in range(5))
    ref[30_000:30_000 + len(blk)] = list(blk)
    ref[39_800:40_100] = list(R[:300])
    ref[40_300:40_600] = list(R[:300])
    F = rnd(r, 400)
    ref[99_000:99_400] = list(substitute(r, F, 26))     # old copy of F at 93.5% identity, 1 kb before the RECUR place
    ref = "".join(ref)[:ref_len]

    def clean(p):
        while "A" * 12 in ref[p - 150:p + 150] or "T" * 12 in ref[p - 150:p + 150]:
            p += 37
        return p

    px_tir = rnd(r, 31)
    px = px_tir + rnd(r, 1200) + P.rc(px_tir)
    lx_body = rnd(r, 1500)
    rep_at = 40_200
    sites = [dict(kind="DIV", pos=clean(10_000), el=px, tsd=8, strand="+"),
             dict(kind="LOWCOV", pos=clean(24_000), el=px, tsd=8, strand="-"),
             dict(kind="REP", pos=rep_at, el=px, tsd=8, strand="+"),
             dict(kind="POLYA", pos=clean(52_000), el=None, tsd=10, strand="+"),
             dict(kind="POLYA", pos=clean(62_000), el=None, tsd=10, strand="-"),
             dict(kind="FAM", pos=clean(72_000), el=F, tsd=8, strand="+"),
             dict(kind="FAM", pos=clean(80_000), el=F, tsd=8, strand="-"),
             dict(kind="FAM", pos=clean(88_000), el=F, tsd=8, strand="+"),
             dict(kind="RECUR", pos=clean(100_400), el=F, tsd=8, strand="+"),
             dict(kind="LOCALDUP", pos=clean(120_000), el=None, tsd=8, strand="+")]
    for s in sites:
        if s["kind"] == "LOCALDUP":
            s["el"] = ref[s["pos"] - 1400:s["pos"] - 1000]
    div = (6_000, 14_000)
    low = (20_000, 28_000)
    ref_div = ref[:div[0]] + mutate(r, ref[div[0]:div[1]], 1 / 42) + ref[div[1]:]

    def variant(base, rr):
        """(sample, offsets of the DIV and LOWCOV blocks) with each poly-A tail drawn at 34-54 letters."""
        out, last, shift, offs = [], 0, 0, {}
        for s in sites:
            el = s["el"]
            if s["kind"] == "POLYA":
                el = lx_body + "A" * rr.randint(34, 54)
            el = el if s["strand"] == "+" else P.rc(el)
            p = s["pos"]
            for name, (a, b) in (("div", div), ("low", low)):
                if last <= a < p:
                    offs[name] = (a + shift, b + shift + (len(el) + s["tsd"] if b > p else 0))
            out.append(base[last:p + s["tsd"]])
            out.append(el)
            shift += len(el) + s["tsd"]
            last = p
        out.append(base[last:])
        return "".join(out), offs

    variants = [variant(ref if v == 0 else ref_div, random.Random(seed * 100 + v)) for v in range(n_var)]
    n = int(depth * len(variants[1][0]) / rlen)
    reads = []
    for _ in range(n):
        sample, offs = variants[0 if r.random() < undiverged else r.randrange(1, n_var)]
        a = r.randrange(0, len(sample) - rlen)
        lo, hi = offs["low"]
        if lo - rlen < a < hi and r.random() < 0.8:
            continue
        w = sample[a:a + rlen]
        reads.append(w if r.random() < 0.5 else P.rc(w))
    for s in sites:
        s.pop("el")
    return dict(ref=ref, reads=reads, sites=sites, F=F)
