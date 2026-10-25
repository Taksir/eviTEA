"""Synthetic sample for the anchored-extension tests: a random reference, a P-like element with a perfect 31-bp terminal
inverted repeat, a LINE-like element with 5' truncations and a poly-A tail, error-free 125-bp reads from both strands."""
import random
from evitea.library_free import _p4_lib as P

ACGT = "ACGT"


def rnd(r, n):
    return "".join(r.choice(ACGT) for _ in range(n))


class SubCounter:
    """Count a k-mer and its reverse complement in the reads (same semantics as P.Counter).

    Uses an exact table for 30- and 31-mers and an overlapping text search for any other length."""
    def __init__(self, reads):
        self.text = "|".join(reads)
        self.occ = {}
        for rd in reads:
            for k in (30, 31):
                for i in range(len(rd) - k + 1):
                    w = rd[i:i + k]
                    self.occ[w] = self.occ.get(w, 0) + 1
        self.cache, self.use, self.n = {}, True, 0

    def _occ(self, w):
        if len(w) in (30, 31):
            return self.occ.get(w, 0)
        n, j = 0, self.text.find(w)
        while j >= 0:
            n += 1
            j = self.text.find(w, j + 1)
        return n

    def cnt(self, w):
        v = self.cache.get(w)
        if v is None:
            v = 0 if not P.is_acgt(w) else self._occ(w) + self._occ(P.rc(w))
            self.n += 2
            self.cache[w] = v
        return v

    def clear(self):
        self.cache.clear()


def build(seed=7, ref_len=120_000, n_px=12, n_lx=8, n_uq=3, n_st=2, polya=36, n_arun=150, depth=30, rlen=125):
    r = random.Random(seed)
    ref = list(rnd(r, ref_len))
    # A runs of 32-40 letters over the reference, so that 30 A in a row are multi-copy as in a real genome
    for _ in range(n_arun):
        x = r.randrange(0, ref_len - 50)
        ref[x:x + 0] = []
        L = r.randint(32, 40)
        ref[x:x + L] = list("A" * L if r.random() < 0.5 else "T" * L)
    ref = "".join(ref)[:ref_len]
    tir = rnd(r, 31)
    px = tir + rnd(r, 1438) + P.rc(tir)
    lx_body = rnd(r, 3000)
    sites, x = [], 3000
    kinds = ["PX"] * n_px + ["LX"] * n_lx + ["UQ"] * n_uq + ["ST"] * n_st
    r.shuffle(kinds)
    for k in kinds:
        x += r.randint(3000, 4500)
        while "A" * 12 in ref[x - 150:x + 150] or "T" * 12 in ref[x - 150:x + 150]:
            x += 300
        sites.append(dict(kind=k, pos=x))
    assert sites[-1]["pos"] < ref_len - 3000, sites[-1]
    sample, last, ins = [], 0, []
    for s in sites:
        if s["kind"] == "PX":
            el, tsd = px, 8
        elif s["kind"] == "LX":
            el, tsd = lx_body[-r.randint(400, 2500):] + "A" * polya, 10
        elif s["kind"] == "UQ":                       # single-copy non-TE insertion: must stay exactly as before
            el, tsd = rnd(r, 300), 6
        else:                                         # satellite insertion: must stay rejected
            el, tsd = "AATAT" * 40, 5
        strand = r.choice("+-")
        el = el if strand == "+" else P.rc(el)
        p = s["pos"]
        sample.append(ref[last:p + tsd])
        sample.append(el)
        last = p                      # the TSD ref[p:p+tsd] appears on both sides
        s.update(strand=strand, insert=el, tsd=tsd)
    sample.append(ref[last:])
    sample = "".join(sample)
    n = int(depth * len(sample) / rlen)
    reads = []
    for _ in range(n):
        a = r.randrange(0, len(sample) - rlen)
        w = sample[a:a + rlen]
        reads.append(w if r.random() < 0.5 else P.rc(w))
    return dict(ref=ref, sample=sample, reads=reads, sites=sites, px=px, tir=tir, lx_body=lx_body)
