"""Synthetic reference and reads for the library-free tests (the scan itself requires no msBWT or mappy)."""
import random

from evitea.library_free import _p4_lib as P

BASES = "ACGT"


class NoFirewall:
    """Segment firewall that rejects nothing."""

    def local(self, *a):
        return False

    def numt(self, *a):
        return False

    def sd(self, *a):
        return False


def rnd(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def build(seed=7, glen=36000, sites=((3000, 1.0), (6000, 1.0), (9000, 0.5), (12000, 0.5), (15000, 0.5),
                                      (18000, 0.25), (21000, 0.25), (24000, 0.25), (27000, 0.25), (30000, 0.5)),
          tsd=5, ins_len=400, depth=30, rl=125):
    """Return ref, sites and a counter over mixed-haplotype reads. The allele with an insertion at x is ref[:x] + INS + ref[x-tsd:]."""
    rng = random.Random(seed)
    ref = rnd(rng, glen)
    inserts = {x: rnd(rng, ins_len) for x, _ in sites}
    reads = []
    for _ in range(glen * depth // rl):
        p = rng.randrange(0, glen - rl)
        r = ref[p:p + rl]
        for x, f in sites:
            if p < x < p + rl:
                if rng.random() < f:
                    r = (ref[:x] + inserts[x] + ref[x - tsd:])[p:p + rl]
                break
        reads.append(P.rc(r) if rng.random() < 0.5 else r)
    for x, f in sites:
        allele = ref[:x] + inserts[x] + ref[x - tsd:]
        for _ in range(int((ins_len - 1) * depth * f / rl)):
            q = rng.randrange(x + 1, x + ins_len)
            r = allele[q:q + rl]
            reads.append(P.rc(r) if rng.random() < 0.5 else r)
    C = P.HashCounter(reads, ks=(30, 31, 41))
    C.use = False
    C.cache = {}
    return ref, sites, C


def build_snp2(seed=21, glen=40000, specs=(), tsd=5, ins_len=400, depth=30, rl=125):
    """specs: dicts(x, f, snps=[offsets before x], ins=bool, no_left=bool). The alt haplotype (fraction f) carries SNPs at x - offset and,
    if ins, an insertion at x (ref[:x] + INS + ref[x-tsd:]). no_left: no alt read spans the left junction. Returns ref, counter."""
    rng = random.Random(seed)
    ref = rnd(rng, glen)
    alleles = []
    for sp in specs:
        x, f = sp["x"], sp["f"]
        a = ref
        for off in sp.get("snps", ()):
            alt = rng.choice([b for b in BASES if b != ref[x - off]])
            a = a[:x - off] + alt + a[x - off + 1:]
        ins = rnd(rng, ins_len) if sp.get("ins") else ""
        if ins:
            a = a[:x] + ins + a[x - tsd:]
        offs = list(sp.get("snps", ()))
        lo = x - max(offs) if offs else x
        hi = x + 1 if ins else (x - min(offs) + 1)
        alleles.append((x, f, a, lo, hi, bool(ins), sp.get("no_left", False)))
    reads = []
    for _ in range(glen * depth // rl):
        p = rng.randrange(0, glen - rl)
        r = ref[p:p + rl]
        for (x, f, a, lo, hi, wi, nl) in alleles:
            if p < hi and p + rl > lo:
                if rng.random() < f and not (nl and p < x < p + rl):
                    r = a[p:p + rl]
                break
        reads.append(P.rc(r) if rng.random() < 0.5 else r)
    for (x, f, a, lo, hi, wi, nl) in alleles:
        if wi:
            for _ in range(int((ins_len - 1) * depth * f / rl)):
                q = rng.randrange(x + 1, x + ins_len)
                r = a[q:q + rl]
                reads.append(P.rc(r) if rng.random() < 0.5 else r)
    C = P.HashCounter(reads, ks=(30, 31, 41))
    C.use = False
    C.cache = {}
    return ref, C
