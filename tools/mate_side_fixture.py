"""Synthetic fixture of the mate-anchored search (no bowtie2, no msBWT): the same loci and reads run under either side rule.

    PYTHONPATH=<src> python tools/mate_side_fixture.py [--mate-side-rule strand] [--junction-side-rule split]

Prints one JSON document: hashes of the evidence, the counts (in order), the ledger and the calls, and the numbers the rule
changes (pairs with a side, calls). Run with the src folder of 18bae44 (no options) and with the current src (no options) to
compare the default byte for byte; the opt-in runs show the before/after count.

Loci (envelope = the span of the host mates, as the candidate step makes it): 6 with forward host mates only, 3 with reverse host
mates only, 3 with both strands (forward left of the insertion, reverse right of it). Every pair has a host mate placed uniquely in
the reference and a mate cut from the family's TE consensus.
"""

import hashlib
import json
import random
import sys

from evitea.incumbent import mate_rescue as MR
from evitea.incumbent.candidates import CandidateLocus
from evitea.incumbent.sequence import revcomp
from evitea.incumbent.te_library import LocalTEIndex, TERecord

SEED = 11
SPLIT_ROW = getattr(MR, "LocusSplit", ())     # 18bae44 has none


def build():
    rng = random.Random(SEED)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    te = rand(900)
    reference = {"chrS": rand(24000)}
    kmers = {te[i:i + 21] for i in range(len(te) - 20)} | {revcomp(te)[i:i + 21] for i in range(len(te) - 20)}
    loci, reads = [], {}
    pid = 0
    kinds = ["F"] * 6 + ["R"] * 3 + ["B"] * 3
    for index, kind in enumerate(kinds):
        site = 1500 + index * 1800                                  # the insertion point in the reference
        pairs, mids = [], []
        n = 8 if kind != "B" else 10
        for k in range(n):
            pid += 1
            strand = "F" if kind == "F" else "R" if kind == "R" else ("F" if k % 2 == 0 else "R")
            if strand == "F":
                start = site - 40 - rng.randrange(0, 260)          # forward host mates lie left of the site
            else:
                start = site + 40 + rng.randrange(0, 260)          # reverse host mates lie right of it
            host = reference["chrS"][start:start + 100]
            host = host if strand == "F" else revcomp(host)
            t0 = rng.randrange(0, len(te) - 100)
            tread = te[t0:t0 + 100]
            tread = tread if rng.random() < 0.5 else revcomp(tread)
            reads[pid] = (host, tread) if k % 2 == 0 else (tread, host)
            pairs.append(pid)
            mids.append(start + 50)
        loci.append(CandidateLocus(f"L{index:03d}", "chrS", "F1", float(min(mids)), float(max(mids)), tuple(pairs)))
    return loci, reads, reference, {"F1": kmers}, LocalTEIndex([TERecord("T1", "F1", te)])


def main(argv):
    side_rule = argv[argv.index("--mate-side-rule") + 1] if "--mate-side-rule" in argv else None
    junction = argv[argv.index("--junction-side-rule") + 1] if "--junction-side-rule" in argv else None
    loci, reads, reference, kmers, index = build()
    kw = {}
    if side_rule not in (None, "midpoint"):
        kw["side_rule"] = side_rule
    if junction not in (None, "midpoint"):
        kw["record_split"] = True
    evidence, counts = MR.mate_rescue_evidence(loci=loci, reads=reads, reference=reference, family_kmers=kmers, te_index=index, **kw)
    result = MR.finalize_mate_rescue("S", loci, [], evidence, counts)
    sides = {}
    for row in evidence:
        if isinstance(row, MR.MateEvidence):
            sides[(row.locus_id, row.pair_id, row.te_mate)] = row.side
    out = {
        "evidence_sha256": hashlib.sha256(repr([r for r in evidence if not isinstance(r, SPLIT_ROW)]).encode()).hexdigest(),
        "counts_in_order": list(counts.items()),
        "ledger_sha256": hashlib.sha256(repr(result.ledger).encode()).hexdigest(),
        "calls_sha256": hashlib.sha256(repr(result.calls).encode()).hexdigest(),
        "pairs_with_a_side": len(sides),
        "pairs_left": sum(1 for v in sides.values() if v == "LEFT"),
        "pairs_right": sum(1 for v in sides.values() if v == "RIGHT"),
        "loci_total": len(loci),
        "loci_with_two_sides_each_2plus": sum(1 for row in result.ledger if row.passes),
        "calls": len(result.calls),
    }
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
