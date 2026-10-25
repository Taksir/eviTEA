"""Small synthetic data set for `evitea call`: a reference, a TE library with its family map and reference copies,
and paired reads of a sample that carries non-reference TE insertions.

    python tools/synthetic/make_call_data.py OUTPUT_DIR [--read-length 125] [--fragment-mean 350] [--fragment-sd 40]
                                                        [--depth 30] [--seed 11] [--naming flybase|plain]
                                                        [--contigs chr2L,chr2R,chr3L,chr3R,chrX] [--mito chrM]

Files written to OUTPUT_DIR:
  reference.fa            contigs (default the five fly arms, 200 kb each) and a mitochondrial contig; each
                          contig holds 4 diverged TE copies and 14 near-identical copies of one high-copy record
  te_library.fa           120 TE records in 60 families (the last ten share their ends with the first)
  family_map.tsv          FlyBase-style table (FBte_ID, FlyBase_Family), or record<TAB>family with --naming plain
  repeat_annotation.bed   BED6 of the TE copies placed in the reference (name = TE record id)
  sample_R1.fastq, sample_R2.fastq
  truth_insertions.tsv    contig, position (0-based, before the TSD), TE record, family, strand

Everything follows from the seed; the same arguments give the same files byte for byte.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def rnd(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def te_names(naming: str, n: int) -> list[str]:
    if naming == "flybase":
        return [f"FBte{i:07d}" for i in range(1, n + 1)]
    stems = ("B1_Mm", "L1MdA_I", "Ty1", "IAPEz", "MERVL", "Ty3", "Gypsy", "Copia", "hAT", "Mariner")
    return [f"{stems[i % len(stems)]}_{i // len(stems) + 1}" for i in range(n)]


def family_names(naming: str, n: int) -> list[str]:
    if naming == "flybase":
        return [f"Fam{i}" for i in range(1, n + 1)]
    return [f"fam_{i}" for i in range(1, n + 1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", type=Path)
    parser.add_argument("--read-length", type=int, default=125)
    parser.add_argument("--fragment-mean", type=float, default=350.0)
    parser.add_argument("--fragment-sd", type=float, default=40.0)
    parser.add_argument("--depth", type=float, default=30.0)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--naming", choices=("flybase", "plain"), default="flybase")
    parser.add_argument("--contigs", default="chr2L,chr2R,chr3L,chr3R,chrX")
    parser.add_argument("--contig-length", type=int, default=200_000)
    parser.add_argument("--mito", default="chrM", help="name of the mitochondrial contig ('' for none)")
    parser.add_argument("--insertions", type=int, default=30)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    # TE library: 60 families of two related records each (the second record differs at ~3% of its bases).
    records, families = te_names(args.naming, 120), family_names(args.naming, 60)
    te_seq: dict[str, str] = {}
    te_family: dict[str, str] = {}
    first_base = ""
    for f, family in enumerate(families):
        base = rnd(rng, rng.randrange(900, 2600))
        if f == 0:
            first_base = base
        elif f >= 50:
            # the last ten families share 120 bp ends with the first, so their own end k-mers are not family-specific
            base = first_base[:120] + base[120:-120] + first_base[-120:]
        variant = "".join(rng.choice(BASES) if rng.random() < 0.03 else b for b in base)
        for name, seq in ((records[2 * f], base), (records[2 * f + 1], variant)):
            te_seq[name], te_family[name] = seq, family

    # Reference: random contigs with 4 resident TE copies each (diverged ~8%), annotated in BED6.
    contigs = [c for c in args.contigs.split(",") if c]
    reference: dict[str, str] = {}
    annotation: list[tuple[str, int, int, str, str]] = []
    for contig in contigs:
        seq = rnd(rng, args.contig_length)
        for slot in range(4):
            name = rng.choice(records[:40])
            copy = "".join(rng.choice(BASES) if rng.random() < 0.08 else b for b in te_seq[name])
            strand = rng.choice("+-")
            copy = copy if strand == "+" else rc(copy)
            at = (slot + 1) * args.contig_length // 5
            seq = seq[:at] + copy + seq[at + len(copy):]
            annotation.append((contig, at, at + len(copy), name, strand))
        # 14 near-identical copies of one record per contig: a high-copy family, as real genomes have (its k-mers then
        # occur more than 512 times in the reads, which the library layer's widest panel is meant for)
        slots = [start for where, start, _, _, _ in annotation if where == contig]
        spots = [at for at in range(2_000, args.contig_length - 3_000, 5_000) if all(abs(at - s) > 6_000 for s in slots)]
        for at in spots[:14]:
            name = records[2]
            copy = "".join(rng.choice(BASES) if rng.random() < 0.003 else b for b in te_seq[name])
            strand = rng.choice("+-")
            copy = copy if strand == "+" else rc(copy)
            seq = seq[:at] + copy + seq[at + len(copy):]
            annotation.append((contig, at, at + len(copy), name, strand))
        reference[contig] = seq
    if args.mito:
        reference[args.mito] = rnd(rng, 16_000)

    # Sample: non-reference insertions (5 bp TSD) on the nuclear contigs, away from the resident copies.
    sample = dict(reference)
    truth: list[tuple[str, int, str, str, str]] = []
    per_contig: dict[str, list[int]] = {c: [] for c in contigs}
    for i in range(args.insertions):
        contig = contigs[i % len(contigs)]
        while True:
            pos = rng.randrange(3_000, args.contig_length - 3_000)
            near_copy = any(c == contig and s - 2_000 < pos < e + 2_000 for c, s, e, _, _ in annotation)
            near_other = any(abs(pos - p) < 4_000 for p in per_contig[contig])
            if not near_copy and not near_other:
                break
        per_contig[contig].append(pos)
        name = rng.choice(records[40:])
        strand = rng.choice("+-")
        truth.append((contig, pos, name, te_family[name], strand))
    for contig in contigs:
        seq = reference[contig]
        for _, pos, name, _, strand in sorted((t for t in truth if t[0] == contig), key=lambda t: -t[1]):
            element = te_seq[name] if strand == "+" else rc(te_seq[name])
            seq = seq[:pos] + element + seq[pos - 5:]
        sample[contig] = seq

    # Paired reads: fragments from the sample genome, read 1 forward, read 2 reverse complement.
    length = args.read_length
    total = sum(len(s) for s in sample.values())
    pairs = int(total * args.depth / (2 * length))
    names = list(sample)
    weights = [len(sample[c]) for c in names]
    with (out / "sample_R1.fastq").open("w", newline="\n") as r1, (out / "sample_R2.fastq").open("w", newline="\n") as r2:
        for i in range(pairs):
            contig = rng.choices(names, weights)[0]
            seq = sample[contig]
            span = max(length, int(rng.gauss(args.fragment_mean, args.fragment_sd)))
            if len(seq) <= span:
                continue
            start = rng.randrange(0, len(seq) - span)
            fragment = seq[start:start + span]
            if rng.random() < 0.5:
                fragment = rc(fragment)
            mate1, mate2 = fragment[:length], rc(fragment[-length:])
            quality = "I" * length
            r1.write(f"@pair{i}/1\n{mate1}\n+\n{quality}\n")
            r2.write(f"@pair{i}/2\n{mate2}\n+\n{quality}\n")

    with (out / "reference.fa").open("w", newline="\n") as handle:
        for contig, seq in reference.items():
            handle.write(f">{contig}\n{wrap(seq)}\n")
    with (out / "te_library.fa").open("w", newline="\n") as handle:
        for name in records:
            handle.write(f">{name} name={te_family[name]}\n{wrap(te_seq[name])}\n")
    with (out / "family_map.tsv").open("w", newline="\n") as handle:
        if args.naming == "flybase":
            handle.write("FBte_ID\tFlyBase_Family\n")
        for name in records:
            handle.write(f"{name}\t{te_family[name]}\n")
    with (out / "repeat_annotation.bed").open("w", newline="\n") as handle:
        for contig, start, end, name, strand in sorted(annotation):
            handle.write(f"{contig}\t{start}\t{end}\t{name}\t0\t{strand}\n")
    with (out / "truth_insertions.tsv").open("w", newline="\n") as handle:
        handle.write("contig\tposition0\tte_record\tfamily\tstrand\n")
        for row in sorted(truth):
            handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
