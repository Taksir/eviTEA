"""``evitea build-repeat-map``: the repeated part of a reference, from the reference's alignment to itself.

A run without a TE library has no curated annotation of the reference's TE copies to read, and the resident-copy test and the
record guard need to know where old copies sit.  The map answers that from the sequence alone.  Windows of 500 bp (step 250) of
the chosen contigs are aligned to the whole reference with BLAST+ (``makeblastdb``, ``blastn``; both on PATH).  A base is
repeated when a hit of at least 100 aligned bases at 80% identity or more covers it, other than the window's own locus (same
contig, subject start within 300 bp of the query position).  Query sides and, on the chosen contigs, subject sides are marked
and merged.

Output: ``<output>/repeat_map.bed`` (contig, 0-based start, end), read by ``--repeat-map`` of ``library-free`` and
``build-family-records``; ``<output>/repeat_map_bed6.bed``, the same intervals as BED6 (name ``repeat``, score 0, strand
``.``), which ``prepare-reference --repeat-annotation`` accepts, so that a run without a library requires no curated annotation;
``<output>/repeat_map_summary.json`` holds the constants, the counts and the checksums.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from collections import defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path

from .errors import InputContractError, OperationalResourceError

__all__ = [
    "BED", "BED6", "build_repeat_map", "default_contigs", "intervals_from_hits", "merge_intervals", "read_lengths",
    "window_records",
]

WINDOW = 500
STEP = 250
MIN_IDENTITY = 80.0     # percent
MIN_ALIGNED = 100       # aligned bases of a hit
SELF_BP = 300           # a hit at the window's own locus: subject start within this many bases of the query position
EVALUE = "1e-5"
MAX_SUBJECTS = 10
MAX_HSPS = 10
HIT_FORMAT = "6 qseqid sseqid pident length qstart qend sstart send"
BED = "repeat_map.bed"
BED6 = "repeat_map_bed6.bed"
BED6_NAME = "repeat"
SUMMARY = "repeat_map_summary.json"
IDENTITY = "build_identity.json"
TOOLS = ("makeblastdb", "blastn")


def _fasta(path: Path) -> Iterator[tuple[str, str]]:
    name, buf = None, []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                fields = line[1:].split()
                if not fields:
                    raise InputContractError(f"reference {path} has a sequence with an empty name")
                name, buf = fields[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def read_lengths(reference_fasta: Path) -> dict[str, int]:
    """Return contig lengths by name, in reference order."""

    lengths: dict[str, int] = {}
    for name, seq in _fasta(reference_fasta):
        if name in lengths:
            raise InputContractError(f"reference {reference_fasta} has the contig name {name!r} twice")
        lengths[name] = len(seq)
    return lengths


def default_contigs(reference_fasta: Path) -> tuple[str, ...]:
    """Return the contigs ``library-free`` scans by default: those of at least 1 Mb, in reference order."""

    from .library_free.runtime import MIN_DEFAULT_CONTIG

    return tuple(name for name, n in read_lengths(reference_fasta).items() if n >= MIN_DEFAULT_CONTIG)


def window_records(reference_fasta: Path, contigs: Iterable[str]) -> Iterator[tuple[str, int, str]]:
    """Yield (contig, start, window sequence) of each 500 bp window (step 250) of the chosen contigs, in reference order.
    A window of more than half N is left out."""

    chosen = set(contigs)
    for name, seq in _fasta(reference_fasta):
        if name not in chosen:
            continue
        seq = seq.upper()
        for start in range(0, len(seq) - WINDOW + 1, STEP):
            window = seq[start:start + WINDOW]
            if window.count("N") > WINDOW // 2:
                continue
            yield name, start, window


def merge_intervals(intervals: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for contig, ivs in intervals.items():
        merged: list[list[int]] = []
        for a, b in sorted(ivs):
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        out[contig] = [(a, b) for a, b in merged]
    return out


def subject_name(subject: str, names: set[str] | None) -> str:
    """Return the reference's own name of a blastn subject id.  With ``-parse_seqids`` blastn prints an id that looks like a database
    accession with its wrappers (``ref|NW_001845059.1|``, ``gb|CM000001.1|``, ``gi|123|ref|NW_001845059.1|``) and a plain name as
    is or as ``lcl|name``.  The name is the last field of the id that is a contig of the reference; with ``names`` None only a
    ``lcl|`` prefix is removed.  An id with no such field is an error."""

    if names is None:
        return subject[4:] if subject.startswith("lcl|") else subject
    if subject in names:
        return subject
    for field in reversed(subject.split("|")):
        if field in names:
            return field
    raise InputContractError(f"blastn reported the subject {subject!r}, which is not a contig of the reference")


def intervals_from_hits(lines: Iterable[str], contigs: Iterable[str], known: Iterable[str] | None = None) -> dict[str, list[tuple[int, int]]]:
    """Return repeated intervals (0-based, end exclusive, not yet merged) from blastn lines in ``HIT_FORMAT``.  Query ids are
    ``<contig>|<window start>``.  ``known`` (all contig names of the reference) turns the subject ids back into the reference's names (see ``subject_name``)."""

    chosen = set(contigs)
    names = set(known) if known is not None else None
    iv: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for line in lines:
        f = line.split()
        if len(f) < 8:
            continue
        q, s, pident, length, qs, qe, ss, se = f[:8]
        s = subject_name(s, names)
        if float(pident) < MIN_IDENTITY or int(length) < MIN_ALIGNED:
            continue
        contig, window = q.rsplit("|", 1)
        window, qs, qe, ss, se = int(window), int(qs), int(qe), int(ss), int(se)
        lo, hi = min(ss, se), max(ss, se)
        if s == contig and abs(lo - (window + qs)) < SELF_BP:
            continue
        iv[contig].append((window + qs - 1, window + qe))
        if s in chosen:
            iv[s].append((lo - 1, hi))
    return iv


def write_bed(path: Path, intervals: dict[str, list[tuple[int, int]]]) -> tuple[int, int]:
    """Write the merged intervals, contigs sorted by name, and return (intervals, bases)."""

    n = bases = 0
    with Path(path).open("w", encoding="utf-8") as handle:
        for contig in sorted(intervals):
            for a, b in intervals[contig]:
                handle.write(f"{contig}\t{a}\t{b}\n")
                n += 1
                bases += b - a
    return n, bases


def write_bed6(path: Path, intervals: dict[str, list[tuple[int, int]]]) -> None:
    """Write the intervals of ``write_bed`` as BED6 (name ``repeat``, score 0, no strand): a repeat annotation for prepare-reference."""

    with Path(path).open("w", encoding="utf-8") as handle:
        for contig in sorted(intervals):
            for a, b in intervals[contig]:
                handle.write(f"{contig}\t{a}\t{b}\t{BED6_NAME}\t0\t.\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def _require_blast() -> None:
    missing = [t for t in TOOLS if shutil.which(t) is None]
    if missing:
        raise OperationalResourceError(
            "build-repeat-map requires BLAST+ (makeblastdb and blastn) on PATH; not found: " + ", ".join(missing))


def build_repeat_map(*, reference_fasta: Path, output: Path, contigs: Iterable[str] | None = None, parts: int | None = None,
                     progress=print) -> dict:
    """Build ``<output>/repeat_map.bed``.  Saved BLAST work in ``output`` is reused when the reference and the parameters are the
    same (the build resumes).  An existing BED is never replaced."""

    reference_fasta, output = Path(reference_fasta), Path(output)
    lengths = read_lengths(reference_fasta)
    chosen = tuple(contigs) if contigs is not None else default_contigs(reference_fasta)
    bad = [c for c in chosen if c not in lengths]
    if not chosen or bad or len(set(chosen)) != len(chosen):
        raise InputContractError("contigs must be distinct contigs of the reference"
                                 + (f"; not found: {bad}" if bad else "; none given and none of 1 Mb or more"))
    parts = max(1, int(parts if parts is not None else min(os.cpu_count() or 1, 8)))
    _require_blast()
    output.mkdir(parents=True, exist_ok=True)
    if (output / BED).exists():
        raise InputContractError(f"{output / BED} exists; use an output folder of its own")
    identity = {
        "reference_sha256": _sha256(reference_fasta), "contigs": list(chosen),
        "window": WINDOW, "step": STEP, "evalue": EVALUE, "max_subjects": MAX_SUBJECTS, "max_hsps": MAX_HSPS,
    }
    marker = output / IDENTITY
    if marker.exists() and json.loads(marker.read_text(encoding="utf-8")) != identity:
        raise InputContractError(f"{output} holds BLAST work of another reference or other parameters; use an output folder of its own")
    marker.write_text(json.dumps(identity, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    db = output / "db"
    if not (output / "db.done").exists():
        progress("repeat map: makeblastdb")
        subprocess.run(["makeblastdb", "-in", str(reference_fasta), "-dbtype", "nucl", "-out", str(db), "-parse_seqids"],
                       check=True, stdout=subprocess.DEVNULL)
        (output / "db.done").write_text("", encoding="utf-8")
    if not (output / "hits.done").exists():
        files = [(output / f"win_{i}.fa").open("w", encoding="utf-8") for i in range(parts)]
        n = 0
        for contig, start, window in window_records(reference_fasta, chosen):
            files[n % parts].write(f">{contig}|{start}\n{window}\n")
            n += 1
        for handle in files:
            handle.close()
        progress(f"repeat map: blastn of {n} windows in {parts} parts")
        procs = [subprocess.Popen(["blastn", "-task", "blastn", "-query", str(output / f"win_{i}.fa"), "-db", str(db),
                                   "-evalue", EVALUE, "-max_target_seqs", str(MAX_SUBJECTS), "-max_hsps", str(MAX_HSPS),
                                   "-num_threads", "1", "-outfmt", HIT_FORMAT, "-out", str(output / f"hits_{i}.tsv")])
                 for i in range(parts)]
        if any([p.wait() for p in procs]):
            raise OperationalResourceError("blastn failed")
        (output / "hits.done").write_text(f"windows {n} parts {parts}\n", encoding="utf-8")
    done = (output / "hits.done").read_text(encoding="utf-8").split()
    windows, saved_parts = int(done[1]), int(done[3])
    lines = []
    for i in range(saved_parts):
        lines.extend((output / f"hits_{i}.tsv").read_text(encoding="utf-8").splitlines())
    iv = merge_intervals(intervals_from_hits(lines, chosen, known=lengths))
    count, bases = write_bed(output / BED, iv)
    write_bed6(output / BED6, iv)
    summary = {
        "reference_bases": sum(lengths.values()), "contigs": list(chosen), "windows": windows,
        "constants": {"window": WINDOW, "step": STEP, "min_identity": MIN_IDENTITY, "min_aligned": MIN_ALIGNED, "self_bp": SELF_BP,
                      "evalue": EVALUE, "max_subjects": MAX_SUBJECTS, "max_hsps": MAX_HSPS},
        "intervals": count, "bases": bases, "reference_sha256": identity["reference_sha256"], "bed_sha256": _sha256(output / BED),
        "bed6_sha256": _sha256(output / BED6),
    }
    (output / SUMMARY).write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    progress(f"repeat map: intervals {count} bases {bases}")
    return summary
