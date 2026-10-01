"""Frozen external-aligner invocations and SAM readers for the incumbent path.

Bowtie2 seeds its per-read pseudo-random generator from the read name,
sequence and qualities, so query names and option lists are part of the
frozen semantics and are fixed here.
"""

from __future__ import annotations

import re
import subprocess
from collections import defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path

from ..errors import BackendCompatibilityError

__all__ = [
    "write_fasta",
    "run_bowtie2",
    "run_bowtie2_build",
    "parse_sam_records",
    "cigar_ref_span",
    "HOST_GLOBAL_OPTIONS",
    "EXACTISH_OPTIONS",
    "TE_LIBRARY_OPTIONS",
]

HOST_GLOBAL_OPTIONS = ("-f", "--very-sensitive", "--end-to-end", "--no-unal")
EXACTISH_OPTIONS = ("-f", "--very-sensitive", "--end-to-end", "--no-unal", "-N", "0", "-L", "20")
TE_LIBRARY_OPTIONS = (
    "-f", "-a", "--end-to-end", "--very-sensitive", "--no-unal",
    "-N", "1", "-L", "12", "-D", "30", "-R", "5", "--score-min", "L,-70,0",
)
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with Path(path).open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def run_bowtie2(executable: str | Path, arguments: Sequence[str]) -> None:
    command = [str(executable), *[str(value) for value in arguments]]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError as exc:
        raise BackendCompatibilityError(f"failed to start Bowtie2: {exc}") from exc
    if completed.returncode != 0:
        raise BackendCompatibilityError(
            f"Bowtie2 exited with code {completed.returncode}: {completed.stderr[-2000:].strip()}"
        )


def run_bowtie2_build(executable: str | Path, fasta: Path, prefix: Path, extra: Sequence[str] = ()) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [str(executable), *extra, str(fasta), str(prefix)]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise BackendCompatibilityError(
            f"bowtie2-build exited with code {completed.returncode}: {completed.stderr[-2000:].strip()}"
        )


def cigar_ref_span(cigar: str, fallback: int) -> int:
    if not cigar or cigar == "*":
        return fallback
    span = sum(int(count) for count, op in _CIGAR.findall(cigar) if op in {"M", "D", "N", "=", "X"})
    return span or fallback


def parse_sam_records(path: Path) -> dict[str, list[dict[str, object]]]:
    """Mapped SAM rows by query name (unmapped rows skipped)."""

    out: dict[str, list[dict[str, object]]] = defaultdict(list)
    with Path(path).open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 0x4:
                continue
            tags: dict[str, str] = {}
            for token in fields[11:]:
                parts = token.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            out[fields[0]].append({
                "flag": flag,
                "rname": fields[2],
                "pos0": int(fields[3]) - 1,
                "mapq": int(fields[4]),
                "cigar": fields[5],
                "seq": fields[9],
                "nm": int(tags.get("NM", "999999")),
                "rev": bool(flag & 0x10),
                "secondary": bool(flag & 0x100),
                "supplementary": bool(flag & 0x800),
            })
    return out
