"""Checks made at the start of a run, so that unsupported input stops with a plain message before any work.

* Paired-end reads: read 1 and read 2 are two different files whose first records pair up (same read name), or a prepared
  sample / msBWT index that holds both mates of every pair.
* One read length: every read has the same length (the first ``READ_LENGTH_SAMPLE`` reads of each FASTQ are checked here; the
  sample builder checks every read).
* Contig names: the contigs a run names (``--primary-contigs``) are contigs of the reference.
* TE record ids: the TE library's records are named in the family map.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

from .errors import InputContractError, RawInputValidationError

__all__ = [
    "READ_LENGTH_SAMPLE",
    "check_contig_names",
    "check_paired_reads",
    "check_te_records",
    "check_msbwt_index",
]

READ_LENGTH_SAMPLE = 10_000


def _first_records(path: Path, count: int, mate) -> list[tuple[str, int]]:
    """(pair name, length) of the first ``count`` reads, read with the sample builder's FASTQ reader."""

    from .inputs import _iter_fastq_records

    out = []
    for record in _iter_fastq_records(path, mate):
        out.append((record.pair_name, len(record.sequence)))
        if len(out) >= count:
            break
    return out


def check_paired_reads(read1: Path, read2: Path, *, sample: int = READ_LENGTH_SAMPLE) -> int:
    """Refuse single-end or mismatched input and reads of several lengths; returns the read length seen."""

    read1, read2 = Path(read1), Path(read2)
    if read1.resolve() == read2.resolve():
        raise RawInputValidationError(
            f"read 1 and read 2 are the same file ({read1}); EviTEA needs paired-end reads in two files"
        )
    from .contracts import ReadMate

    first, second = _first_records(read1, sample, ReadMate.READ1), _first_records(read2, sample, ReadMate.READ2)
    if not first or not second:
        raise RawInputValidationError("a FASTQ file holds no reads; EviTEA needs paired-end reads in two files")
    if first[0][0] != second[0][0]:
        raise RawInputValidationError(
            f"the first reads of the two files do not pair up ({first[0][0]!r} and {second[0][0]!r}); EviTEA needs "
            "paired-end reads with read 1 and read 2 of each pair at the same place in the two files"
        )
    lengths = sorted({length for _, length in first + second})
    if len(lengths) > 1:
        shown = ", ".join(str(x) for x in lengths[:5]) + (" ..." if len(lengths) > 5 else "")
        raise RawInputValidationError(
            f"the reads have more than one length ({shown} bp); EviTEA needs reads of one length: trim them to one "
            "length first (for example with cutadapt --length, keeping only pairs of the full length)"
        )
    return lengths[0]


def check_msbwt_index(index_dir: Path) -> None:
    """An msBWT index of paired reads (even record count, both mates) of one length."""

    import numpy as np

    from .sample_geometry import read_length_from_index

    about = np.load(Path(index_dir) / "about.npy", mmap_mode="r", allow_pickle=False)
    names = getattr(about.dtype, "names", None)
    if about.ndim != 1 or not names or about.shape[0] % 2:
        raise InputContractError(
            f"the msBWT index {index_dir} does not hold paired-end reads (an even number of reads, read 1 and read 2 of "
            "every pair)"
        )
    if set(np.unique(about[names[0]]).tolist()) != {0, 1}:
        raise InputContractError(f"the msBWT index {index_dir} was not built from two FASTQ files (read 1 and read 2)")
    read_length_from_index(index_dir)


def check_contig_names(names: Iterable[str], reference_contigs: Iterable[str], option: str) -> None:
    known = set(reference_contigs)
    missing = [name for name in names if name not in known]
    if missing:
        hint = [f"chr{m}" for m in missing if f"chr{m}" in known]
        hint += [m[3:] for m in missing if m.startswith("chr") and m[3:] in known]
        raise InputContractError(
            f"{option} names contigs that are not in the reference: {missing[:10]}"
            + (f" (the reference names them {hint[:10]})" if hint else "")
        )


def check_te_records(records: Iterable[str], families: Mapping[str, str]) -> None:
    """At least one library record must be named in the family map."""

    records = list(records)
    named = [r for r in records if r in families]
    if not named:
        raise InputContractError(
            "no TE library record is named in the family map: the map's record ids must be the first word of the TE "
            f"FASTA headers (library: {records[:3]}, map: {list(families)[:3]})"
        )
