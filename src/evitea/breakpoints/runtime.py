"""Real backends for breakpoint refinement: msBWT read counts and a minimap2 segment aligner.

minimap2 is reached through ``mappy`` (optional extra ``evitea[breakpoints]``).
mappy publishes no CPython 3.14 Windows wheel, so the switch currently needs
Linux (or WSL); without mappy it fails closed.
"""

from __future__ import annotations

import csv
import shutil
from collections.abc import Callable
from pathlib import Path

from ..backends.msbwt import load_msbwt
from ..callset import CallsetState
from ..errors import InputContractError, OperationalResourceError
from . import walker as W
from . import MIN_WIDTH
from .refine import BreakpointRow, RefinementResult, breakpoint_row, refine_native, walk_native

__all__ = ["BreakpointRefiner", "read_fasta"]


def read_fasta(path: Path) -> dict[str, str]:
    """name (first word of the header) -> sequence, in file order."""

    records: dict[str, list[str]] = {}
    name = None
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise InputContractError(f"duplicate FASTA record {name!r} in {path}")
                records[name] = []
            elif name is not None:
                records[name].append(line)
    return {key: "".join(parts) for key, parts in records.items()}


def _family_map(path: Path) -> dict[str, str]:
    """FlyBase-style table (FBte_ID, FlyBase_Family) or two columns record<TAB>family; labels kept as written."""

    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if rows and "FBte_ID" in rows[0] and "FlyBase_Family" in rows[0]:
        i, j = rows[0].index("FBte_ID"), rows[0].index("FlyBase_Family")
        return {r[i]: r[j] for r in rows[1:] if len(r) > max(i, j)}
    return {r[0]: r[1] for r in rows if len(r) >= 2}


def _rm_records(path: Path):
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            ch, s, e, i = line.rstrip("\r\n").split("\t")[:4]
            yield ch, int(s), int(e), i


def _mappy_aligner(te_fasta: Path, reference_fasta: Path, work_dir: Path):
    try:
        import mappy
    except ImportError as exc:
        raise OperationalResourceError(
            "breakpoint refinement needs mappy (minimap2); install evitea[breakpoints] on Linux or WSL"
        ) from exc
    work_dir.mkdir(parents=True, exist_ok=True)
    combined = work_dir / "breakpoint_aligner_library.fa"
    with combined.open("wb") as out:  # TE library first, then the host reference (byte concatenation)
        for source in (te_fasta, reference_fasta):
            with Path(source).open("rb") as handle:
                shutil.copyfileobj(handle, out)
    aligner = mappy.Aligner(str(combined), preset="sr")
    if not aligner:
        raise OperationalResourceError("minimap2 failed to index the breakpoint aligner library")
    cons = {name for name, _, _ in mappy.fastx_read(str(te_fasta))}
    return aligner, cons, getattr(mappy, "__version__", "unknown")


class BreakpointRefiner:
    """Callable ``CallsetState -> tuple[BreakpointRow, ...]`` over one sample msBWT and the prepared assets."""

    def __init__(
        self,
        *,
        msbwt_index: Path,
        reference_fasta: Path,
        repeat_annotation: Path,
        te_fasta: Path,
        family_map: Path,
        work_dir: Path,
        progress: Callable[[str], None] | None = None,
        pin_one_sided: bool = False,
        call_annotations: bool = False,
    ) -> None:
        self.pin_one_sided = bool(pin_one_sided)
        self.call_annotations = bool(call_annotations)
        self.msbwt_index = Path(msbwt_index)
        self.reference_fasta = Path(reference_fasta)
        self.repeat_annotation = Path(repeat_annotation)
        self.te_fasta = Path(te_fasta)
        self.family_map = Path(family_map)
        self.work_dir = Path(work_dir)
        self.progress = progress
        self.aligner_version: str | None = None

    def __call__(self, state: CallsetState) -> tuple[BreakpointRow, ...] | RefinementResult:
        """Rows only with both extras off (the switch-only behaviour); a RefinementResult otherwise."""

        aligner, cons, self.aligner_version = _mappy_aligner(self.te_fasta, self.reference_fasta, self.work_dir)
        rm_index = W.build_rm_index(_rm_records(self.repeat_annotation))
        classifier = W.SegmentClassifier(aligner, cons, _family_map(self.family_map), rm_index)
        bwt = load_msbwt(self.msbwt_index)
        counter = W.KmerCounter(bwt.countOccurrencesOfSeq)
        reference = read_fasta(self.reference_fasta)
        if not (self.pin_one_sided or self.call_annotations):
            return refine_native(state.native, reference=reference, counter=counter, classifier=classifier,
                                 rm_index=rm_index, progress=self.progress)
        from .annotate import annotate_walk
        from .pinning import pin_rows

        walks = walk_native(state.native, reference=reference, counter=counter, classifier=classifier,
                            rm_index=rm_index, include_narrow=self.call_annotations, progress=self.progress)
        wide = tuple(walk for walk in walks if walk.e0 - walk.s0 > MIN_WIDTH)
        rows = tuple(breakpoint_row(walk) for walk in wide)
        pins = pin_rows(rows, wide) if self.pin_one_sided else None
        annotations = (tuple(annotate_walk(walk, reference[walk.chrom]) for walk in walks)
                       if self.call_annotations else None)
        return RefinementResult(rows, pins, annotations)
