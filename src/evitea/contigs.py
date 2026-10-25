"""Contig names of the reference: bare names, the mitochondrial contig and the primary contigs.

Bare names
    An annotation may name a contig without the ``chr`` prefix the reference uses (``2L`` for ``chr2L``, ``1`` for ``chr1``).
    ``match_contig_name`` maps a bare name to ``chr`` + name only when that contig exists in the reference and the bare name
    itself does not; every other name is kept as it is.

Mitochondrial contig
    The library-free layer's segment firewall tests segments against the mitochondrial sequence (nuclear copies of
    mitochondrial DNA).  ``--mito-contig`` names that contig; without it the first of ``MITO_CONTIG_NAMES`` (``chrM``, ``MT``,
    ``M``, ``chrMT``) found in the reference is used, and with none of them the mitochondrial test is off.

Primary contigs
    The full-memory proposals of the library layer look only at the primary contigs.  ``--primary-contigs`` lists them (a comma
    list, or a file with one name per line).  Without it they are taken from the reference by this rule:

    * the contig is at least ``PRIMARY_MIN_FRACTION`` (10 %) as long as the longest contig, so the rule scales with the genome:
      it keeps every yeast chromosome (the shortest, I, is about 15 % of the longest) and every mouse or human chromosome, and
      leaves out unplaced scaffolds (at most a few percent of the longest chromosome in those assemblies);
    * the contig is not an organelle genome (``ORGANELLE_CONTIG_NAMES``: mitochondrion and chloroplast names);
    * the name does not mark an unplaced, unlocalised or alternative sequence (``chrUn...``, ``Un...``, or a name with
      ``_random``, ``_alt``, ``_fix``, ``_hap`` or ``_decoy``).

    Contigs are listed in reference order.  On dm6 the rule gives chr2L, chr2R, chr3L, chr3R, chrX and chrY (chr4 is 4 % of
    chr3R); the fly runs of the paper pass ``--primary-contigs chr2L,chr2R,chr3L,chr3R,chrX``, the list the layer was built
    with (``DEVELOPMENT_PRIMARY_CONTIGS``).
"""

from __future__ import annotations

import re
from collections.abc import Collection, Iterable, Mapping
from pathlib import Path

from .errors import InputContractError

__all__ = [
    "DEVELOPMENT_PRIMARY_CONTIGS",
    "MITO_CONTIG_NAMES",
    "ORGANELLE_CONTIG_NAMES",
    "PRIMARY_MIN_FRACTION",
    "default_primary_contigs",
    "find_mito_contig",
    "match_contig_name",
    "parse_contig_list",
    "resolve_primary_contigs",
]

DEVELOPMENT_PRIMARY_CONTIGS = ("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
MITO_CONTIG_NAMES = ("chrM", "MT", "M", "chrMT")
ORGANELLE_CONTIG_NAMES = frozenset(
    name.lower() for name in (*MITO_CONTIG_NAMES, "mito", "mitochondrion", "chrMito", "ChrC", "chrPt", "Pt", "chloroplast")
)
PRIMARY_MIN_FRACTION = 0.10
_UNPLACED = re.compile(r"^(chr)?Un([_\W]|$)|_(random|alt|fix|hap\d*|decoy)\b", re.IGNORECASE)


def match_contig_name(name: str, contigs: Collection[str]) -> str:
    """``name``, or ``chr`` + ``name`` when only that one is a contig of the reference (``contigs``)."""

    if contigs and name not in contigs and "chr" + name in contigs:
        return "chr" + name
    return name


def find_mito_contig(contigs: Iterable[str], option: str | None = None) -> str | None:
    """The mitochondrial contig: ``option`` (which must exist), else the first of ``MITO_CONTIG_NAMES`` present, else None."""

    names = set(contigs)
    if option:
        if option not in names:
            raise InputContractError(f"--mito-contig {option!r} is not a contig of the reference")
        return option
    return next((name for name in MITO_CONTIG_NAMES if name in names), None)


def default_primary_contigs(lengths: Mapping[str, int]) -> tuple[str, ...]:
    """The primary contigs of a reference (``{name: length}`` in reference order) by the rule in the module notes."""

    if not lengths:
        return ()
    floor = PRIMARY_MIN_FRACTION * max(lengths.values())
    return tuple(
        name for name, length in lengths.items()
        if length >= floor and name.lower() not in ORGANELLE_CONTIG_NAMES and not _UNPLACED.search(name)
    )


def parse_contig_list(text: str) -> tuple[str, ...]:
    """A comma list of names, or the path of a file with one name per line (``#`` comments allowed)."""

    path = Path(text).expanduser()
    if "," not in text and path.is_file():
        lines = path.read_text(encoding="utf-8").splitlines()
        names = [line.split("#", 1)[0].strip() for line in lines]
    else:
        names = [part.strip() for part in text.split(",")]
    names = [name for name in names if name]
    if not names or len(set(names)) != len(names):
        raise InputContractError(f"contig list must name distinct contigs: {text!r}")
    return tuple(names)


def resolve_primary_contigs(option: Iterable[str] | None, lengths: Mapping[str, int]) -> tuple[tuple[str, ...], str]:
    """(primary contigs, source): the option's names (which must be reference contigs) or the reference rule."""

    if option is not None:
        names = tuple(option)
        missing = [name for name in names if name not in lengths]
        if missing:
            raise InputContractError(f"--primary-contigs names contigs that are not in the reference: {missing}")
        return names, "option"
    return default_primary_contigs(lengths), "reference-rule"
