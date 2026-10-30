"""Command-line parser; the scientific thresholds are fixed and have no option."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from collections.abc import Callable, Sequence
from importlib import metadata
from pathlib import Path

from ._version import __version__
from .assets.reference import load_prepared_reference, prepare_reference, verify_reference_asset
from .assets.te import load_prepared_te_assets, prepare_te_assets, verify_te_assets
from .backends.exact_query import load_prepared_sample, verify_prepared_sample
from .backends.msbwt import build_prepared_msbwt_sample
from .config import (
    EviTEAConfig,
    InputConfig,
    OutputConfig,
    PreparedAssetConfig,
    RefinementConfig,
    RuntimeConfig,
    load_config,
    merge_config,
)
from .constants import ARCHITECTURE_ID
from .contigs import parse_contig_list
from .preflight import check_contig_names, check_msbwt_index, check_paired_reads, check_te_records
from .sample_geometry import validate_read_length
from .errors import (
    ArchitectureInvariantError,
    AssetCompatibilityError,
    BackendCompatibilityError,
    EviTEAError,
    InputContractError,
    OperationalResourceError,
    OutputFinalizationError,
    OutputStateError,
    PairSynchronizationError,
    RawInputValidationError,
    ResumeStateError,
)
from .inputs import (
    is_flybase_family_table,
    load_family_map,
    load_repeat_annotation,
    validate_reference_fasta,
    validate_te_fasta,
)
from .pipeline import (
    BREAKPOINTS_OUTPUT,
    CANONICAL_OUTPUTS,
    OPTIONAL_OUTPUTS,
    PreparedAssets,
    RawSampleBuilder,
    RunResult,
    derive_sample_id,
    prepare_assets,
    run_sample,
)
from .provenance import (
    RUN_MANIFEST_SCHEMA,
    RunManifestContext,
    file_identity,
    hash_file,
)

CallEngine = Callable[..., object]


def _resolved_cli_path(path: Path | None) -> Path | None:
    return path.expanduser().resolve() if path is not None else None


def _distribution_version(*names: str) -> str | None:
    for name in names:
        try:
            return metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    return None


def _msbwt_distribution() -> tuple[str, str | None]:
    """Return (name, version) of the installed msBWT package: pymsbwt first, then the other package names."""

    for name in ("pymsbwt", "msbwt-modern3", "msbwt"):
        version = _distribution_version(name)
        if version is not None:
            return name, version
    return "pymsbwt", None


def _executable_version(executable: str | Path) -> str | None:
    resolved = shutil.which(str(executable))
    if resolved is None:
        return None
    try:
        completed = subprocess.run(
            (resolved, "--version"),
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    lines = (completed.stdout or completed.stderr).splitlines()
    return lines[0].strip() if completed.returncode == 0 and lines else None


def _source_tree_sha256() -> str:
    """sha256 over the relative path and the sha256 of every .py file of this package, sorted by path: which code ran, also where the
    container has no git. ``tools/compare_scan.py treehash`` computes the same value for a checkout."""

    package = Path(__file__).resolve().parent
    digest = hashlib.sha256()
    for path in sorted(p for p in package.rglob("*.py") if "__pycache__" not in p.parts):
        digest.update(path.relative_to(package).as_posix().encode() + b"\0" + hashlib.sha256(path.read_bytes()).hexdigest().encode() + b"\0")
    return digest.hexdigest()


def _source_identity() -> dict[str, object]:
    identity: dict[str, object] = {"package_version": __version__, "source_tree_sha256": _source_tree_sha256()}
    root = Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.run(
            ("git", "-C", str(root), "rev-parse", "HEAD"),
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ("git", "-C", str(root), "status", "--porcelain"),
                check=True,
                capture_output=True,
                text=True,
                timeout=10,
            ).stdout
        )
    except (OSError, subprocess.SubprocessError):
        return identity
    identity.update({"git_commit": commit, "git_dirty": dirty})
    return identity


def _add_engineering_options(
    parser: argparse.ArgumentParser, *, use_defaults: bool = True
) -> None:
    default_threads = 0 if use_defaults else None
    default_log_level = "INFO" if use_defaults else None
    parser.add_argument("--threads", type=int, default=default_threads)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument(
        "--log-level",
        choices=("ERROR", "WARNING", "INFO", "DEBUG"),
        default=default_log_level,
    )


def _add_breakpoint_extras(parser: argparse.ArgumentParser, *, use_defaults: bool = False) -> None:
    default = False if use_defaults else None
    parser.add_argument(
        "--pin-one-sided",
        action="store_true",
        default=default,
        help="with breakpoint refinement: report the one-sided TE junction of a call left wide "
        "(breakpoints.tsv pin_side / pin_pos; the call itself is unchanged)",
    )
    parser.add_argument(
        "--call-annotations",
        action="store_true",
        default=default,
        help="with breakpoint refinement: write call_annotations.tsv (zygosity label and TSD letters; "
        "output only)",
    )


KEEP_RESIDENT_COPIES_HELP = (
    "keep the calls inside a reference copy of their own TE family (default: they are left out of the final call tables "
    "and listed in resident_copy_dropped.tsv; use this switch to reproduce a run without that filter)"
)


def _add_keep_resident_copies(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--keep-resident-copies", action="store_true", help=KEEP_RESIDENT_COPIES_HELP)


def _add_depth_cap(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--depth-cap",
        default="none",
        metavar="N|none",
        help="absolute upper bound of the local-depth window of the breakpoint walk and the library-free layer "
        "(default none: only the lower bound of 10 reads applies; --depth-cap 100 caps it at 100)",
    )


def _apply_depth_cap(args: argparse.Namespace) -> None:
    from . import depth_window

    try:
        depth_window.set_depth_cap(args.depth_cap)
    except ValueError as exc:
        raise InputContractError(str(exc)) from exc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="evitea")
    commands = parser.add_subparsers(dest="command", required=True)

    reference = commands.add_parser("prepare-reference")
    reference.add_argument("--reference", type=Path, required=True)
    reference.add_argument("--output", type=Path, required=True)
    reference.add_argument("--repeat-annotation", type=Path)
    reference.add_argument("--bowtie2-build", default="bowtie2-build")
    reference.add_argument(
        "--bowtie2-index",
        type=Path,
        help="adopt an existing Bowtie2 index prefix of this reference instead of building one",
    )
    _add_engineering_options(reference)

    te = commands.add_parser("prepare-te")
    te.add_argument("--te-library", type=Path, required=True)
    te.add_argument(
        "--family-map",
        type=Path,
        required=True,
        help="TE record to family map: a FlyBase table (FBte_ID, FlyBase_Family columns) or two tab-separated "
        "columns record<TAB>family with any record names (the first word of each TE FASTA header)",
    )
    te.add_argument("--g1a-archive", type=Path, help=argparse.SUPPRESS)
    te.add_argument("--family-crosswalk", type=Path, help=argparse.SUPPRESS)
    te.add_argument(
        "--terminal-codebook",
        type=Path,
        help="final safe terminal 21-mer codebook (with --internal-codebook)",
    )
    te.add_argument(
        "--internal-codebook",
        type=Path,
        help="internal-only safe 21-mer codebook (with --terminal-codebook)",
    )
    te.add_argument(
        "--build-languages",
        action="store_true",
        help="build the terminal and internal 21-mer codebooks from the TE library, family map, reference and "
        "RepeatMasker annotation (evitea.languages)",
    )
    te.add_argument(
        "--read-length",
        type=int,
        default=125,
        help="with --build-languages: read length (bp) of the block column of the internal 21-mer audit (default 125)",
    )
    te.add_argument("--reference-assets", type=Path, required=True)
    te.add_argument("--output", type=Path, required=True)
    _add_engineering_options(te)

    sample = commands.add_parser("prepare-sample")
    sample.add_argument("--read1", type=Path, required=True)
    sample.add_argument("--read2", type=Path, required=True)
    sample.add_argument("--output", type=Path, required=True)
    sample.add_argument("--sample-id")
    sample.add_argument("--exact-query-backend", choices=("msbwt",), default="msbwt")
    _add_engineering_options(sample)

    call = commands.add_parser("call")
    sample_mode = call.add_mutually_exclusive_group()
    sample_mode.add_argument("--sample-assets", type=Path)
    sample_mode.add_argument("--read1", type=Path)
    sample_mode.add_argument(
        "--msbwt-index",
        type=Path,
        help="existing msBWT index directory (with about.npy); requires --sample-id",
    )
    call.add_argument("--read2", type=Path)
    call.add_argument("--reference-assets", type=Path)
    call.add_argument("--te-assets", type=Path)
    call.add_argument("--output", type=Path)
    call.add_argument("--sample-id")
    call.add_argument("--config", type=Path)
    call.add_argument("--resume", action="store_true", default=None)
    call.add_argument("--keep-intermediates", action="store_true", default=None)
    call.add_argument("--exact-query-backend", choices=("msbwt",))
    call.add_argument(
        "--profile",
        help="'production' (default), or AKA017 or ISO1 for the read names and event labels of those samples",
    )
    call.add_argument("--bowtie2", type=Path)
    call.add_argument("--bowtie2-build", type=Path)
    call.add_argument(
        "--edlib-python",
        type=Path,
        help="interpreter with edlib for the fragment decoder, used only if this one lacks it",
    )
    call.add_argument("--workers", type=int, help="locus-batch worker processes")
    call.add_argument(
        "--refine-breakpoints",
        action="store_true",
        default=None,
        help="narrow wide calls to the TE junction found in the reads (off by default; requires mappy)",
    )
    call.add_argument(
        "--complete-routing",
        action="store_true",
        default=None,
        help="first stage: route reads through other TE anchor k-mers (not low-complexity, one family, at most 1024 "
        "occurrences in the sample); keep an extra pair only when its TE-side read has no alignment within q99 of the "
        "host-side read (off by default)",
    )
    _add_breakpoint_extras(call)
    _add_depth_cap(call)
    _add_keep_resident_copies(call)
    call.add_argument(
        "--max-rescue-interval-mass",
        type=int,
        help="wide-anchor work limit in msBWT interval rows (default 1500000)",
    )
    call.add_argument(
        "--read-length",
        type=int,
        help="read length in bp (default: taken from the sample; the value used is written to run_manifest.json)",
    )
    call.add_argument(
        "--fragment-q95",
        type=float,
        help="95th percentile of the fragment size in bp (default: estimated from the sample's properly paired "
        "reads, 494 when too few; the value used is written to run_manifest.json)",
    )
    call.add_argument(
        "--fragment-q99",
        type=int,
        help="99th percentile of the fragment size in bp, also Bowtie2's maximum fragment length (default: "
        "estimated from the sample, 805 when too few pairs)",
    )
    call.add_argument(
        "--primary-contigs",
        metavar="LIST|FILE",
        help="contigs of the full-memory proposals: a comma list or a file with one name per line (default: "
        "contigs at least 10%% as long as the longest, without organelle, unplaced or alternative contigs; the fly "
        "runs pass chr2L,chr2R,chr3L,chr3R,chrX)",
    )
    _add_engineering_options(call, use_defaults=False)

    refine = commands.add_parser(
        "refine-breakpoints",
        help="apply breakpoint refinement to an existing call output directory",
    )
    refine.add_argument("--run-dir", type=Path, required=True, help="completed `evitea call` output")
    refine_sample = refine.add_mutually_exclusive_group(required=True)
    refine_sample.add_argument("--msbwt-index", type=Path)
    refine_sample.add_argument("--sample-assets", type=Path)
    refine.add_argument("--reference-assets", type=Path, required=True)
    refine.add_argument("--te-assets", type=Path, required=True)
    refine.add_argument("--output", type=Path, required=True)
    refine.add_argument("--work-dir", type=Path)
    _add_breakpoint_extras(refine, use_defaults=True)
    _add_depth_cap(refine)
    _add_keep_resident_copies(refine)

    lf = commands.add_parser(
        "library-free",
        help="library-free detector on an existing call output: family calls with no EviTEA call "
        "within 500 bp are added (Linux/WSL; requires mappy and parasail)",
    )
    lf.add_argument("--run-dir", type=Path, required=True, help="completed `evitea call` output")
    lf_sample = lf.add_mutually_exclusive_group()
    lf_sample.add_argument("--msbwt-index", type=Path)
    lf_sample.add_argument("--sample-assets", type=Path)
    lf.add_argument("--reference-assets", type=Path, required=True)
    lf.add_argument("--output", type=Path, required=True)
    lf.add_argument("--work-dir", type=Path)
    lf.add_argument("--workers", type=int, default=1, metavar="N",
                    help="processes that scan a contig (default 1; 0 = every CPU of the job). A contig is cut into tiles of 250 kb "
                    "that the processes walk, with the same places as one process; the scan pieces and every output table are the "
                    "same for any N")
    lf.add_argument("--counter", choices=("plain", "fast"), default="fast",
                    help="how the read index is searched (default fast). fast keeps and extends the search intervals of the "
                    "words it has counted, so it makes fewer index searches for the same counts; it checks 2,000 words against the "
                    "plain counts before the scan and stops on any difference; plain is the search of the basic scan, one word at a time")
    lf.add_argument("--reference-cache", type=Path, metavar="DIR",
                    help="folder that keeps the reference word tables (multi-copy words, walk-out words) under a key of the "
                    "reference and the detection code, so that a later job loads them instead of building them (default: "
                    "refcache/ in the work directory; one folder can serve several work directories)")
    lf.add_argument("--contigs", help="comma-separated contigs to scan, in this order (default: the primary chromosomes and the "
                    "contigs >= 1 Mb, reference order)")
    lf.add_argument("--rule-c", nargs="?", choices=("on", "off"), const="on", default="on",
                    help="the family-identity filter: also drop family calls whose members have a median pairwise "
                    "identity of <= 98.31%% (default on; a bare --rule-c means on)")
    lf.add_argument("--rule-c-tsd-skip", choices=("on", "off"), default="on",
                    help="TSD skip of the family-identity filter: a family that the filter would drop is kept when at least "
                    "--rule-c-tsd-min of its trusted members (calls near a library call, or passing the one-site check) have "
                    "a TSD of 2 to 20 bp at a share of at least --rule-c-tsd-share (default on)")
    lf.add_argument("--rule-c-tsd-share", type=float, default=0.3, metavar="SHARE",
                    help="share of the trusted members of a family with a TSD at which the TSD skip keeps the family "
                    "(default 0.3: families the filter should drop sit at 0.11 to 0.16, real families at 0.39 or more on the "
                    "training runs)")
    lf.add_argument("--rule-c-tsd-min", type=int, default=5, metavar="N",
                    help="trusted members a family needs for the TSD skip to judge it (default 5)")
    lf.add_argument("--mito-contig",
                    help="the reference's mitochondrial contig; the segment filter drops inserted segments that match it (default: "
                    "first of chrM, MT, M, chrMT found; none: this test is off)")
    lf.add_argument("--adapter-test", choices=("on", "off"), default="on",
                    help="set aside the places whose inserted segment holds a sequencing-adapter core (TruSeq or Nextera, also "
                    "with the junction up to 4 letters inside the adapter) before the grouping, so they form no family "
                    "(default on; not run with --detection-dir, whose grouping stands)")
    lf.add_argument("--single-copy-contigs", default="auto", metavar="auto|all|none|LIST",
                    help="contigs the sample carries in one copy (the X and Y of a male, every contig of a haploid genome): there "
                    "a library-free call that would be added is dropped when it comes from the haplotype scan or has the "
                    "heterozygous genotype (reference-path fraction above 0.15). auto: contigs at no more than 0.625 of the "
                    "genome depth whose library-confirmed calls (at least 10) are at most 25%% heterozygous; all: every contig; "
                    "none: off; or comma-separated contig names (default auto)")
    lf.add_argument("--het-guard", choices=("on", "off"), default="on",
                    help="local heterozygous guard: a call that would be added with the reference-path fraction above 0.15 is dropped "
                    "when fewer than a quarter of the 10 confirmed library-free places nearest to it on the contig (the call "
                    "itself counts) are heterozygous-looking; a contig with fewer than 5 confirmed places uses the share of the "
                    "sample, a sample with fewer than 10 has no guard (default on). Turn it off for a screen for new insertions "
                    "in an inbred line: a recent private insertion there is heterozygous, and the guard drops it among homozygous "
                    "neighbours")
    lf.add_argument("--te-assets", type=Path,
                    help="prepared TE assets of the library run (the TE record to family map): calls.combined.tsv and "
                    "calls.combined.vcf then leave out the library calls inside a reference copy of their own family "
                    "(resident_copy_dropped.tsv); the detection and the added calls read all library calls. Required when "
                    "the run has library calls, unless --keep-resident-copies is given")
    _add_keep_resident_copies(lf)
    lf.add_argument("--detection-dir", type=Path,
                    help="reuse library_free_places.tsv and library_free_segments.fa from an earlier detection")
    _add_depth_cap(lf)
    from .library_free.options import add_arguments as _add_lf_options

    _add_lf_options(lf)

    gt = commands.add_parser(
        "genotype",
        help="hom / het genotype of every call of a library-free output (genotypes.tsv); the call files are only read",
    )
    gt.add_argument("--library-free-dir", type=Path, required=True,
                    help="completed `evitea library-free` output (calls.combined.tsv, calls.library_free.tsv, library_free_family_calls.tsv)")
    gt_sample = gt.add_mutually_exclusive_group()
    gt_sample.add_argument("--msbwt-index", type=Path)
    gt_sample.add_argument("--sample-assets", type=Path)
    gt.add_argument("--reference-assets", type=Path,
                    help="prepared reference of the run; needed with the msBWT index for the library calls (Step B)")
    gt.add_argument("--output", type=Path, help="folder for genotypes.tsv and genotype_summary.json (default: the library-free folder)")
    gt.add_argument("--counter", choices=("plain", "fast"), default="fast")

    from .family_records import run as _fr

    fr_build = commands.add_parser(
        "build-family-records",
        help="library records for library-free families missing from the library, with their reference copies and "
        "TE-end and TE-interior k-mer codebooks (Linux/WSL; requires mappy)",
    )
    _fr.add_arguments_build(fr_build)
    fr_add = commands.add_parser(
        "add-record-calls",
        help="add the calls of the library run on the record assets to the existing call set",
    )
    _fr.add_arguments_add(fr_add)

    rmap = commands.add_parser(
        "build-repeat-map",
        help="the repeated part of the reference, from its alignment to itself (BLAST+ on PATH): the map that stands in for a TE "
        "annotation in runs without a library (see --repeat-map)",
    )
    rmap_ref = rmap.add_mutually_exclusive_group(required=True)
    rmap_ref.add_argument("--reference-fasta", type=Path)
    rmap_ref.add_argument("--reference-assets", type=Path)
    rmap.add_argument("--output", type=Path, required=True, help="output folder; the map is <output>/repeat_map.bed")
    rmap.add_argument("--contigs", help="comma-separated contigs whose windows are aligned (default: the primary chromosomes and "
                      "the contigs >= 1 Mb, the contigs that library-free scans by default)")
    rmap.add_argument("--parts", type=int, help="blastn processes (default: the CPUs, at most 8); the map does not depend on it")

    rprof = commands.add_parser(
        "repeat-profile",
        help="repeat content of the reference windows that library-free's --repeat-cutoff genome-share reads: the window counts of an "
        "anchor reference, and the cutoff a reference gets",
    )
    rprof_ref = rprof.add_mutually_exclusive_group(required=True)
    rprof_ref.add_argument("--reference-fasta", type=Path)
    rprof_ref.add_argument("--reference-assets", type=Path)
    rprof.add_argument("--output", type=Path, help="JSON file for the profile (default: print only)")

    tiers = commands.add_parser(
        "call-tiers",
        help="confidence tier of every call of a finished run (and of its library-free calls); output only",
    )
    tiers.add_argument("--run-dir", type=Path, required=True, help="completed `evitea call` output")
    tiers.add_argument("--library-free-dir", type=Path, help="completed `evitea library-free` output: tier its added calls too")
    tiers.add_argument("--reference-assets", type=Path, required=True)
    tiers.add_argument("--te-assets", type=Path, required=True)
    tiers.add_argument("--output", type=Path, required=True)

    resident = commands.add_parser(
        "filter-resident-copies",
        help="drop core calls that touch or overlap a reference copy of their own TE family, from a finished call "
        "table, with the copies read from a RepeatMasker .out (`call` and `library-free` apply the same filter in the run, "
        "from the BED6 repeat annotation; this command filters tables of runs made without it)",
    )
    resident.add_argument("--calls", type=Path, required=True, help="calls.native.tsv or calls.combined.tsv")
    resident.add_argument("--reference-te-copies", type=Path, required=True, help="RepeatMasker .out of the reference")
    resident.add_argument("--te-assets", type=Path, required=True, help="prepared TE assets (the TE record to family map)")
    resident.add_argument("--reference-assets", type=Path,
                          help="prepared reference assets: bare copy-table contig names (2L, 1) are matched to its chr "
                          "contigs (default: to the contigs of the call table)")
    resident.add_argument("--output", type=Path, required=True)
    resident.add_argument("--allow-unlisted-copies", action="store_true",
                          help="accept a copy table whose sha256 is not listed in the tool")

    from . import run_command as _run

    run = commands.add_parser(
        "run",
        help="one command from FASTQ to calls: the prepare steps and call (and library-free), or, with --no-library, the "
        "run without a TE library (repeat map, library-free, family records, record pass)",
    )
    _run.add_arguments(run)

    validate = commands.add_parser("validate")
    targets = validate.add_subparsers(dest="validation_target", required=True)
    targets.add_parser("install")
    reference_validation = targets.add_parser("reference")
    reference_validation.add_argument("--reference-assets", type=Path, required=True)
    te_validation = targets.add_parser("te")
    te_validation.add_argument("--te-assets", type=Path, required=True)
    sample_validation = targets.add_parser("sample")
    sample_validation.add_argument("--sample-assets", type=Path, required=True)
    run_validation = targets.add_parser("run")
    run_validation.add_argument("--output", type=Path, required=True)

    check = commands.add_parser(
        "check-install",
        help="list every program and Python module EviTEA uses, with its version, or say that it is missing",
    )
    check.add_argument("--json", action="store_true")

    version = commands.add_parser("version")
    version.add_argument("--json", action="store_true")
    return parser


def _prepare_reference(args: argparse.Namespace) -> None:
    reference = validate_reference_fasta(args.reference)
    repeat = (
        load_repeat_annotation(args.repeat_annotation, reference)
        if args.repeat_annotation
        else None
    )
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise OutputStateError("prepared-reference output must be absent or empty")
    bowtie_dir = output / "bowtie2"
    bowtie_dir.mkdir(parents=True, exist_ok=True)
    prefix = bowtie_dir / "reference"
    if args.bowtie2_index is not None:
        source = args.bowtie2_index.expanduser().resolve()
        index_files = sorted(
            {*source.parent.glob(source.name + ".*.bt2"), *source.parent.glob(source.name + ".*.bt2l")}
        )
        if not index_files:
            raise InputContractError(f"no Bowtie2 index files found for prefix {source}")
        for path in index_files:
            shutil.copyfile(path, bowtie_dir / ("reference" + path.name[len(source.name):]))
        executable = f"adopted:{source}"
    else:
        executable = shutil.which(str(args.bowtie2_build))
        if executable is None:
            raise OperationalResourceError("bowtie2-build executable is unavailable")
        completed = subprocess.run(
            (executable, str(reference.path), str(prefix)),
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise OperationalResourceError(
                f"bowtie2-build failed with exit code {completed.returncode}"
            )
    index_paths = tuple(
        path.relative_to(output).as_posix()
        for path in sorted(bowtie_dir.iterdir())
        if path.is_file()
    )
    prepare_reference(
        reference,
        output,
        mapper_profile={
            "tool": "bowtie2-build" if args.bowtie2_index is None else "bowtie2-index-adopted",
            "executable": executable,
            "index_prefix": prefix.relative_to(output).as_posix(),
        },
        index_relative_paths=index_paths,
        repeat_context=repeat,
    )


def _prepare_te(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise OutputStateError("prepared-TE output must be absent or empty")
    te = validate_te_fasta(args.te_library)
    family_map = load_family_map(args.family_map)
    check_te_records((record.record_id for record in te.records), family_map.assignments)
    reference = load_prepared_reference(args.reference_assets)
    prepare_te_assets(
        te,
        family_map,
        args.family_map,
        reference,
        output,
        g1a_archive_path=args.g1a_archive,
        family_crosswalk_path=args.family_crosswalk,
        terminal_codebook_path=args.terminal_codebook,
        internal_codebook_path=args.internal_codebook,
        allow_unassigned_records=is_flybase_family_table(args.family_map),
        build_languages=args.build_languages,
        language_read_length=validate_read_length(args.read_length),
        progress=lambda message: print(f"[evitea prepare-te] {message}", file=sys.stderr, flush=True),
    )


def _prepare_sample(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise OutputStateError("prepared-sample output must be absent or empty")
    check_paired_reads(args.read1, args.read2)
    sample_id = derive_sample_id(args.read1, args.sample_id)
    threads = RuntimeConfig(threads=args.threads).resolved_threads
    try:
        build_prepared_msbwt_sample(
            args.read1,
            args.read2,
            output,
            sample_id=sample_id,
            threads=threads,
        )
    except PairSynchronizationError:
        raise
    except InputContractError as exc:
        raise RawInputValidationError(str(exc)) from exc


def _validate_run(output: Path) -> None:
    manifest_path = output / "run_manifest.json"
    try:
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AssetCompatibilityError("run manifest is missing or malformed") from exc
    if document.get("architecture_id") != ARCHITECTURE_ID:
        raise AssetCompatibilityError("run architecture identity mismatch")
    if document.get("schema_version") != RUN_MANIFEST_SCHEMA:
        raise AssetCompatibilityError("run manifest schema mismatch")
    values = document.get("values")
    if not isinstance(values, dict):
        raise AssetCompatibilityError("run manifest values are malformed")
    output_hashes = values.get("output_hashes")
    if not isinstance(output_hashes, dict):
        raise AssetCompatibilityError("run output hashes are malformed")
    required = set(CANONICAL_OUTPUTS) - {"run_manifest.json"}
    if set(output_hashes) - set(OPTIONAL_OUTPUTS) != required:
        raise AssetCompatibilityError(
            "run manifest must bind every canonical non-manifest output"
        )
    for filename, expected in output_hashes.items():
        try:
            actual = hash_file(output / filename)
        except InputContractError as exc:
            raise AssetCompatibilityError(
                f"run output is missing or unreadable: {filename}"
            ) from exc
        if actual != expected:
            raise AssetCompatibilityError(f"run output hash mismatch: {filename}")


def _validate(args: argparse.Namespace) -> None:
    target = args.validation_target
    if target == "install":
        return
    if target == "reference":
        verify_reference_asset(load_prepared_reference(args.reference_assets))
    elif target == "te":
        verify_te_assets(load_prepared_te_assets(args.te_assets))
    elif target == "sample":
        verify_prepared_sample(load_prepared_sample(args.sample_assets))
    elif target == "run":
        _validate_run(args.output)


# (name, what uses it, required, how to find the version)
_INSTALL_ITEMS = (
    ("python", "everything (CPython 3.14)", True, "python"),
    ("pymsbwt", "the msBWT of the sample (prepare-sample, call)", True, "module:MUSCython.MultiStringBWTCython"),
    ("numpy", "everything", True, "module:numpy"),
    ("bowtie2", "call (read mapping)", True, "program:bowtie2"),
    ("bowtie2-build", "prepare-reference, prepare-te (indexes)", True, "program:bowtie2-build"),
    ("edlib", "call: the paired-end fragment decoder (or --edlib-python); extra 'anchorpe'", False, "module:edlib"),
    ("mappy", "library-free, build-family-records, --refine-breakpoints; extra 'library-free'", False, "module:mappy"),
    ("parasail", "library-free; extra 'library-free'", False, "module:parasail"),
    ("makeblastdb", "build-repeat-map (BLAST+)", False, "program:makeblastdb"),
    ("blastn", "build-repeat-map (BLAST+)", False, "program:blastn"),
)


def _install_status(name: str, probe: str) -> str | None:
    """Return the version of an item, ``installed`` when it has none to report, or None when it is missing."""

    import importlib

    if probe == "python":
        version = sys.version.split()[0]
        return version if sys.version_info[:2] == (3, 14) else None
    kind, target = probe.split(":", 1)
    if kind == "module":
        try:
            importlib.import_module(target)
        except Exception:
            return None
        if name == "pymsbwt":
            return _msbwt_distribution()[1] or "installed"
        return _distribution_version(name) or "installed"
    found = shutil.which(target)
    if found is None:
        return None
    for flag in ("--version", "-version"):
        try:
            completed = subprocess.run((found, flag), check=False, capture_output=True, text=True, timeout=20)
        except (OSError, subprocess.SubprocessError):
            continue
        lines = (completed.stdout or completed.stderr).splitlines()
        if completed.returncode == 0 and lines:
            number = re.search(r"(\d+\.\d+[\w.+-]*)\s*$", lines[0])
            return number.group() if number else lines[0].strip()
    return "installed"


def _check_install(args: argparse.Namespace) -> int:
    rows = [(name, need, required, _install_status(name, probe)) for name, need, required, probe in _INSTALL_ITEMS]
    if args.json:
        print(json.dumps([{"name": n, "used_by": u, "required": r, "version": v} for n, u, r, v in rows], sort_keys=True))
    else:
        print(f"EviTEA {__version__} ({sys.platform})")
        for name, need, required, version in rows:
            state = version if version is not None else ("MISSING (required)" if required else "missing")
            print(f"  {name:<14} {state:<45} {need}")
        if sys.version_info[:2] != (3, 14):
            print(f"  EviTEA requires CPython 3.14; this is {sys.version.split()[0]}")
        if not sys.platform.startswith("linux"):
            print("  the library-free layer and the family-record commands run on Linux (or WSL) only")
    return 0 if all(version is not None for _, _, required, version in rows if required) else 8


def _version(args: argparse.Namespace) -> None:
    exact_query_package, exact_query_version = _msbwt_distribution()
    numpy_version = _distribution_version("numpy")
    native_version = _distribution_version("evitea-native")
    document = {
        "architecture_id": ARCHITECTURE_ID,
        "bowtie2_version": _executable_version("bowtie2"),
        "evitea_version": __version__,
        "exact_query_backend": {
            "name": "msbwt",
            "package": exact_query_package,
            "version": exact_query_version,
        },
        "native_backend": {
            "available": native_version is not None,
            "version": native_version,
        },
        "numpy_version": numpy_version,
        "python_version": sys.version.split()[0],
    }
    if args.json:
        print(json.dumps(document, sort_keys=True, separators=(",", ":")))
    else:
        print(f"EviTEA {__version__}")
        print(f"architecture: {ARCHITECTURE_ID}")
        print(f"Python: {document['python_version']}")
        print(f"NumPy: {numpy_version or 'unavailable'}")
        print(f"exact-query backend: msbwt ({exact_query_package} {exact_query_version or 'unavailable'})")
        print(f"native backend: {native_version or 'unavailable'}")
        print(f"Bowtie2: {document['bowtie2_version'] or 'unavailable'}")


def _call_config(args: argparse.Namespace) -> EviTEAConfig:
    base = load_config(args.config) if args.config else EviTEAConfig()
    if args.read2 is not None and args.read1 is None:
        raise InputContractError("--read2 requires --read1")
    msbwt_index = None
    if args.sample_assets is not None:
        read1 = None
        read2 = None
        sample_assets = _resolved_cli_path(args.sample_assets)
    elif args.read1 is not None:
        read1 = _resolved_cli_path(args.read1)
        read2 = _resolved_cli_path(args.read2)
        sample_assets = None
    elif args.msbwt_index is not None:
        read1 = read2 = sample_assets = None
        msbwt_index = _resolved_cli_path(args.msbwt_index)
    else:
        read1 = base.inputs.read1
        read2 = base.inputs.read2
        sample_assets = base.inputs.sample_assets
        msbwt_index = base.inputs.msbwt_index
    if read1 is not None and read2 is None:
        raise InputContractError("--read1 requires --read2")
    if sample_assets is not None and read2 is not None:
        raise InputContractError("--read2 is invalid with --sample-assets")
    if sample_assets is None and read1 is None and msbwt_index is None:
        raise InputContractError(
            "call requires --read1 with --read2, --sample-assets, --msbwt-index, "
            "or configured inputs"
        )
    return merge_config(
        base,
        inputs=InputConfig(
            sample_id=args.sample_id or base.inputs.sample_id,
            read1=read1,
            read2=read2,
            sample_assets=sample_assets,
            msbwt_index=msbwt_index,
            profile=args.profile or base.inputs.profile,
        ),
        assets=PreparedAssetConfig(
            _resolved_cli_path(args.reference_assets) or base.assets.reference_assets,
            _resolved_cli_path(args.te_assets) or base.assets.te_assets,
        ),
        runtime=RuntimeConfig(
            threads=(
                args.threads if args.threads is not None else base.runtime.threads
            ),
            work_dir=_resolved_cli_path(args.work_dir) or base.runtime.work_dir,
            exact_query_backend=(
                args.exact_query_backend or base.runtime.exact_query_backend
            ),
            resume=args.resume if args.resume is not None else base.runtime.resume,
            keep_intermediates=(
                args.keep_intermediates
                if args.keep_intermediates is not None
                else base.runtime.keep_intermediates
            ),
            bowtie2=_resolved_cli_path(args.bowtie2) or base.runtime.bowtie2,
            bowtie2_build=_resolved_cli_path(args.bowtie2_build) or base.runtime.bowtie2_build,
            edlib_python=_resolved_cli_path(args.edlib_python) or base.runtime.edlib_python,
            workers=args.workers if args.workers is not None else base.runtime.workers,
            max_rescue_interval_mass=(
                args.max_rescue_interval_mass
                if getattr(args, "max_rescue_interval_mass", None) is not None
                else base.runtime.max_rescue_interval_mass
            ),
            read_length=(
                args.read_length
                if getattr(args, "read_length", None) is not None
                else base.runtime.read_length
            ),
            fragment_q95=(
                args.fragment_q95
                if getattr(args, "fragment_q95", None) is not None
                else base.runtime.fragment_q95
            ),
            fragment_q99=(
                args.fragment_q99
                if getattr(args, "fragment_q99", None) is not None
                else base.runtime.fragment_q99
            ),
            primary_contigs=(
                parse_contig_list(args.primary_contigs)
                if getattr(args, "primary_contigs", None) is not None
                else base.runtime.primary_contigs
            ),
        ),
        output=OutputConfig(
            _resolved_cli_path(args.output) or base.output.directory,
            args.log_level or base.output.log_level,
        ),
        refinement=RefinementConfig(
            breakpoints=(
                args.refine_breakpoints
                if args.refine_breakpoints is not None
                else base.refinement.breakpoints
            ),
            complete_routing=(
                args.complete_routing
                if getattr(args, "complete_routing", None) is not None
                else base.refinement.complete_routing
            ),
            pin_one_sided=(
                args.pin_one_sided
                if getattr(args, "pin_one_sided", None) is not None
                else base.refinement.pin_one_sided
            ),
            call_annotations=(
                args.call_annotations
                if getattr(args, "call_annotations", None) is not None
                else base.refinement.call_annotations
            ),
        ),
    )


def _breakpoint_refiner(
    *,
    msbwt_index: Path | None,
    assets: PreparedAssets,
    work_dir: Path,
    progress: Callable[[str], None] | None = None,
    pin_one_sided: bool = False,
    call_annotations: bool = False,
):
    """Return the breakpoint refiner, imported only when breakpoint refinement is on."""

    from .breakpoints.runtime import BreakpointRefiner

    if msbwt_index is None:
        if assets.sample is None:
            raise InputContractError(
                "breakpoint refinement requires --msbwt-index or --sample-assets "
                "(run prepare-sample first for raw reads)"
            )
        profile = assets.sample.backend_profile
        msbwt_index = assets.sample.root / str(profile.get("index_directory", "index"))
    return BreakpointRefiner(
        msbwt_index=Path(msbwt_index),
        reference_fasta=assets.reference.reference_fasta_path,
        repeat_annotation=assets.reference.repeat_annotation_path,
        te_fasta=assets.te.source_path("te_library"),
        family_map=assets.te.source_path("family_map"),
        work_dir=work_dir,
        progress=progress,
        pin_one_sided=pin_one_sided,
        call_annotations=call_annotations,
    )


def _refine_breakpoints(args: argparse.Namespace) -> None:
    """Refine an existing run: calls.native.tsv (unchanged copy), calls.tsv, calls.bed, breakpoints.tsv."""

    _apply_depth_cap(args)
    from .breakpoints.annotate import call_annotations_tsv_bytes
    from .breakpoints.refine import (
        RefinementResult,
        breakpoints_tsv_bytes,
        refined_coordinates,
        refinement_summary,
    )
    from .pipeline import ANNOTATIONS_OUTPUT
    from .callset import CallsetState, calls_bed_bytes, calls_tsv_bytes, calls_vcf_bytes
    from .outputs import write_bytes, write_run_summary
    from .pipeline import VCF_OUTPUT, prepare_output_directory
    from .resident_filter import RESIDENT_DROPPED_OUTPUT, RESIDENT_SUMMARY_KEY

    run_dir = _resolved_cli_path(args.run_dir)
    native = (run_dir / "calls.native.tsv").read_bytes()
    try:
        summary = json.loads((run_dir / "run_summary.json").read_text(encoding="utf-8"))
        sample_id = summary["sample_id"]
    except (OSError, json.JSONDecodeError, KeyError) as exc:
        raise InputContractError(f"run directory has no readable run_summary.json: {run_dir}") from exc
    config = EviTEAConfig(
        inputs=InputConfig(sample_assets=_resolved_cli_path(args.sample_assets)),
        assets=PreparedAssetConfig(
            _resolved_cli_path(args.reference_assets), _resolved_cli_path(args.te_assets)
        ),
    )
    assets = prepare_assets(config)
    resident_filter = None if args.keep_resident_copies else _own_family_copy_filter(assets.reference, assets.te)
    output = prepare_output_directory(_resolved_cli_path(args.output), resume=False)
    started = time.time()

    def progress(message: str) -> None:
        print(f"[evitea {time.time() - started:9.1f}s] {message}", file=sys.stderr, flush=True)

    refiner = _breakpoint_refiner(
        msbwt_index=_resolved_cli_path(args.msbwt_index),
        assets=assets,
        work_dir=_resolved_cli_path(args.work_dir) or (output / ".evitea-work"),
        progress=progress,
        pin_one_sided=args.pin_one_sided,
        call_annotations=args.call_annotations,
    )
    state = CallsetState(sample_id=sample_id, native=native)
    result = RefinementResult.of(refiner(state))
    rows, pins, annotations = result.rows, result.pins, result.annotations
    refined = refined_coordinates(rows)
    shown, resident = state, None
    if resident_filter is not None:  # the user-facing tables only; calls.native.tsv stays the full table
        from .resident_filter import dropped_tsv_bytes

        resident = resident_filter(native)
        shown = CallsetState(sample_id=sample_id, native=resident.filtered)
        dropped = frozenset(row.event_id for row in resident.dropped)
        refined = {call_id: span for call_id, span in refined.items() if call_id not in dropped}
        rows = [row for row in rows if row.event_id not in dropped]
        pins = {call_id: pin for call_id, pin in pins.items() if call_id not in dropped} if pins is not None else None
        annotations = [row for row in annotations if row.event_id not in dropped] if annotations is not None else None
    write_bytes(output / "calls.native.tsv", native)
    write_bytes(output / "calls.tsv", calls_tsv_bytes(shown, refined))
    write_bytes(
        output / "calls.bed", calls_bed_bytes(shown, dict(assets.reference.contigs), refined)
    )
    write_bytes(output / BREAKPOINTS_OUTPUT, breakpoints_tsv_bytes(rows, pins=pins))
    write_bytes(
        output / VCF_OUTPUT,
        calls_vcf_bytes(shown, assets.reference.contigs, assets.reference.reference_fasta_path, refined),
    )
    written = ["calls.native.tsv", "calls.tsv", "calls.bed", VCF_OUTPUT, BREAKPOINTS_OUTPUT]
    if resident is not None:
        write_bytes(output / RESIDENT_DROPPED_OUTPUT, dropped_tsv_bytes(resident.dropped))
        written.append(RESIDENT_DROPPED_OUTPUT)
    if annotations is not None:
        write_bytes(output / ANNOTATIONS_OUTPUT, call_annotations_tsv_bytes(annotations))
        written.append(ANNOTATIONS_OUTPUT)
    record = {
        "sample_id": sample_id,
        "source_run": str(run_dir),
        "native_callset_sha256": state.native_sha256,
        "aligner": f"mappy {refiner.aligner_version}",
        "depth_cap": args.depth_cap,
        "seconds": round(time.time() - started, 1),
        **refinement_summary(rows),
        **({"pinned_calls": len(pins)} if pins is not None else {}),
        **({"annotated_calls": len(annotations)} if annotations is not None else {}),
        **({RESIDENT_SUMMARY_KEY: resident_filter.summary(resident)} if resident is not None else {}),
        "output_hashes": {name: hash_file(output / name) for name in written},
    }
    write_run_summary(output / "refinement_summary.json", record)


def _build_family_records(args: argparse.Namespace) -> None:
    from .family_records.run import build_family_records

    if args.reference_fasta is not None:
        reference_fasta = _resolved_cli_path(args.reference_fasta)
    else:
        reference_fasta = load_prepared_reference(_resolved_cli_path(args.reference_assets)).reference_fasta_path
    started = time.time()
    build_family_records(
        family_calls=[_resolved_cli_path(p) for p in args.family_calls],
        segments=_resolved_cli_path(args.segments),
        second_segments=_resolved_cli_path(args.second_segments),
        reference_fasta=reference_fasta,
        output=_resolved_cli_path(args.output),
        flag_rule=args.flag_rule,
        record_guard=args.record_guard == "on",
        repeat_annotation=_resolved_cli_path(args.repeat_annotation),
        repeat_map=_resolved_cli_path(args.repeat_map),
        te_share=args.te_share,
        copy_report=args.copy_report == "on",
        progress=lambda m: print(f"[evitea {time.time() - started:9.1f}s] {m}", file=sys.stderr, flush=True),
    )


def _build_repeat_map(args: argparse.Namespace) -> None:
    from .repeat_map import build_repeat_map

    if args.reference_fasta is not None:
        reference_fasta = _resolved_cli_path(args.reference_fasta)
    else:
        reference_fasta = load_prepared_reference(_resolved_cli_path(args.reference_assets)).reference_fasta_path
    started = time.time()
    summary = build_repeat_map(
        reference_fasta=reference_fasta,
        output=_resolved_cli_path(args.output),
        contigs=tuple(c for c in (args.contigs or "").split(",") if c) or None,
        parts=args.parts,
        progress=lambda m: print(f"[evitea {time.time() - started:9.1f}s] {m}", file=sys.stderr, flush=True),
    )
    print(f"intervals {summary['intervals']} bases {summary['bases']}")


def _repeat_profile(args: argparse.Namespace) -> None:
    import hashlib

    import numpy as np

    from .breakpoints.runtime import read_fasta
    from .library_free import _p4_lib as P
    from .library_free import cross as X

    if args.reference_fasta is not None:
        reference_fasta = _resolved_cli_path(args.reference_fasta)
    else:
        reference_fasta = load_prepared_reference(_resolved_cli_path(args.reference_assets)).reference_fasta_path
    digest = hashlib.sha256()
    with Path(reference_fasta).open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    genome = read_fasta(Path(reference_fasta))
    rm = P.RefMulti(P.multi_codes(genome, P.REP_K), np.zeros(0, dtype=np.uint64))
    reps = X.repeat_profile(genome, rm)
    out = {"reference_fasta": str(reference_fasta), "reference_sha256": digest.hexdigest(),
           "contigs": [name for name, seq in genome.items() if len(seq) >= X.MIN_PROFILE_CONTIG],
           "window_bp": 2 * P.REP_FLANK, "rep_max": P.REP_MAX, "windows": int(len(reps)),
           "at_cutoff": int(np.count_nonzero(reps >= P.REP_MAX)),
           "quantiles": {str(q): round(float(np.quantile(reps, q)), 4) for q in (0.5, 0.75, 0.9, 0.95, 0.99)} if len(reps) else {}}
    if X.REPEAT_ANCHOR is not None:
        out["cutoff"] = X.repeat_cutoff(reps)[1]
    text = json.dumps(out, indent=2)
    if args.output is not None:
        Path(_resolved_cli_path(args.output)).write_text(text + "\n", encoding="utf-8")
    print(text)


def _add_record_calls(args: argparse.Namespace) -> None:
    from .family_records.run import add_record_calls

    root = _resolved_cli_path(args.records)
    report = add_record_calls(
        calls=_resolved_cli_path(args.calls),
        base_calls=_resolved_cli_path(args.base_calls),
        annotation=root / "annotation.bed",
        family_map=root / "family_map.tsv",
        output=_resolved_cli_path(args.output),
    )
    print(json.dumps(report, sort_keys=True))


def _library_free(args: argparse.Namespace) -> None:
    """Run the library-free detector on an existing run: places, family calls, added calls, combined calls."""

    _apply_depth_cap(args)
    from .library_free.runtime import run_library_free

    reference = load_prepared_reference(_resolved_cli_path(args.reference_assets))
    msbwt_index = _resolved_cli_path(args.msbwt_index)
    if msbwt_index is None and args.sample_assets is not None:
        sample = load_prepared_sample(_resolved_cli_path(args.sample_assets))
        msbwt_index = sample.root / str(sample.backend_profile.get("index_directory", "index"))
    if msbwt_index is None and args.detection_dir is None:
        raise InputContractError("library-free detection requires --msbwt-index or --sample-assets (or --detection-dir)")
    resident_filter = None
    if not args.keep_resident_copies:
        if args.te_assets is not None:
            resident_filter = _own_family_copy_filter(reference, load_prepared_te_assets(_resolved_cli_path(args.te_assets)))
        else:
            native_path = _resolved_cli_path(args.run_dir) / "calls.native.tsv"
            if native_path.is_file() and len(native_path.read_bytes().splitlines()) > 1:
                raise InputContractError(
                    "the run has library calls: pass --te-assets (the family map of the own-family reference copy filter) "
                    "or --keep-resident-copies"
                )
    started = time.time()
    from .library_free.options import from_args as _lf_options
    from .library_free.single_copy import Choice as _single_copy
    from .library_free.tsd_skip import Skip as _tsd_skip

    try:
        options = _lf_options(args)
        single_copy = _single_copy.parse(args.single_copy_contigs)
        tsd_skip = _tsd_skip(on=args.rule_c_tsd_skip == "on", share=args.rule_c_tsd_share, minimum=args.rule_c_tsd_min)
    except ValueError as exc:
        raise InputContractError(str(exc)) from exc

    def progress(message: str) -> None:
        print(f"[evitea {time.time() - started:9.1f}s] {message}", file=sys.stderr, flush=True)

    run_library_free(
        run_dir=_resolved_cli_path(args.run_dir),
        reference=reference,
        msbwt_index=msbwt_index,
        output=_resolved_cli_path(args.output),
        work_dir=_resolved_cli_path(args.work_dir),
        contigs=tuple(c for c in (args.contigs or "").split(",") if c) or None,
        rule_c=args.rule_c == "on",
        detection_dir=_resolved_cli_path(args.detection_dir),
        command_line=args.command_line,
        source_identity=_source_identity(),
        progress=progress,
        options=options,
        repeat_annotation=_resolved_cli_path(Path(args.repeat_annotation)) if args.repeat_annotation else None,
        repeat_map=_resolved_cli_path(Path(args.repeat_map)) if args.repeat_map else None,
        mito_contig=args.mito_contig,
        resident_filter=resident_filter,
        adapter_test=args.adapter_test == "on",
        single_copy=single_copy,
        het_guard=args.het_guard,
        tsd_skip=tsd_skip,
        workers=args.workers,
        counter=args.counter,
        reference_cache=_resolved_cli_path(args.reference_cache),
    )


def _genotype(args: argparse.Namespace) -> None:
    """Write genotypes.tsv for a finished library-free output: Step A from the recorded walk fractions, Step B from the msBWT."""

    from .genotype import run_genotype

    msbwt_index = _resolved_cli_path(args.msbwt_index)
    if msbwt_index is None and args.sample_assets is not None:
        sample = load_prepared_sample(_resolved_cli_path(args.sample_assets))
        msbwt_index = sample.root / str(sample.backend_profile.get("index_directory", "index"))
    reference_fasta = None
    if args.reference_assets is not None:
        reference_fasta = Path(load_prepared_reference(_resolved_cli_path(args.reference_assets)).reference_fasta_path)
    if (msbwt_index is None) != (reference_fasta is None):
        raise InputContractError("the library calls need both --msbwt-index (or --sample-assets) and --reference-assets")
    started = time.time()
    library_free_dir = _resolved_cli_path(args.library_free_dir)

    def progress(message: str) -> None:
        print(f"[evitea {time.time() - started:9.1f}s] {message}", file=sys.stderr, flush=True)

    run_genotype(
        library_free_dir=library_free_dir,
        output=_resolved_cli_path(args.output) or library_free_dir,
        reference_fasta=reference_fasta,
        msbwt_index=msbwt_index,
        progress=progress,
        command_line=args.command_line,
        source_identity=_source_identity(),
        counter=args.counter,
    )


def _tier_resident_paths(args: argparse.Namespace) -> tuple[Path, Path]:
    """Return (repeat annotation, TE family map) paths of the prepared assets."""

    reference = load_prepared_reference(_resolved_cli_path(args.reference_assets))
    te = load_prepared_te_assets(_resolved_cli_path(args.te_assets))
    return reference.repeat_annotation_path, te.source_path("family_map")


def _call_tiers(args: argparse.Namespace) -> None:
    """Write call_tiers.tsv, one confidence tier per call of the combined call table; the run files are only read."""

    from collections import Counter

    from .outputs import write_bytes, write_run_summary
    from .pipeline import prepare_output_directory
    from .tiers import MIN_SUPPORT, NEAR_BP, STRICT_TIERS, call_tiers_tsv_bytes, load_resident_index, tier_rows

    run_dir = _resolved_cli_path(args.run_dir)
    native_path = run_dir / "calls.native.tsv"
    added_path = _resolved_cli_path(args.library_free_dir) / "calls.library_free.tsv" if args.library_free_dir else None
    if not native_path.is_file():
        raise InputContractError(f"run directory has no calls.native.tsv: {run_dir}")
    if added_path is not None and not added_path.is_file():
        raise InputContractError(f"library-free directory has no calls.library_free.tsv: {added_path.parent}")
    annotation, family_map = _tier_resident_paths(args)
    output = prepare_output_directory(_resolved_cli_path(args.output), resume=False)
    rows = tier_rows(native_path.read_bytes(), added_path.read_bytes() if added_path else None,
                     load_resident_index(annotation, family_map))
    write_bytes(output / "call_tiers.tsv", call_tiers_tsv_bytes(rows))
    inputs = {"calls.native.tsv": native_path, **({"calls.library_free.tsv": added_path} if added_path else {}),
              "repeat_annotation": annotation, "family_map": family_map}
    write_run_summary(output / "call_tiers_summary.json", {
        "source_run": str(run_dir),
        "library_free_dir": str(added_path.parent) if added_path else None,
        "calls": len(rows),
        "tiers": {str(tier): n for tier, n in sorted(Counter(r.tier for r in rows).items())},
        "parameters": {"min_support": MIN_SUPPORT, "near_bp": NEAR_BP, "strict_tiers": list(STRICT_TIERS)},
        "input_hashes": {name: hash_file(path) for name, path in inputs.items()},
        "output_hashes": {"call_tiers.tsv": hash_file(output / "call_tiers.tsv")},
    })


def _own_family_copy_filter(reference, te):
    """Return the in-run own-family copy filter of prepared reference and TE assets (BED6 annotation, family map)."""

    from .resident_filter import OwnFamilyCopyFilter

    return OwnFamilyCopyFilter(reference.repeat_annotation_path, te.source_path("family_map"),
                               tuple(name for name, _ in reference.contigs))


def _resident_family_map_path(args: argparse.Namespace) -> Path:
    """Return the TE record to family map of the prepared TE assets."""

    return load_prepared_te_assets(_resolved_cli_path(args.te_assets)).source_path("family_map")


def _filter_resident_copies(args: argparse.Namespace) -> None:
    """Write the call table without the core calls inside a reference copy of their own family (an ablation)."""

    from .outputs import write_bytes, write_run_summary
    from .pipeline import prepare_output_directory
    from .resident_filter import (
        CopyIndex,
        FamilyIds,
        call_contigs,
        dropped_tsv_bytes,
        filter_calls,
        load_copy_table,
        load_family_ids,
        verify_copy_table,
    )

    calls_path = _resolved_cli_path(args.calls)
    copies_path = _resolved_cli_path(args.reference_te_copies)
    if not calls_path.is_file():
        raise InputContractError(f"call table not found: {calls_path}")
    if not copies_path.is_file():
        raise InputContractError(f"reference TE copy table not found: {copies_path}")
    copies_sha, listed = verify_copy_table(copies_path, allow_unlisted=args.allow_unlisted_copies)
    family_map = _resident_family_map_path(args)
    started = time.time()
    calls = calls_path.read_bytes()
    contigs = (
        tuple(name for name, _ in load_prepared_reference(_resolved_cli_path(args.reference_assets)).contigs)
        if getattr(args, "reference_assets", None) is not None
        else call_contigs(calls)
    )
    index = CopyIndex(load_copy_table(copies_path), contigs)
    result = filter_calls(calls, index, load_family_ids(family_map))
    output = prepare_output_directory(_resolved_cli_path(args.output), resume=False)
    write_bytes(output / calls_path.name, result.filtered)
    write_bytes(output / "resident_copy_dropped.tsv", dropped_tsv_bytes(result.dropped))
    write_run_summary(output / "resident_copy_filter_summary.json", {
        "calls": str(calls_path),
        "kept": result.kept,
        "dropped": len(result.dropped),
        "unmapped_family": result.unmapped_family,
        "library_free_passed": result.library_free_passed,
        "copies": len(index),
        "copy_table_listed": listed,
        "parameters": {"gap_bp": 0, "scope": "core calls (event_id not starting LF_)"},
        "input_hashes": {"calls": hash_file(calls_path), "reference_te_copies": copies_sha,
                         "family_map": hash_file(family_map)},
        "output_hashes": {calls_path.name: hash_file(output / calls_path.name),
                          "resident_copy_dropped.tsv": hash_file(output / "resident_copy_dropped.tsv")},
        "wall_seconds": round(time.time() - started, 3),
    })


def _call(
    args: argparse.Namespace,
    call_engine: CallEngine | None,
    sample_builder: RawSampleBuilder | None = None,
) -> RunResult:
    config = _call_config(args)
    _apply_depth_cap(args)
    if config.inputs.read1 is not None:
        check_paired_reads(config.inputs.read1, config.inputs.read2)
    elif config.inputs.msbwt_index is not None and (Path(config.inputs.msbwt_index) / "about.npy").is_file():
        check_msbwt_index(Path(config.inputs.msbwt_index))
    if call_engine is None:
        call_engine = _default_engine_with_progress(config)
    assets = prepare_assets(config)
    resident_filter = None if args.keep_resident_copies else _own_family_copy_filter(assets.reference, assets.te)
    if config.runtime.primary_contigs is not None:
        check_contig_names(config.runtime.primary_contigs, dict(assets.reference.contigs), "--primary-contigs")
    input_assets = {}
    if config.inputs.read1 is not None:
        input_assets = {
            "read1": file_identity(config.inputs.read1),
            "read2": file_identity(config.inputs.read2),
        }
    elif config.inputs.msbwt_index is not None:
        index = Path(config.inputs.msbwt_index)
        files = sorted(path for path in index.iterdir() if path.is_file() and path.suffix == ".npy")
        if not any(path.name == "about.npy" for path in files):
            raise InputContractError(f"msBWT index has no about.npy: {index}")
        input_assets = {
            "msbwt_index": {path.name: file_identity(path) for path in files},
        }
    context = RunManifestContext(
        normalized_config=config.normalized(),
        command_line=args.command_line,
        source_identity=_source_identity(),
        input_assets=input_assets,
        prepared_identities={
            "reference": assets.reference.identity_sha256,
            "te": assets.te.identity_sha256,
            **(
                {"sample": assets.sample.identity_sha256}
                if assets.sample is not None
                else {}
            ),
        },
        backend_versions={
            "exact_query": (
                f"{config.runtime.exact_query_backend} "
                f"{' '.join(str(part or 'unavailable') for part in _msbwt_distribution())}"
            ),
            "bowtie2": (
                _executable_version(config.runtime.bowtie2 or "bowtie2")
                or "unavailable"
            ),
        },
        resolved_threads=config.runtime.resolved_threads,
    )
    refiner = None
    if config.refinement.breakpoints:
        if config.output.directory is None:
            raise InputContractError("output directory is required")
        output = Path(config.output.directory)
        refiner = _breakpoint_refiner(
            msbwt_index=config.inputs.msbwt_index,
            assets=assets,
            work_dir=(config.runtime.work_dir or (output / ".evitea-work")) / "breakpoints",
            pin_one_sided=config.refinement.pin_one_sided,
            call_annotations=config.refinement.call_annotations,
        )
    return run_sample(
        config,
        assets,
        call_engine,
        manifest_context=context,
        sample_builder=sample_builder,
        breakpoint_refiner=refiner,
        resident_filter=resident_filter,
    )


def _default_engine_with_progress(config: EviTEAConfig) -> CallEngine:
    """Return the default call engine, reporting progress on stderr."""

    from .pipeline import default_call_engine

    started = time.time()
    quiet = config.output.log_level in ("ERROR", "WARNING")

    def progress(message: str) -> None:
        if not quiet:
            print(f"[evitea {time.time() - started:9.1f}s] {message}", file=sys.stderr, flush=True)

    def engine(run_config, assets, checkpoints):
        return default_call_engine(run_config, assets, checkpoints, progress=progress)

    return engine


def exit_code_for_error(error: BaseException) -> int:
    if isinstance(error, OutputFinalizationError):
        return 10
    if isinstance(error, ResumeStateError):
        return 6
    if isinstance(error, OutputStateError):
        return 7
    if isinstance(error, PairSynchronizationError):
        return 3
    if isinstance(error, RawInputValidationError):
        return 3
    if isinstance(error, AssetCompatibilityError):
        return 4
    if isinstance(error, BackendCompatibilityError):
        return 5
    if isinstance(error, OperationalResourceError):
        return 8
    if isinstance(error, ArchitectureInvariantError):
        return 9
    if isinstance(error, (InputContractError, argparse.ArgumentError)):
        return 2
    if isinstance(error, OSError):
        return 7
    return 9


def main(
    argv: Sequence[str] | None = None,
    *,
    call_engine: CallEngine | None = None,
    sample_builder: RawSampleBuilder | None = None,
) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        args.command_line = (
            ("evitea", *argv) if argv is not None else tuple(sys.argv)
        )
        if args.command == "prepare-reference":
            _prepare_reference(args)
        elif args.command == "prepare-te":
            _prepare_te(args)
        elif args.command == "prepare-sample":
            _prepare_sample(args)
        elif args.command == "call":
            _call(args, call_engine, sample_builder)
        elif args.command == "refine-breakpoints":
            _refine_breakpoints(args)
        elif args.command == "library-free":
            _library_free(args)
        elif args.command == "genotype":
            _genotype(args)
        elif args.command == "build-family-records":
            _build_family_records(args)
        elif args.command == "build-repeat-map":
            _build_repeat_map(args)
        elif args.command == "repeat-profile":
            _repeat_profile(args)
        elif args.command == "add-record-calls":
            _add_record_calls(args)
        elif args.command == "call-tiers":
            _call_tiers(args)
        elif args.command == "filter-resident-copies":
            _filter_resident_copies(args)
        elif args.command == "run":
            from .run_command import run as _run_steps

            started = time.time()
            return _run_steps(
                args, main,
                lambda m: print(f"[evitea {time.time() - started:9.1f}s] {m}", file=sys.stderr, flush=True),
            )
        elif args.command == "validate":
            _validate(args)
        elif args.command == "version":
            _version(args)
        elif args.command == "check-install":
            return _check_install(args)
        return 0
    except EviTEAError as error:
        print(f"evitea: {error}", file=sys.stderr)
        return exit_code_for_error(error)
    except OSError as error:
        print(f"evitea: {error}", file=sys.stderr)
        return exit_code_for_error(error)
