"""Thin public command parser; frozen scientific thresholds are absent."""

from __future__ import annotations

import argparse
import json
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
    RuntimeConfig,
    load_config,
    merge_config,
)
from .constants import ARCHITECTURE_ID
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
    CANONICAL_OUTPUTS,
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


def _source_identity() -> dict[str, object]:
    identity: dict[str, object] = {"package_version": __version__}
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
    te.add_argument("--family-map", type=Path, required=True)
    te.add_argument("--g1a-archive", type=Path)
    te.add_argument("--family-crosswalk", type=Path)
    te.add_argument(
        "--terminal-codebook",
        type=Path,
        help="final safe terminal 21-mer codebook (v3 language; with --internal-codebook)",
    )
    te.add_argument(
        "--internal-codebook",
        type=Path,
        help="internal-only safe 21-mer codebook (v3 language; with --terminal-codebook)",
    )
    te.add_argument(
        "--build-languages",
        action="store_true",
        help="build both v3 languages from the TE library, family map, reference and "
        "RepeatMasker annotation (evitea.languages)",
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
        help="'production' (default) or a reproduction-only profile (AKA017, ISO1)",
    )
    call.add_argument("--bowtie2", type=Path)
    call.add_argument("--bowtie2-build", type=Path)
    call.add_argument(
        "--edlib-python",
        type=Path,
        help="interpreter with edlib for the fragment decoder, used only if this one lacks it",
    )
    call.add_argument("--workers", type=int, help="locus-batch worker processes")
    _add_engineering_options(call, use_defaults=False)

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
    reference = load_prepared_reference(args.reference_assets)
    te = validate_te_fasta(args.te_library)
    family_map = load_family_map(args.family_map)
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
        progress=lambda message: print(f"[evitea prepare-te] {message}", file=sys.stderr, flush=True),
    )


def _prepare_sample(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise OutputStateError("prepared-sample output must be absent or empty")
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
    if set(output_hashes) != required:
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


def _version(args: argparse.Namespace) -> None:
    exact_query_version = _distribution_version("msbwt-modern3", "msbwt")
    numpy_version = _distribution_version("numpy")
    native_version = _distribution_version("evitea-native")
    document = {
        "architecture_id": ARCHITECTURE_ID,
        "bowtie2_version": _executable_version("bowtie2"),
        "evitea_version": __version__,
        "exact_query_backend": {
            "name": "msbwt",
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
        print(f"exact-query backend: msbwt {exact_query_version or 'unavailable'}")
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
        ),
        output=OutputConfig(
            _resolved_cli_path(args.output) or base.output.directory,
            args.log_level or base.output.log_level,
        ),
    )


def _call(
    args: argparse.Namespace,
    call_engine: CallEngine | None,
    sample_builder: RawSampleBuilder | None = None,
) -> RunResult:
    config = _call_config(args)
    if call_engine is None:
        call_engine = _default_engine_with_progress(config)
    assets = prepare_assets(config)
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
                f"{_distribution_version('msbwt-modern3', 'msbwt') or 'unavailable'}"
            ),
            "bowtie2": (
                _executable_version(config.runtime.bowtie2 or "bowtie2")
                or "unavailable"
            ),
        },
        resolved_threads=config.runtime.resolved_threads,
    )
    return run_sample(
        config,
        assets,
        call_engine,
        manifest_context=context,
        sample_builder=sample_builder,
    )


def _default_engine_with_progress(config: EviTEAConfig) -> CallEngine:
    """The installed default engine, reporting progress on stderr."""

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
        elif args.command == "validate":
            _validate(args)
        elif args.command == "version":
            _version(args)
        return 0
    except EviTEAError as error:
        print(f"evitea: {error}", file=sys.stderr)
        return exit_code_for_error(error)
    except OSError as error:
        print(f"evitea: {error}", file=sys.stderr)
        return exit_code_for_error(error)
