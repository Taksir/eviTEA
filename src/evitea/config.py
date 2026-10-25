"""Typed engineering configuration with immutable scientific constants."""

from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .constants import (
    ARCHITECTURE_ID,
    EVENT_EQUIV_DISTANCE,
    MIN_FOLDS,
    MIN_PAIRS,
    N_FOLDS,
    Q99,
)
from .contracts import stable_to_primitive
from .errors import InputContractError

_LOG_LEVELS = frozenset({"ERROR", "WARNING", "INFO", "DEBUG"})
_SCIENCE_KEYS = frozenset(
    {
        "q99",
        "min_pairs",
        "min_folds",
        "n_folds",
        "equivalence_distance",
        "edit_threshold",
        "family_alias",
        "truth",
    }
)


@dataclass(frozen=True, slots=True)
class ScienceProfile:
    architecture_id: str = ARCHITECTURE_ID
    q99: int = Q99
    min_pairs: int = MIN_PAIRS
    min_folds: int = MIN_FOLDS
    n_folds: int = N_FOLDS
    equivalence_distance: int = EVENT_EQUIV_DISTANCE


@dataclass(frozen=True, slots=True)
class InputConfig:
    sample_id: str | None = None
    read1: Path | None = None
    read2: Path | None = None
    sample_assets: Path | None = None
    msbwt_index: Path | None = None
    profile: str = "production"

    def __post_init__(self) -> None:
        raw_mode = self.read1 is not None or self.read2 is not None
        prepared_mode = self.sample_assets is not None
        index_mode = self.msbwt_index is not None
        if raw_mode and (self.read1 is None or self.read2 is None):
            raise InputContractError("raw sample mode requires both read1 and read2")
        if raw_mode + prepared_mode + index_mode > 1:
            raise InputContractError(
                "choose paired reads, prepared sample assets, or an msBWT index, not several"
            )
        if index_mode and not self.sample_id:
            raise InputContractError("an msBWT index input requires an explicit sample_id")
        if not isinstance(self.profile, str) or not self.profile:
            raise InputContractError("profile must be non-empty text")


@dataclass(frozen=True, slots=True)
class PreparedAssetConfig:
    reference_assets: Path | None = None
    te_assets: Path | None = None


@dataclass(frozen=True, slots=True)
class RuntimeConfig:
    threads: int = 0
    work_dir: Path | None = None
    exact_query_backend: str = "msbwt"
    resume: bool = False
    keep_intermediates: bool = False
    bowtie2: Path | None = None
    bowtie2_build: Path | None = None
    edlib_python: Path | None = None
    workers: int = 0
    # recall-repair resource limit in msBWT interval rows (frozen default 1,500,000)
    max_rescue_interval_mass: int = 1_500_000
    # read length of the sample; None = taken from the sample (see ``sample_geometry``)
    read_length: int | None = None
    # fragment-size 95th / 99th percentiles; None = estimated from the sample (see ``sample_geometry``)
    fragment_q95: float | None = None
    fragment_q99: int | None = None
    # contigs of the full-memory proposals; None = taken from the reference (see ``evitea.contigs``)
    primary_contigs: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        for name in ("threads", "workers"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise InputContractError(f"{name} must be an integer")
            if value < 0:
                raise InputContractError(f"{name} must be non-negative")
        value = self.max_rescue_interval_mass
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise InputContractError("max_rescue_interval_mass must be a positive integer")
        if self.read_length is not None:
            from .sample_geometry import validate_read_length

            validate_read_length(self.read_length)
        from .sample_geometry import validate_fragment_size

        if self.fragment_q95 is not None:
            validate_fragment_size(self.fragment_q95, "fragment_q95")
        if self.fragment_q99 is not None:
            validate_fragment_size(self.fragment_q99, "fragment_q99")
            if not isinstance(self.fragment_q99, int):
                raise InputContractError("fragment_q99 must be an integer")
        if self.primary_contigs is not None and (
            not isinstance(self.primary_contigs, tuple)
            or not self.primary_contigs
            or not all(isinstance(name, str) and name for name in self.primary_contigs)
        ):
            raise InputContractError("primary_contigs must be a non-empty list of contig names")
        if self.exact_query_backend != "msbwt":
            raise InputContractError("the frozen profile supports only msbwt")
        if not isinstance(self.resume, bool):
            raise InputContractError("resume must be a boolean")
        if not isinstance(self.keep_intermediates, bool):
            raise InputContractError("keep_intermediates must be a boolean")

    @property
    def resolved_threads(self) -> int:
        return self.threads or max(1, os.cpu_count() or 1)

    @property
    def resolved_workers(self) -> int:
        """Locus-batch worker processes (default: one fewer than the CPUs)."""

        return self.workers or max(1, (os.cpu_count() or 2) - 1)


@dataclass(frozen=True, slots=True)
class OutputConfig:
    directory: Path | None = None
    log_level: str = "INFO"

    def __post_init__(self) -> None:
        if self.log_level not in _LOG_LEVELS:
            raise InputContractError(
                "log_level must be ERROR, WARNING, INFO, or DEBUG"
            )


@dataclass(frozen=True, slots=True)
class RefinementConfig:
    """Optional post-call layers; all off by default (off = the frozen v1 outputs)."""

    breakpoints: bool = False
    # first-stage layer: complete routing with a partner-read gate
    complete_routing: bool = False
    # breakpoint-layer extras (need ``breakpoints``): one-sided junction positions for calls left
    # wide, and output-only zygosity / TSD annotations; neither adds, removes or moves a call
    pin_one_sided: bool = False
    call_annotations: bool = False

    def __post_init__(self) -> None:
        for name in ("breakpoints", "complete_routing", "pin_one_sided", "call_annotations"):
            if not isinstance(getattr(self, name), bool):
                raise InputContractError(f"refinement.{name} must be a boolean")
        for name in ("pin_one_sided", "call_annotations"):
            if getattr(self, name) and not self.breakpoints:
                raise InputContractError(f"refinement.{name} needs refinement.breakpoints")


@dataclass(frozen=True, slots=True)
class EviTEAConfig:
    inputs: InputConfig = InputConfig()
    assets: PreparedAssetConfig = PreparedAssetConfig()
    runtime: RuntimeConfig = RuntimeConfig()
    output: OutputConfig = OutputConfig()
    science: ScienceProfile = ScienceProfile()
    refinement: RefinementConfig = RefinementConfig()

    def normalized(self) -> dict[str, Any]:
        return stable_to_primitive(self)


_SECTION_KEYS = {
    "inputs": frozenset(
        {
            "sample_id",
            "read1",
            "read2",
            "sample_assets",
            "msbwt_index",
            "profile",
            "reference_assets",
            "te_assets",
        }
    ),
    "runtime": frozenset(
        {
            "threads",
            "work_dir",
            "exact_query_backend",
            "resume",
            "keep_intermediates",
            "bowtie2",
            "bowtie2_build",
            "edlib_python",
            "workers",
            "max_rescue_interval_mass",
            "read_length",
            "fragment_q95",
            "fragment_q99",
            "primary_contigs",
        }
    ),
    "output": frozenset({"directory", "log_level"}),
    "refinement": frozenset({"breakpoints", "complete_routing", "pin_one_sided", "call_annotations"}),
}


def _contig_list(value: Any) -> tuple[str, ...] | None:
    """``primary_contigs`` from TOML: a list of names, or a comma list / file path as on the command line."""

    if value is None:
        return None
    if isinstance(value, list):
        return tuple(value)
    if isinstance(value, str):
        from .contigs import parse_contig_list

        return parse_contig_list(value)
    raise InputContractError("primary_contigs must be a list of names or text")


def _optional_path(value: Any, field: str, base_dir: Path) -> Path | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value:
        raise InputContractError(f"{field} must be non-empty path text")
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base_dir / path).resolve()


def _validate_document(document: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(document, dict):
        raise InputContractError("configuration root must be a TOML table")
    unknown_sections = set(document) - set(_SECTION_KEYS)
    if unknown_sections:
        raise InputContractError(
            f"unknown configuration section(s): {sorted(unknown_sections)}"
        )
    for section, values in document.items():
        if not isinstance(values, dict):
            raise InputContractError(f"configuration section {section} must be a table")
        forbidden = set(values) & _SCIENCE_KEYS
        if forbidden:
            raise InputContractError(
                f"frozen scientific setting is not configurable: {sorted(forbidden)}"
            )
        unknown = set(values) - set(_SECTION_KEYS[section])
        if unknown:
            raise InputContractError(
                f"unknown {section} configuration key(s): {sorted(unknown)}"
            )
    return document


def load_config(path: str | Path) -> EviTEAConfig:
    """Load TOML; relative paths resolve against the TOML file's directory."""

    config_path = Path(path).expanduser().resolve()
    try:
        with config_path.open("rb") as handle:
            document = _validate_document(tomllib.load(handle))
    except OSError as exc:
        raise InputContractError(f"cannot read configuration: {config_path}") from exc
    inputs = document.get("inputs", {})
    runtime = document.get("runtime", {})
    output = document.get("output", {})
    refinement = document.get("refinement", {})
    base = config_path.parent
    return EviTEAConfig(
        inputs=InputConfig(
            sample_id=inputs.get("sample_id"),
            read1=_optional_path(inputs.get("read1"), "read1", base),
            read2=_optional_path(inputs.get("read2"), "read2", base),
            sample_assets=_optional_path(
                inputs.get("sample_assets"), "sample_assets", base
            ),
            msbwt_index=_optional_path(inputs.get("msbwt_index"), "msbwt_index", base),
            profile=inputs.get("profile", "production"),
        ),
        assets=PreparedAssetConfig(
            reference_assets=_optional_path(
                inputs.get("reference_assets"), "reference_assets", base
            ),
            te_assets=_optional_path(inputs.get("te_assets"), "te_assets", base),
        ),
        runtime=RuntimeConfig(
            threads=runtime.get("threads", 0),
            work_dir=_optional_path(runtime.get("work_dir"), "work_dir", base),
            exact_query_backend=runtime.get("exact_query_backend", "msbwt"),
            resume=runtime.get("resume", False),
            keep_intermediates=runtime.get("keep_intermediates", False),
            bowtie2=_optional_path(runtime.get("bowtie2"), "bowtie2", base),
            bowtie2_build=_optional_path(runtime.get("bowtie2_build"), "bowtie2_build", base),
            edlib_python=_optional_path(runtime.get("edlib_python"), "edlib_python", base),
            workers=runtime.get("workers", 0),
            max_rescue_interval_mass=runtime.get("max_rescue_interval_mass", 1_500_000),
            read_length=runtime.get("read_length"),
            fragment_q95=runtime.get("fragment_q95"),
            fragment_q99=runtime.get("fragment_q99"),
            primary_contigs=_contig_list(runtime.get("primary_contigs")),
        ),
        output=OutputConfig(
            directory=_optional_path(output.get("directory"), "directory", base),
            log_level=output.get("log_level", "INFO"),
        ),
        refinement=RefinementConfig(
            breakpoints=refinement.get("breakpoints", False),
            complete_routing=refinement.get("complete_routing", False),
            pin_one_sided=refinement.get("pin_one_sided", False),
            call_annotations=refinement.get("call_annotations", False),
        ),
    )


def merge_config(
    base: EviTEAConfig,
    *,
    inputs: InputConfig | None = None,
    assets: PreparedAssetConfig | None = None,
    runtime: RuntimeConfig | None = None,
    output: OutputConfig | None = None,
    refinement: RefinementConfig | None = None,
) -> EviTEAConfig:
    """Apply legal explicit overrides without exposing science settings."""

    return EviTEAConfig(
        inputs=inputs or base.inputs,
        assets=assets or base.assets,
        runtime=runtime or base.runtime,
        output=output or base.output,
        science=base.science,
        refinement=refinement or base.refinement,
    )
