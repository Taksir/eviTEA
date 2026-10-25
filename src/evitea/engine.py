"""Default scientific engine: the library layer followed by the extension layers.

One call runs, on one prepared sample msBWT::

    library layer  -> bilateral -> asymmetric -> full memory -> internal 21-mer

``evitea.incumbent.run_incumbent_caller`` derives the repeat-aware library-layer
callset; ``evitea.extension_engine.run_extension`` appends the four extension
layers to its exact bytes.  The final native callset has 838 rows for AKA017
and 194 rows for ISO1 on the development samples.

Profiles
--------
Only Bowtie2 query names (they seed Bowtie2's per-read tie-breaking) and the
sample labels in extension event identifiers differ between the default
profile and a byte-for-byte reproduction of a development sample:

=================  ==========================  ===================  ============
profile            paired-end host names       wide-anchor names    labels
=================  ==========================  ===================  ============
production         ``P{pid}``                  ``pair_{pid}``       from sample
AKA017 (repro)     ``p{pid}|h{mate}``          ``pair_{pid}``       fixed
ISO1 (repro)       ``P{pid}``                  ``{pid}``            fixed
=================  ==========================  ===================  ============

A reproduction profile accepts only its own sample id.

No scientific threshold is a profile setting.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

from .anchor_pe.alignments import PRODUCTION_HOST_READ_NAMES
from .anchor_pe.engine import AnchorPEInputs
from .errors import InputContractError
from .extension_engine import ExtensionRun, run_extension
from .contigs import DEVELOPMENT_PRIMARY_CONTIGS
from .sample_geometry import (
    DEFAULT_READ_LENGTH,
    FragmentSizes,
    estimate_fragment_sizes,
    resolve_fragment_sizes,
    validate_fragment_size,
    validate_read_length,
)
from .incumbent import IncumbentInputs, IncumbentResult, MsbwtSample, run_incumbent_caller
from .shadow_layers import ShadowCallset, ShadowLabels, UnresolvedShadowEvent

__all__ = [
    "PRODUCTION_PROFILE",
    "REPRODUCTION_PROFILES",
    "EngineInputs",
    "EngineResult",
    "RunProfile",
    "default_bowtie2_build",
    "resolve_profile",
    "run_engine",
]

PRODUCTION_PROFILE = "production"


@dataclass(frozen=True, slots=True)
class RunProfile:
    name: str
    host_read_names: str
    labels: ShadowLabels
    reproduction_only: bool
    rescue_read_names: str = "pair"


REPRODUCTION_PROFILES: Mapping[str, RunProfile] = MappingProxyType(
    {
        "AKA017": RunProfile(
            "AKA017", "aka-historical", ShadowLabels("", "AKA", "AKA", "AKA"), True
        ),
        "ISO1": RunProfile(
            "ISO1", PRODUCTION_HOST_READ_NAMES, ShadowLabels("ISO1_", "ISO1", "ISO", "ISO"), True,
            rescue_read_names="numeric",
        ),
    }
)


def resolve_profile(name: str | None, sample_id: str) -> RunProfile:
    """Return the default profile (``production``) or a named reproduction profile."""

    if not sample_id:
        raise InputContractError("sample_id must be non-empty text")
    if name in (None, "", PRODUCTION_PROFILE):
        return RunProfile(
            PRODUCTION_PROFILE, PRODUCTION_HOST_READ_NAMES, ShadowLabels.for_sample(sample_id), False
        )
    try:
        profile = REPRODUCTION_PROFILES[name]
    except KeyError as exc:
        raise InputContractError(
            f"unknown run profile {name!r}; use {PRODUCTION_PROFILE!r} or one of "
            f"{sorted(REPRODUCTION_PROFILES)}"
        ) from exc
    if sample_id != profile.name:
        raise InputContractError(
            f"reproduction profile {profile.name} reproduces only sample {profile.name}, "
            f"not {sample_id!r}"
        )
    return profile


@dataclass(frozen=True, slots=True)
class EngineInputs:
    sample_id: str
    msbwt_index: Path
    te_fasta: Path
    te_family_map: Path
    reference_fasta: Path
    repeatmasker_bed: Path
    host_index_prefix: Path
    terminal_codebook: Path
    internal_codebook: Path
    bowtie2: Path
    bowtie2_build: Path | None = None
    te_index_prefix: Path | None = None
    profile: str = PRODUCTION_PROFILE
    threads: int = 1
    workers: int = 1
    edlib_python: Path | None = None
    complete_routing: bool = False
    max_rescue_interval_mass: int = 1_500_000
    read_length: int = DEFAULT_READ_LENGTH
    # where the read length came from (see ``sample_geometry.resolve_read_length``); not part of the identity
    read_length_source: str = "default"
    # fragment-size 95th / 99th percentile overrides; None = estimated from the sample (``sample_geometry``)
    fragment_q95: float | None = None
    fragment_q99: int | None = None
    # contigs of the full-memory proposals (``evitea.contigs``) and where the list came from
    primary_contigs: tuple[str, ...] = DEVELOPMENT_PRIMARY_CONTIGS
    primary_contigs_source: str = "default"

    def __post_init__(self) -> None:
        validate_read_length(self.read_length)
        if self.fragment_q95 is not None:
            validate_fragment_size(self.fragment_q95, "fragment_q95")
        if self.fragment_q99 is not None:
            validate_fragment_size(self.fragment_q99, "fragment_q99")
            if not isinstance(self.fragment_q99, int):
                raise InputContractError("fragment_q99 must be an integer")
        for name in ("threads", "workers"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise InputContractError(f"{name} must be a positive integer")
        resolve_profile(self.profile, self.sample_id)

    @property
    def run_profile(self) -> RunProfile:
        return resolve_profile(self.profile, self.sample_id)

    def identity(self) -> dict[str, object]:
        """Return the scientific identity (inputs and profile), without threads and workers."""

        profile = self.run_profile
        return {
            "sample_id": self.sample_id,
            "profile": profile.name,
            "host_read_names": profile.host_read_names,
            "rescue_read_names": profile.rescue_read_names,
            "labels": [
                profile.labels.bilateral_prefix,
                profile.labels.asymmetric,
                profile.labels.full_memory,
                profile.labels.internal21,
            ],
            **{
                name: None if getattr(self, name) is None else str(getattr(self, name))
                for name in (
                    "msbwt_index", "te_fasta", "te_family_map", "reference_fasta",
                    "repeatmasker_bed", "host_index_prefix", "terminal_codebook",
                    "internal_codebook", "te_index_prefix",
                )
            },
            # recorded only when not at the default, so default identities do not change
            **({"complete_routing": True} if self.complete_routing else {}),
            **(
                {"max_rescue_interval_mass": self.max_rescue_interval_mass}
                if self.max_rescue_interval_mass != 1_500_000
                else {}
            ),
            **({"read_length": self.read_length} if self.read_length != DEFAULT_READ_LENGTH else {}),
            **({"fragment_q95": self.fragment_q95} if self.fragment_q95 is not None else {}),
            **({"fragment_q99": self.fragment_q99} if self.fragment_q99 is not None else {}),
            **({"primary_contigs": list(self.primary_contigs)}
               if tuple(self.primary_contigs) != DEVELOPMENT_PRIMARY_CONTIGS else {}),
        }

    def run_parameters(self) -> dict[str, object]:
        """Return the sample-derived values used by this run (recorded in the run manifest)."""

        return {"read_length": self.read_length, "read_length_source": self.read_length_source,
                "primary_contigs": list(self.primary_contigs), "primary_contigs_source": self.primary_contigs_source}


@dataclass(slots=True)
class EngineResult:
    sample_id: str
    profile: RunProfile
    incumbent: IncumbentResult
    extension: ExtensionRun
    fragment: FragmentSizes | None = None

    @property
    def incumbent_bytes(self) -> bytes:
        return self.incumbent.legacy_tsv()

    @property
    def final(self) -> ShadowCallset:
        return self.extension.layers[-1].callset

    @property
    def unresolved(self) -> tuple[UnresolvedShadowEvent, ...]:
        return tuple(row for layer in self.extension.layers for row in layer.unresolved)

    @property
    def layer_rows(self) -> dict[str, int]:
        rows = {"incumbent": ShadowCallset.from_bytes(self.incumbent_bytes).row_count}
        for layer in self.extension.layers:
            rows[layer.layer.lower()] = layer.callset.row_count
        return rows

    @property
    def layer_sha256(self) -> dict[str, str]:
        """Return the SHA-256 of the callset after each layer."""

        hashes = {"incumbent": hashlib.sha256(self.incumbent_bytes).hexdigest()}
        for layer in self.extension.layers:
            hashes[layer.layer.lower()] = hashlib.sha256(layer.callset.data).hexdigest()
        return hashes

    @property
    def counts(self) -> dict[str, dict[str, int]]:
        return {
            "incumbent": dict(self.incumbent.counts),
            "anchor_pe": dict(self.extension.anchor_pe.counts),
            "full_memory": dict(self.extension.full_memory.counts),
            "internal21": dict(self.extension.internal21.counts),
        }


def default_bowtie2_build(bowtie2: Path) -> Path | None:
    """Return ``bowtie2-build`` next to ``bowtie2`` (same suffix), else on PATH."""

    candidate = Path(bowtie2)
    sibling = candidate.with_name("bowtie2-build" + "".join(candidate.suffixes))
    if sibling.is_file():
        return sibling
    found = shutil.which("bowtie2-build")
    return Path(found) if found else None


def _fragment_sizes(inputs: EngineInputs, sample, work: Path, say: Callable[[str], None]) -> FragmentSizes:
    """Return the fragment-size percentiles of this run: options, else the sample's estimate (saved for resume)."""

    saved = work / "fragment_sizes.json"
    key = {"msbwt_index": str(inputs.msbwt_index), "host_index_prefix": str(inputs.host_index_prefix),
           "fragment_q95": inputs.fragment_q95, "fragment_q99": inputs.fragment_q99}
    if saved.is_file():
        record = json.loads(saved.read_text(encoding="utf-8"))
        if record.get("key") == key:
            return FragmentSizes(**record["value"])

    def estimate() -> FragmentSizes:
        say("fragment size: mapping evenly spaced read pairs to the reference")
        return estimate_fragment_sizes(
            sample, bowtie2=inputs.bowtie2, host_index_prefix=inputs.host_index_prefix,
            work_directory=work / "fragment_sizes", threads=inputs.threads,
        )

    fragment = resolve_fragment_sizes(inputs.fragment_q95, inputs.fragment_q99, estimate, warn=say)
    value = {"q95": fragment.q95, "q99": fragment.q99, "pairs": fragment.pairs, "source": fragment.source}
    saved.write_text(json.dumps({"key": key, "value": value}, sort_keys=True), encoding="utf-8")
    return fragment


def run_engine(
    inputs: EngineInputs,
    work_directory: str | Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> EngineResult:
    """Run the library layer and the extension layers; checkpoints go under ``work_directory``."""

    profile = inputs.run_profile
    work = Path(work_directory)
    work.mkdir(parents=True, exist_ok=True)
    bowtie2_build = inputs.bowtie2_build or default_bowtie2_build(inputs.bowtie2)
    if inputs.te_index_prefix is None and bowtie2_build is None:
        raise InputContractError("bowtie2-build is required to index the TE library")
    sample = MsbwtSample(inputs.msbwt_index)

    def say(message: str) -> None:
        if progress is not None:
            progress(message)

    fragment = _fragment_sizes(inputs, sample, work, say)
    say(f"read length {inputs.read_length} bp ({inputs.read_length_source}); fragment size 95th / 99th percentile "
        f"{fragment.q95:g} / {fragment.q99} bp ({fragment.source}, {fragment.pairs} pairs)")
    say(f"incumbent: sample {inputs.sample_id}, profile {profile.name}")
    incumbent = run_incumbent_caller(
        IncumbentInputs(
            sample_id=inputs.sample_id,
            te_fasta=inputs.te_fasta,
            te_family_map=inputs.te_family_map,
            reference_fasta=inputs.reference_fasta,
            host_index_prefix=inputs.host_index_prefix,
            bowtie2=inputs.bowtie2,
            bowtie2_build=bowtie2_build,
            te_index_prefix=inputs.te_index_prefix,
            threads=inputs.threads,
            workers=inputs.workers,
            rescue_read_names=profile.rescue_read_names,
            complete_routing=inputs.complete_routing,
            max_rescue_interval_mass=inputs.max_rescue_interval_mass,
            read_length=inputs.read_length,
            fragment_q95=fragment.q95,
            fragment_q99=fragment.q99,
        ),
        sample,
        work / "incumbent",
        progress=progress,
    )
    say(f"incumbent calls: {len(incumbent.calls)}; v3 extension")
    extension = run_extension(
        AnchorPEInputs(
            sample_id=inputs.sample_id,
            safe_codebook=inputs.terminal_codebook,
            te_fasta=inputs.te_fasta,
            te_family_map=inputs.te_family_map,
            repeatmasker_bed=inputs.repeatmasker_bed,
            host_index_prefix=inputs.host_index_prefix,
            bowtie2=inputs.bowtie2,
            host_read_names=profile.host_read_names,
            threads=inputs.threads,
            edlib_python=inputs.edlib_python,
            fragment_q99=fragment.q99,
            primary_contigs=tuple(inputs.primary_contigs),
        ),
        inputs.internal_codebook,
        sample,
        incumbent.legacy_tsv(),
        profile.labels,
        work / "extension",
        progress=progress,
    )
    result = EngineResult(inputs.sample_id, profile, incumbent, extension, fragment)
    say(f"layers: {result.layer_rows}")
    say(f"layer sha256: {result.layer_sha256}")
    return result
