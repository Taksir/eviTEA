"""Runtime provenance identities and deterministic output finalization."""

from __future__ import annotations

import hashlib
import platform
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ._version import __version__
from .constants import ARCHITECTURE_ID
from .contracts import PairIdentityMode, ProvenanceRecord, stable_to_primitive
from .errors import AssetCompatibilityError, InputContractError

RUN_MANIFEST_SCHEMA = "evitea.run-manifest.v1"


def hash_file(path: str | Path, chunk_size: int = 8 << 20) -> str:
    candidate = Path(path)
    digest = hashlib.sha256()
    try:
        with candidate.open("rb") as handle:
            while chunk := handle.read(chunk_size):
                digest.update(chunk)
    except OSError as exc:
        raise InputContractError(f"cannot hash file: {candidate}") from exc
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class RunManifestContext:
    normalized_config: Mapping[str, Any]
    command_line: tuple[str, ...]
    source_identity: Mapping[str, Any]
    input_assets: Mapping[str, Mapping[str, Any]]
    prepared_identities: Mapping[str, str]
    backend_versions: Mapping[str, str]
    resolved_threads: int
    stage_history: tuple[Mapping[str, Any], ...] = ()
    censor_counts: Mapping[str, int] = field(default_factory=dict)
    warning_counts: Mapping[str, int] = field(default_factory=dict)
    started_at: str | None = None
    ended_at: str | None = None
    pair_identity_mode: PairIdentityMode = PairIdentityMode.EVITEA_ORDINAL0
    run_parameters: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "censor_counts", dict(self.censor_counts))
        object.__setattr__(self, "warning_counts", dict(self.warning_counts))
        if self.resolved_threads <= 0:
            raise InputContractError("resolved thread count must be positive")


def file_identity(path: str | Path) -> dict[str, Any]:
    candidate = Path(path).resolve()
    try:
        size = candidate.stat().st_size
    except OSError as exc:
        raise InputContractError(f"cannot inspect file: {candidate}") from exc
    return {"path": str(candidate), "size": size, "sha256": hash_file(candidate)}


def build_run_manifest(context: RunManifestContext) -> ProvenanceRecord:
    values = {
        "evitea_version": __version__,
        "source_identity": context.source_identity,
        "normalized_config": context.normalized_config,
        "command_line": context.command_line,
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "python_version": sys.version.split()[0],
        "backend_versions": context.backend_versions,
        "input_assets": context.input_assets,
        "prepared_identities": context.prepared_identities,
        "pair_identity_mode": context.pair_identity_mode,
        "pairing_contract_version": "evitea.paired-fastq.v1",
        "resolved_threads": context.resolved_threads,
        "run_parameters": dict(context.run_parameters),
        "stage_history": context.stage_history,
        "censor_counts": context.censor_counts,
        "warning_counts": context.warning_counts,
        "started_at": context.started_at,
        "ended_at": context.ended_at,
        "output_hashes": {},
        "output_counts": {},
    }
    return ProvenanceRecord(
        schema_version=RUN_MANIFEST_SCHEMA,
        architecture_id=ARCHITECTURE_ID,
        values=stable_to_primitive(values),
    )


def verify_prepared_asset(
    path: str | Path,
    expected_sha256: str,
) -> None:
    actual = hash_file(path)
    if actual != expected_sha256:
        raise AssetCompatibilityError(
            f"prepared asset hash mismatch: expected {expected_sha256}, got {actual}"
        )


def finalize_output_hashes(
    manifest: ProvenanceRecord,
    output_directory: str | Path,
    filenames: Sequence[str],
    *,
    output_counts: Mapping[str, int] | None = None,
) -> ProvenanceRecord:
    """Record the hash of each output file in the manifest, which never hashes itself."""

    if "run_manifest.json" in filenames:
        raise InputContractError("run manifest cannot recursively hash itself")
    root = Path(output_directory)
    hashes = {name: hash_file(root / name) for name in sorted(set(filenames))}
    values = dict(manifest.values)
    values["output_hashes"] = hashes
    values["output_counts"] = dict(sorted((output_counts or {}).items()))
    return ProvenanceRecord(
        schema_version=manifest.schema_version,
        architecture_id=manifest.architecture_id,
        values=values,
    )
