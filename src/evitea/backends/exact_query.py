"""Prepared-sample identity maps and exact-query manifests."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

from ..assets.context import (
    hash_relative_files,
    load_manifest,
    manifest_sha256,
    sha256_file,
    verify_relative_files,
    write_manifest,
)
from ..constants import ARCHITECTURE_ID
from ..contracts import ReadMate
from ..errors import AssetCompatibilityError, InputContractError
from .protocols import ReadIdentity

READ_IDENTITY_SCHEMA = "evitea.read-identity-map.v1"
SAMPLE_MANIFEST_SCHEMA = "evitea.prepared-sample.v1"
SAMPLE_MANIFEST_NAME = "sample_manifest.json"


def _nonnegative_int(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise InputContractError(f"{field} must be a non-negative integer")
    return value


@dataclass(frozen=True, slots=True)
class ReadIdentityMap:
    """Explicit backend-read to EviTEA pair/mate relation."""

    identities: Mapping[int, ReadIdentity]

    def __post_init__(self) -> None:
        normalized: dict[int, ReadIdentity] = {}
        for key, identity in self.identities.items():
            backend_id = _nonnegative_int(key, "backend_read_id")
            if identity.backend_read_id != backend_id:
                raise InputContractError("identity-map key and backend_read_id differ")
            _nonnegative_int(identity.pair_id, "pair_id")
            if not isinstance(identity.mate, ReadMate):
                raise InputContractError("mate must be ReadMate.READ1 or READ2")
            if backend_id in normalized:
                raise InputContractError("duplicate backend_read_id")
            normalized[backend_id] = identity
        object.__setattr__(
            self, "identities", MappingProxyType(dict(sorted(normalized.items())))
        )

    def resolve(self, backend_read_id: int) -> ReadIdentity:
        try:
            return self.identities[backend_read_id]
        except KeyError as exc:
            raise AssetCompatibilityError(
                f"backend read identity is absent from prepared sample: {backend_read_id}"
            ) from exc

    def to_document(self) -> dict[str, Any]:
        return {
            "schema_version": READ_IDENTITY_SCHEMA,
            "identities": [
                {
                    "backend_read_id": identity.backend_read_id,
                    "pair_id": identity.pair_id,
                    "mate": int(identity.mate),
                }
                for identity in self.identities.values()
            ],
        }


def make_read_identity_map(identities: Iterable[ReadIdentity]) -> ReadIdentityMap:
    materialized = tuple(identities)
    mapping = {identity.backend_read_id: identity for identity in materialized}
    if len(mapping) != len(materialized):
        raise InputContractError("duplicate backend_read_id")
    return ReadIdentityMap(mapping)


def write_read_identity_map(path: str | Path, identity_map: ReadIdentityMap) -> Path:
    return write_manifest(Path(path), identity_map.to_document())


def load_read_identity_map(path: str | Path) -> ReadIdentityMap:
    document = load_manifest(path)
    if document.get("schema_version") != READ_IDENTITY_SCHEMA:
        raise AssetCompatibilityError("unsupported read-identity-map schema")
    try:
        rows = document["identities"]
        if not isinstance(rows, list):
            raise TypeError("identities must be a list")
        identities = []
        for row in rows:
            if not isinstance(row, dict):
                raise TypeError("identity row must be an object")
            backend_read_id = row["backend_read_id"]
            pair_id = row["pair_id"]
            mate = row["mate"]
            if any(
                isinstance(value, bool) or not isinstance(value, int)
                for value in (backend_read_id, pair_id, mate)
            ):
                raise TypeError("identity integer fields must be JSON integers")
            identities.append(
                ReadIdentity(
                    backend_read_id=backend_read_id,
                    pair_id=pair_id,
                    mate=ReadMate(mate),
                )
            )
        return make_read_identity_map(identities)
    except (KeyError, TypeError, ValueError) as exc:
        raise AssetCompatibilityError("malformed read identity map") from exc


@dataclass(frozen=True, slots=True)
class PreparedSample:
    root: Path
    architecture_id: str
    sample_id: str
    pair_count: int
    read1_raw_sha256: str
    read2_raw_sha256: str
    identity_map_path: str
    identity_map_sha256: str
    backend_profile: Mapping[str, Any]
    file_hashes: Mapping[str, str]
    identity_sha256: str
    schema_version: str = SAMPLE_MANIFEST_SCHEMA

    def __post_init__(self) -> None:
        if not self.sample_id:
            raise InputContractError("sample_id must be non-empty text")
        _nonnegative_int(self.pair_count, "pair_count")
        object.__setattr__(self, "root", Path(self.root))
        object.__setattr__(
            self, "backend_profile", MappingProxyType(dict(self.backend_profile))
        )
        object.__setattr__(
            self, "file_hashes", MappingProxyType(dict(sorted(self.file_hashes.items())))
        )

    def identity_document(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "architecture_id": self.architecture_id,
            "sample_id": self.sample_id,
            "pair_count": self.pair_count,
            "read1_raw_sha256": self.read1_raw_sha256,
            "read2_raw_sha256": self.read2_raw_sha256,
            "identity_map_path": self.identity_map_path,
            "identity_map_sha256": self.identity_map_sha256,
            "backend_profile": dict(self.backend_profile),
            "file_hashes": dict(self.file_hashes),
        }

    def to_document(self) -> dict[str, Any]:
        return dict(self.identity_document(), identity_sha256=self.identity_sha256)


def prepare_sample_manifest(
    output_dir: str | Path,
    *,
    sample_id: str,
    pair_count: int,
    read1_path: str | Path,
    read2_path: str | Path,
    identity_map_relative_path: str | Path,
    backend_profile: Mapping[str, Any],
    index_relative_paths: Sequence[str | Path],
) -> PreparedSample:
    """Bind source reads, explicit identity map, and exact-query index files."""

    if not isinstance(sample_id, str) or not sample_id:
        raise InputContractError("sample_id must be non-empty text")
    _nonnegative_int(pair_count, "pair_count")
    if not index_relative_paths:
        raise InputContractError("prepared sample requires explicit index files")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / SAMPLE_MANIFEST_NAME
    if manifest_path.exists():
        raise InputContractError(f"prepared sample already exists: {manifest_path}")
    identity_relative = Path(identity_map_relative_path).as_posix()
    identity_path = root / identity_relative
    identity_map = load_read_identity_map(identity_path)
    _verify_identity_population(identity_map, pair_count)
    file_hashes = hash_relative_files(root, [identity_relative, *index_relative_paths])
    base = {
        "schema_version": SAMPLE_MANIFEST_SCHEMA,
        "architecture_id": ARCHITECTURE_ID,
        "sample_id": sample_id,
        "pair_count": pair_count,
        "read1_raw_sha256": sha256_file(read1_path),
        "read2_raw_sha256": sha256_file(read2_path),
        "identity_map_path": identity_relative,
        "identity_map_sha256": file_hashes[identity_relative],
        "backend_profile": dict(backend_profile),
        "file_hashes": file_hashes,
    }
    identity = manifest_sha256(base)
    write_manifest(manifest_path, dict(base, identity_sha256=identity))
    return load_prepared_sample(root)


def load_prepared_sample(path: str | Path) -> PreparedSample:
    candidate = Path(path)
    manifest_path = candidate if candidate.is_file() else candidate / SAMPLE_MANIFEST_NAME
    document = load_manifest(manifest_path)
    if document.get("schema_version") != SAMPLE_MANIFEST_SCHEMA:
        raise AssetCompatibilityError("unsupported prepared-sample schema")
    if document.get("architecture_id") != ARCHITECTURE_ID:
        raise AssetCompatibilityError("prepared sample architecture mismatch")
    identity = document.get("identity_sha256")
    base = dict(document)
    base.pop("identity_sha256", None)
    if not isinstance(identity, str) or manifest_sha256(base) != identity:
        raise AssetCompatibilityError("prepared-sample manifest identity mismatch")
    try:
        return PreparedSample(
            root=manifest_path.parent,
            architecture_id=str(document["architecture_id"]),
            sample_id=str(document["sample_id"]),
            pair_count=int(document["pair_count"]),
            read1_raw_sha256=str(document["read1_raw_sha256"]),
            read2_raw_sha256=str(document["read2_raw_sha256"]),
            identity_map_path=str(document["identity_map_path"]),
            identity_map_sha256=str(document["identity_map_sha256"]),
            backend_profile=document["backend_profile"],
            file_hashes=document["file_hashes"],
            identity_sha256=identity,
        )
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise AssetCompatibilityError("malformed prepared-sample manifest") from exc


def verify_prepared_sample(asset: PreparedSample) -> ReadIdentityMap:
    if asset.architecture_id != ARCHITECTURE_ID:
        raise AssetCompatibilityError("prepared sample architecture mismatch")
    if manifest_sha256(asset.identity_document()) != asset.identity_sha256:
        raise AssetCompatibilityError("prepared-sample in-memory identity mismatch")
    verify_relative_files(asset.root, asset.file_hashes)
    identity_path = asset.root / asset.identity_map_path
    if sha256_file(identity_path) != asset.identity_map_sha256:
        raise AssetCompatibilityError("prepared-sample identity-map hash mismatch")
    identity_map = load_read_identity_map(identity_path)
    _verify_identity_population(identity_map, asset.pair_count)
    return identity_map


def _verify_identity_population(identity_map: ReadIdentityMap, pair_count: int) -> None:
    pair_mates: dict[int, set[ReadMate]] = {}
    for identity in identity_map.identities.values():
        pair_mates.setdefault(identity.pair_id, set()).add(identity.mate)
    expected_pairs = set(range(pair_count))
    if (
        len(identity_map.identities) != pair_count * 2
        or set(pair_mates) != expected_pairs
        or any(
            mates != {ReadMate.READ1, ReadMate.READ2}
            for mates in pair_mates.values()
        )
    ):
        raise AssetCompatibilityError(
            "identity map must contain exactly READ1 and READ2 for ordinal-zero pair IDs"
        )
