"""Prepared-reference manifests with exact source and index identities."""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

from ..constants import ARCHITECTURE_ID
from ..contracts import ReferenceInput, RepeatContext
from ..errors import AssetCompatibilityError, InputContractError
from .context import (
    hash_relative_files,
    load_manifest,
    manifest_sha256,
    sha256_file,
    verify_relative_files,
    write_manifest,
)

REFERENCE_MANIFEST_SCHEMA = "evitea.prepared-reference.v1"
REFERENCE_MANIFEST_NAME = "reference_manifest.json"


@dataclass(frozen=True, slots=True)
class PreparedReference:
    root: Path
    architecture_id: str
    reference_raw_sha256: str
    reference_semantic_sha256: str
    contigs: tuple[tuple[str, int], ...]
    repeat_raw_sha256: str | None
    mapper_profile: Mapping[str, Any]
    file_hashes: Mapping[str, str]
    identity_sha256: str
    schema_version: str = REFERENCE_MANIFEST_SCHEMA

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", Path(self.root))
        object.__setattr__(
            self, "mapper_profile", MappingProxyType(dict(self.mapper_profile))
        )
        object.__setattr__(
            self, "file_hashes", MappingProxyType(dict(sorted(self.file_hashes.items())))
        )

    def identity_document(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "architecture_id": self.architecture_id,
            "reference_raw_sha256": self.reference_raw_sha256,
            "reference_semantic_sha256": self.reference_semantic_sha256,
            "contigs": [
                {"contig_id": contig_id, "length": length}
                for contig_id, length in self.contigs
            ],
            "repeat_raw_sha256": self.repeat_raw_sha256,
            "mapper_profile": dict(self.mapper_profile),
            "file_hashes": dict(self.file_hashes),
        }

    def _source(self, prefix: str) -> Path | None:
        matches = [path for path in self.file_hashes if path.startswith(f"sources/{prefix}")]
        if len(matches) > 1:
            raise AssetCompatibilityError(f"prepared reference binds several {prefix} sources")
        return self.root / matches[0] if matches else None

    @property
    def reference_fasta_path(self) -> Path:
        path = self._source("reference")
        if path is None:
            raise AssetCompatibilityError("prepared reference binds no reference FASTA")
        return path

    @property
    def repeat_annotation_path(self) -> Path:
        path = self._source("repeat_annotation")
        if path is None:
            raise AssetCompatibilityError(
                "prepared reference binds no RepeatMasker annotation; re-run "
                "prepare-reference with --repeat-annotation"
            )
        return path

    @property
    def bowtie2_index_prefix(self) -> Path:
        prefix = self.mapper_profile.get("index_prefix")
        if not isinstance(prefix, str) or not prefix:
            raise AssetCompatibilityError("prepared reference mapper profile has no index_prefix")
        return self.root / prefix

    def to_document(self) -> dict[str, Any]:
        document = self.identity_document()
        document["identity_sha256"] = self.identity_sha256
        return document


def _source_copy_name(prefix: str, path: Path) -> str:
    suffix = "".join(path.suffixes)
    return f"sources/{prefix}{suffix}"


def prepare_reference(
    reference: ReferenceInput,
    output_dir: str | Path,
    *,
    mapper_profile: Mapping[str, Any],
    index_relative_paths: Sequence[str | Path],
    repeat_context: RepeatContext | None = None,
) -> PreparedReference:
    """Copy source inputs and bind already-built mapper index files."""

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / REFERENCE_MANIFEST_NAME
    if manifest_path.exists():
        raise InputContractError(f"prepared reference already exists: {manifest_path}")
    if not index_relative_paths:
        raise InputContractError("prepared reference requires explicit mapper index files")
    if sha256_file(reference.path) != reference.raw_sha256:
        raise AssetCompatibilityError("reference source changed after validation")
    if repeat_context is not None and (
        repeat_context.reference_semantic_sha256 != reference.semantic_sha256
    ):
        raise AssetCompatibilityError("repeat context is bound to another reference")
    if repeat_context is not None and sha256_file(repeat_context.path) != repeat_context.raw_sha256:
        raise AssetCompatibilityError("repeat annotation changed after validation")

    reference_relative = _source_copy_name("reference", reference.path)
    reference_copy = root / reference_relative
    reference_copy.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(reference.path, reference_copy)
    bound_paths: list[str | Path] = [reference_relative, *index_relative_paths]

    repeat_hash: str | None = None
    if repeat_context is not None:
        repeat_relative = _source_copy_name("repeat_annotation", repeat_context.path)
        shutil.copyfile(repeat_context.path, root / repeat_relative)
        bound_paths.append(repeat_relative)
        repeat_hash = repeat_context.raw_sha256

    file_hashes = hash_relative_files(root, bound_paths)
    base = {
        "schema_version": REFERENCE_MANIFEST_SCHEMA,
        "architecture_id": ARCHITECTURE_ID,
        "reference_raw_sha256": reference.raw_sha256,
        "reference_semantic_sha256": reference.semantic_sha256,
        "contigs": [
            {"contig_id": record.record_id, "length": len(record.sequence)}
            for record in reference.records
        ],
        "repeat_raw_sha256": repeat_hash,
        "mapper_profile": dict(mapper_profile),
        "file_hashes": file_hashes,
    }
    identity = manifest_sha256(base)
    document = dict(base, identity_sha256=identity)
    write_manifest(manifest_path, document)
    return load_prepared_reference(root)


def load_prepared_reference(path: str | Path) -> PreparedReference:
    """Load and structurally verify a prepared-reference manifest."""

    candidate = Path(path)
    manifest_path = candidate if candidate.is_file() else candidate / REFERENCE_MANIFEST_NAME
    document = load_manifest(manifest_path)
    if document.get("schema_version") != REFERENCE_MANIFEST_SCHEMA:
        raise AssetCompatibilityError("unsupported prepared-reference schema")
    if document.get("architecture_id") != ARCHITECTURE_ID:
        raise AssetCompatibilityError("prepared reference architecture mismatch")
    identity = document.get("identity_sha256")
    base = dict(document)
    base.pop("identity_sha256", None)
    if not isinstance(identity, str) or manifest_sha256(base) != identity:
        raise AssetCompatibilityError("prepared-reference manifest identity mismatch")
    try:
        contigs = tuple(
            (str(record["contig_id"]), int(record["length"]))
            for record in document["contigs"]
        )
        return PreparedReference(
            root=manifest_path.parent,
            architecture_id=str(document["architecture_id"]),
            reference_raw_sha256=str(document["reference_raw_sha256"]),
            reference_semantic_sha256=str(document["reference_semantic_sha256"]),
            contigs=contigs,
            repeat_raw_sha256=(
                None
                if document.get("repeat_raw_sha256") is None
                else str(document["repeat_raw_sha256"])
            ),
            mapper_profile=document["mapper_profile"],
            file_hashes=document["file_hashes"],
            identity_sha256=identity,
        )
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise AssetCompatibilityError("malformed prepared-reference manifest") from exc


def verify_reference_asset(asset: PreparedReference) -> None:
    """Verify manifest identity and every bound source/index file."""

    if asset.architecture_id != ARCHITECTURE_ID:
        raise AssetCompatibilityError("prepared reference architecture mismatch")
    if manifest_sha256(asset.identity_document()) != asset.identity_sha256:
        raise AssetCompatibilityError("prepared-reference in-memory identity mismatch")
    verify_relative_files(asset.root, asset.file_hashes)
