"""Prepared TE-library manifests and the v3 detection languages they carry.

A prepared TE asset binds the TE consensus library and its family map to one
prepared reference.  Since schema v3 it also carries the two languages the
default engine consumes:

* ``terminal_final`` - the final safe terminal 21-mer codebook (AnchorPE and
  full-memory branches);
* ``internal_only`` - the internal-only safe 21-mer codebook (INTERNAL21).

Each language records its SHA-256 and its ``source``: ``built:<builder>``
(``evitea.languages``, from the TE library, family map, reference and
RepeatMasker annotation; every intermediate table is kept under
``languages/audit/``) or ``supplied`` (a file given to ``prepare-te``).
"""

from __future__ import annotations

import shutil
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from ..constants import ARCHITECTURE_ID
from ..contracts import ProbeClass, TELibraryInput, TEProbe
from ..errors import AssetCompatibilityError, InputContractError
from ..families import FamilyMap, validate_family_map
from ..evidence.frozen_candidates import build_frozen_candidate_codebooks
from ..inputs import load_repeat_annotation, validate_reference_fasta
from .context import (
    hash_relative_files,
    load_manifest,
    manifest_sha256,
    sha256_file,
    verify_relative_files,
    write_manifest,
)
from .reference import PreparedReference, verify_reference_asset

TE_MANIFEST_SCHEMA = "evitea.prepared-te.v3"
LANGUAGE_KEYS = ("terminal_final", "internal_only")
LANGUAGE_PATHS = {
    "terminal_final": "languages/final_safe_terminal_codebook.tsv",
    "internal_only": "languages/safe_internal21_internal_only_codebook.tsv",
}
TE_MANIFEST_NAME = "te_manifest.json"

if TYPE_CHECKING:
    from ..evidence.codebook import SafeProbeCandidate


@dataclass(frozen=True, slots=True)
class PreparedTEAssets:
    root: Path
    architecture_id: str
    te_raw_sha256: str
    te_semantic_sha256: str
    family_assignments: Mapping[str, str]
    reference_identity_sha256: str
    file_hashes: Mapping[str, str]
    resources: Mapping[str, Mapping[str, str]]
    codebooks: Mapping[str, Mapping[str, str]]
    identity_sha256: str
    schema_version: str = TE_MANIFEST_SCHEMA
    languages: Mapping[str, Mapping[str, str]] = MappingProxyType({})

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", Path(self.root))
        object.__setattr__(
            self,
            "family_assignments",
            MappingProxyType(dict(sorted(self.family_assignments.items()))),
        )
        object.__setattr__(
            self, "file_hashes", MappingProxyType(dict(sorted(self.file_hashes.items())))
        )
        object.__setattr__(
            self,
            "resources",
            MappingProxyType(
                {
                    key: MappingProxyType(dict(sorted(value.items())))
                    for key, value in sorted(self.resources.items())
                }
            ),
        )
        object.__setattr__(
            self,
            "codebooks",
            MappingProxyType(
                {
                    key: MappingProxyType(dict(sorted(value.items())))
                    for key, value in sorted(self.codebooks.items())
                }
            ),
        )
        object.__setattr__(
            self,
            "languages",
            MappingProxyType(
                {
                    key: MappingProxyType(dict(sorted(value.items())))
                    for key, value in sorted(self.languages.items())
                }
            ),
        )

    def language_path(self, key: str) -> Path:
        """Absolute path of a bound language (fails closed when absent)."""

        try:
            return self.root / self.languages[key]["path"]
        except KeyError as exc:
            raise AssetCompatibilityError(
                f"prepared TE assets carry no {key!r} language; re-run prepare-te with "
                "--terminal-codebook and --internal-codebook"
            ) from exc

    def source_path(self, prefix: str) -> Path:
        """Absolute path of the one bound ``sources/<prefix>*`` file."""

        matches = [path for path in self.file_hashes if path.startswith(f"sources/{prefix}")]
        if len(matches) != 1:
            raise AssetCompatibilityError(f"prepared TE assets must bind exactly one {prefix} source")
        return self.root / matches[0]

    def identity_document(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "architecture_id": self.architecture_id,
            "te_raw_sha256": self.te_raw_sha256,
            "te_semantic_sha256": self.te_semantic_sha256,
            "family_assignments": dict(self.family_assignments),
            "reference_identity_sha256": self.reference_identity_sha256,
            "file_hashes": dict(self.file_hashes),
            "resources": {key: dict(value) for key, value in self.resources.items()},
            "codebooks": {key: dict(value) for key, value in self.codebooks.items()},
            "languages": {key: dict(value) for key, value in self.languages.items()},
        }

    def to_document(self) -> dict[str, Any]:
        return dict(self.identity_document(), identity_sha256=self.identity_sha256)


def prepare_te_assets(
    te_library: TELibraryInput,
    family_map: FamilyMap,
    family_map_path: str | Path,
    prepared_reference: PreparedReference,
    output_dir: str | Path,
    *,
    asset_relative_paths: Sequence[str | Path] = (),
    g1a_archive_path: str | Path | None = None,
    family_crosswalk_path: str | Path | None = None,
    terminal_codebook_path: str | Path | None = None,
    internal_codebook_path: str | Path | None = None,
    allow_unassigned_records: bool = False,
    build_languages: bool = False,
    progress=None,
    language_read_length: int = 125,
) -> PreparedTEAssets:
    """Copy TE/family sources and bind them to one prepared reference."""

    verify_reference_asset(prepared_reference)
    validated_map = validate_family_map(
        family_map,
        (record.record_id for record in te_library.records),
        allow_unassigned=allow_unassigned_records,
    )
    source_map_path = Path(family_map_path)
    if not source_map_path.is_file():
        raise InputContractError(f"family map does not exist: {source_map_path}")
    if sha256_file(te_library.path) != te_library.raw_sha256:
        raise AssetCompatibilityError("TE source changed after validation")
    if (terminal_codebook_path is None) != (internal_codebook_path is None):
        raise InputContractError(
            "the terminal and internal-only codebooks must be supplied together"
        )
    if build_languages and terminal_codebook_path is not None:
        raise InputContractError("either build the v3 languages or supply them, not both")
    if (g1a_archive_path is None) != (family_crosswalk_path is None):
        raise InputContractError(
            "frozen G1A archive and FBte family crosswalk must be supplied together"
        )
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / TE_MANIFEST_NAME
    if manifest_path.exists():
        raise InputContractError(f"prepared TE assets already exist: {manifest_path}")

    te_relative = "sources/te_library" + "".join(te_library.path.suffixes)
    map_relative = "sources/family_map" + "".join(source_map_path.suffixes)
    (root / "sources").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(te_library.path, root / te_relative)
    shutil.copyfile(source_map_path, root / map_relative)
    resources: dict[str, dict[str, str]] = {}
    codebooks: dict[str, dict[str, str]] = {}
    additional_paths: list[str] = []
    if g1a_archive_path is not None and family_crosswalk_path is not None:
        archive_source = Path(g1a_archive_path)
        crosswalk_source = Path(family_crosswalk_path)
        archive_relative = "sources/G1A_FROZEN_RESULTS_FOR_CODEBOOK.zip"
        crosswalk_relative = "sources/fbte_to_flybase_family_FROZEN.tsv"
        shutil.copyfile(archive_source, root / archive_relative)
        shutil.copyfile(crosswalk_source, root / crosswalk_relative)

        reference_sources = [
            path
            for path in prepared_reference.file_hashes
            if path.startswith("sources/reference")
        ]
        repeat_sources = [
            path
            for path in prepared_reference.file_hashes
            if path.startswith("sources/repeat_annotation")
        ]
        if len(reference_sources) != 1 or len(repeat_sources) != 1:
            raise AssetCompatibilityError(
                "frozen G1A preparation requires one bound reference and RepeatMasker source"
            )
        reference_input = validate_reference_fasta(
            prepared_reference.root / reference_sources[0]
        )
        if (
            reference_input.raw_sha256 != prepared_reference.reference_raw_sha256
            or reference_input.semantic_sha256 != prepared_reference.reference_semantic_sha256
        ):
            raise AssetCompatibilityError("prepared reference source identity mismatch")
        repeat_context = load_repeat_annotation(
            prepared_reference.root / repeat_sources[0], reference_input
        )
        if repeat_context.raw_sha256 != prepared_reference.repeat_raw_sha256:
            raise AssetCompatibilityError("prepared RepeatMasker source identity mismatch")
        built = build_frozen_candidate_codebooks(
            archive_source,
            crosswalk_source,
            te_library,
            validated_map,
            reference_input,
            repeat_context,
        )
        terminal_relative = "codebooks/terminal_safe_21mers.tsv"
        internal_relative = "codebooks/internal_safe_21mers.tsv"
        internal_only_relative = "codebooks/internal_only_safe_21mers.tsv"
        (root / "codebooks").mkdir(parents=True, exist_ok=True)
        (root / terminal_relative).write_bytes(built.terminal_tsv)
        (root / internal_relative).write_bytes(built.internal_tsv)
        (root / internal_only_relative).write_bytes(built.internal_only_tsv)
        resources = {
            "g1a_archive": {
                "path": archive_relative,
                "sha256": sha256_file(archive_source),
            },
            "family_crosswalk": {
                "path": crosswalk_relative,
                "sha256": sha256_file(crosswalk_source),
            },
        }
        codebooks = {
            "terminal": {
                "path": terminal_relative,
                "sha256": sha256_file(root / terminal_relative),
            },
            "internal_full": {
                "path": internal_relative,
                "sha256": sha256_file(root / internal_relative),
            },
            "internal_only": {
                "path": internal_only_relative,
                "sha256": sha256_file(root / internal_only_relative),
            },
        }
        additional_paths.extend(
            [
                archive_relative,
                crosswalk_relative,
                terminal_relative,
                internal_relative,
                internal_only_relative,
            ]
        )
    languages: dict[str, dict[str, str]] = {}
    if build_languages:
        from ..languages import BUILDER_VERSION
        from ..languages import build_languages as run_language_builders

        built = run_language_builders(
            te_fasta=root / te_relative,
            family_map=root / map_relative,
            reference_fasta=prepared_reference.reference_fasta_path,
            repeatmasker_bed=prepared_reference.repeat_annotation_path,
            progress=progress,
            read_length=language_read_length,
        )
        (root / "languages" / "audit").mkdir(parents=True, exist_ok=True)
        (root / LANGUAGE_PATHS["terminal_final"]).write_bytes(built.terminal_codebook)
        (root / LANGUAGE_PATHS["internal_only"]).write_bytes(built.internal_only_codebook)
        for name, data in built.audit_tables().items():
            relative = f"languages/audit/{name}"
            (root / relative).write_bytes(data)
            additional_paths.append(relative)
        validate_language_files(root / LANGUAGE_PATHS["terminal_final"], root / LANGUAGE_PATHS["internal_only"])
        for key in LANGUAGE_KEYS:
            languages[key] = {
                "path": LANGUAGE_PATHS[key],
                "sha256": sha256_file(root / LANGUAGE_PATHS[key]),
                "source": f"built:{BUILDER_VERSION}"
                + (f";read_length={language_read_length}" if language_read_length != 125 else ""),
            }
            additional_paths.append(LANGUAGE_PATHS[key])
    elif terminal_codebook_path is not None and internal_codebook_path is not None:
        sources = {"terminal_final": Path(terminal_codebook_path), "internal_only": Path(internal_codebook_path)}
        (root / "languages").mkdir(parents=True, exist_ok=True)
        for key, source in sources.items():
            if not source.is_file():
                raise InputContractError(f"{key} codebook does not exist: {source}")
            shutil.copyfile(source, root / LANGUAGE_PATHS[key])
        validate_language_files(root / LANGUAGE_PATHS["terminal_final"], root / LANGUAGE_PATHS["internal_only"])
        for key in LANGUAGE_KEYS:
            languages[key] = {
                "path": LANGUAGE_PATHS[key],
                "sha256": sha256_file(root / LANGUAGE_PATHS[key]),
                "source": "supplied",
            }
            additional_paths.append(LANGUAGE_PATHS[key])
    file_hashes = hash_relative_files(
        root,
        [te_relative, map_relative, *asset_relative_paths, *additional_paths],
    )
    base = {
        "schema_version": TE_MANIFEST_SCHEMA,
        "architecture_id": ARCHITECTURE_ID,
        "te_raw_sha256": te_library.raw_sha256,
        "te_semantic_sha256": te_library.semantic_sha256,
        "family_assignments": dict(validated_map),
        "reference_identity_sha256": prepared_reference.identity_sha256,
        "file_hashes": file_hashes,
        "resources": resources,
        "codebooks": codebooks,
        "languages": languages,
    }
    identity = manifest_sha256(base)
    write_manifest(manifest_path, dict(base, identity_sha256=identity))
    return load_prepared_te_assets(root)


def load_prepared_te_assets(path: str | Path) -> PreparedTEAssets:
    """Load and structurally verify a prepared-TE manifest."""

    candidate = Path(path)
    manifest_path = candidate if candidate.is_file() else candidate / TE_MANIFEST_NAME
    document = load_manifest(manifest_path)
    if document.get("schema_version") != TE_MANIFEST_SCHEMA:
        raise AssetCompatibilityError("unsupported prepared-TE schema")
    if document.get("architecture_id") != ARCHITECTURE_ID:
        raise AssetCompatibilityError("prepared TE architecture mismatch")
    identity = document.get("identity_sha256")
    base = dict(document)
    base.pop("identity_sha256", None)
    if not isinstance(identity, str) or manifest_sha256(base) != identity:
        raise AssetCompatibilityError("prepared-TE manifest identity mismatch")
    try:
        return PreparedTEAssets(
            root=manifest_path.parent,
            architecture_id=str(document["architecture_id"]),
            te_raw_sha256=str(document["te_raw_sha256"]),
            te_semantic_sha256=str(document["te_semantic_sha256"]),
            family_assignments=document["family_assignments"],
            reference_identity_sha256=str(document["reference_identity_sha256"]),
            file_hashes=document["file_hashes"],
            resources=document["resources"],
            codebooks=document["codebooks"],
            identity_sha256=identity,
            languages=document["languages"],
        )
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise AssetCompatibilityError("malformed prepared-TE manifest") from exc


def verify_te_assets(
    asset: PreparedTEAssets,
    prepared_reference: PreparedReference | None = None,
) -> None:
    """Verify TE assets and their optional prepared-reference binding."""

    if asset.architecture_id != ARCHITECTURE_ID:
        raise AssetCompatibilityError("prepared TE architecture mismatch")
    if manifest_sha256(asset.identity_document()) != asset.identity_sha256:
        raise AssetCompatibilityError("prepared-TE in-memory identity mismatch")
    verify_relative_files(asset.root, asset.file_hashes)
    if bool(asset.resources) != bool(asset.codebooks):
        raise AssetCompatibilityError("prepared-TE frozen resource/codebook manifest is incomplete")
    if asset.resources:
        if set(asset.resources) != {"g1a_archive", "family_crosswalk"}:
            raise AssetCompatibilityError("prepared-TE frozen resource manifest is malformed")
        if set(asset.codebooks) != {"terminal", "internal_full", "internal_only"}:
            raise AssetCompatibilityError("prepared-TE codebook manifest is malformed")
        for entries in (asset.resources, asset.codebooks):
            for item in entries.values():
                if set(item) != {"path", "sha256"}:
                    raise AssetCompatibilityError("prepared-TE resource entry is malformed")
                if asset.file_hashes.get(item["path"]) != item["sha256"]:
                    raise AssetCompatibilityError("prepared-TE resource/codebook hash binding mismatch")
    if asset.languages:
        if set(asset.languages) != set(LANGUAGE_KEYS):
            raise AssetCompatibilityError("prepared-TE language manifest is malformed")
        for item in asset.languages.values():
            if set(item) != {"path", "sha256", "source"}:
                raise AssetCompatibilityError("prepared-TE language entry is malformed")
            if asset.file_hashes.get(item["path"]) != item["sha256"]:
                raise AssetCompatibilityError("prepared-TE language hash binding mismatch")
    if prepared_reference is not None:
        verify_reference_asset(prepared_reference)
        if asset.reference_identity_sha256 != prepared_reference.identity_sha256:
            raise AssetCompatibilityError("TE/reference prepared-asset mismatch")


def build_safe_language(
    candidates: Iterable[SafeProbeCandidate],
    probe_class: ProbeClass,
) -> tuple[TEProbe, ...]:
    """Build the safe language through the Phase 4 semantic owner.

    Imports are local so Phase 3 asset loading stays independent of evidence
    construction and its candidate types.
    """

    from ..evidence.codebook import build_safe_probes

    return build_safe_probes(candidates, probe_class)


def validate_probe_independence(probes: Iterable[TEProbe]) -> None:
    from ..evidence.codebook import validate_probe_independence as validate

    validate(probes)


def validate_language_files(terminal_codebook: str | Path, internal_codebook: str | Path) -> None:
    """Parse both languages with the loaders the engines use (fail closed)."""

    from ..anchor_pe.codebook import load_safe_codebook
    from ..internal21.acquisition import load_internal_codebook

    try:
        load_safe_codebook(terminal_codebook)
        load_internal_codebook(internal_codebook)
    except (KeyError, ValueError) as exc:
        raise AssetCompatibilityError(f"malformed v3 language file: {exc}") from exc
