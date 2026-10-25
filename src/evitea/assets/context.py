"""Deterministic hashing and manifest helpers for prepared assets."""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from ..contracts import stable_to_primitive
from ..errors import AssetCompatibilityError, InputContractError


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 digest of the exact file bytes."""

    file_path = Path(path)
    if not file_path.is_file():
        raise InputContractError(f"file does not exist: {file_path}")
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    """Return canonical UTF-8 JSON bytes with one trailing newline."""

    return (
        json.dumps(
            stable_to_primitive(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def manifest_sha256(value: Any) -> str:
    """Hash a manifest's canonical semantic bytes."""

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def write_manifest(path: str | Path, value: Any) -> Path:
    """Atomically write a deterministic manifest."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp-{os.getpid()}")
    try:
        with temporary.open("wb") as handle:
            handle.write(canonical_json_bytes(value))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    return output


def load_manifest(path: str | Path) -> dict[str, Any]:
    """Load a manifest file that holds one JSON object."""

    manifest_path = Path(path)
    try:
        with manifest_path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise AssetCompatibilityError(
            f"cannot load manifest {manifest_path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise AssetCompatibilityError("manifest root must be a JSON object")
    return value


def hash_relative_files(
    root: str | Path,
    relative_paths: Iterable[str | Path],
) -> Mapping[str, str]:
    """Hash explicitly named files under a prepared-asset root."""

    root_path = Path(root).resolve()
    result: dict[str, str] = {}
    for relative in relative_paths:
        relative_path = Path(relative)
        if relative_path.is_absolute() or ".." in relative_path.parts:
            raise InputContractError(f"asset path must be root-relative: {relative}")
        absolute = (root_path / relative_path).resolve()
        try:
            absolute.relative_to(root_path)
        except ValueError as exc:
            raise InputContractError(f"asset path escapes root: {relative}") from exc
        key = relative_path.as_posix()
        if key in result:
            raise InputContractError(f"duplicate asset path: {key}")
        result[key] = sha256_file(absolute)
    return dict(sorted(result.items()))


def verify_relative_files(
    root: str | Path,
    expected: Mapping[str, str],
) -> None:
    """Raise unless every bound file still matches its recorded hash."""

    actual = hash_relative_files(root, expected)
    if actual != dict(expected):
        mismatched = sorted(
            path for path in set(actual) | set(expected) if actual.get(path) != expected.get(path)
        )
        raise AssetCompatibilityError(f"prepared asset hash mismatch: {mismatched}")
