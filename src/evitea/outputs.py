"""Atomic writers for run-level outputs (manifest, summary, log, raw bytes).

Call-level views are rendered by ``evitea.callset`` from the native callset.
"""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from .contracts import ProvenanceRecord, stable_json_dumps, stable_to_primitive
from .errors import InputContractError


def _atomic_write(path: str | Path, data: bytes) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=target.parent, prefix=f".{target.name}.", delete=False
        ) as handle:
            temporary_name = handle.name
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, target)
    except OSError:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
        raise
    return target


def write_bytes(path: str | Path, data: bytes) -> Path:
    """Write exact bytes atomically (no re-encoding, no newline translation)."""

    if not isinstance(data, bytes):
        raise InputContractError("write_bytes requires bytes")
    return _atomic_write(path, data)


def _json_bytes(value: Any) -> bytes:
    return (stable_json_dumps(stable_to_primitive(value)) + "\n").encode("utf-8")


def write_run_manifest(path: str | Path, manifest: ProvenanceRecord) -> Path:
    return _atomic_write(path, _json_bytes(manifest))


def write_run_summary(path: str | Path, summary: Mapping[str, Any]) -> Path:
    forbidden = {key for key in summary if "score" in key.lower() or "truth" in key.lower()}
    if forbidden:
        raise InputContractError(
            f"public run summary cannot contain benchmark fields: {sorted(forbidden)}"
        )
    return _atomic_write(path, _json_bytes(summary))


def write_run_log(path: str | Path, messages: Iterable[str]) -> Path:
    lines = tuple(messages)
    if any("\n" in line or "\r" in line for line in lines):
        raise InputContractError("run log messages must be single-line text")
    return _atomic_write(path, (("\n".join(lines) + "\n") if lines else "").encode("utf-8"))
