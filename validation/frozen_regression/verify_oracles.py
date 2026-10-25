"""Check the recorded callset files, in order, against their expected row counts and SHA-256 hashes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
import sys
from pathlib import Path, PureWindowsPath
from typing import Any, Sequence


ORACLES_PATH = Path(__file__).with_name("oracles.json")


def _load_oracles() -> list[dict[str, Any]]:
    with ORACLES_PATH.open("r", encoding="utf-8") as handle:
        document = json.load(handle)
    if document.get("schema_version") != 1 or not isinstance(
        document.get("checks"), list
    ):
        raise ValueError("unsupported oracles.json schema")
    return document["checks"]


def _relative_oracle_path(value: str) -> Path:
    """Reject absolute paths and ``..`` so checks stay within the given root."""
    path = Path(value)
    windows_path = PureWindowsPath(value)
    if (
        path.is_absolute()
        or windows_path.is_absolute()
        or any(part == ".." for part in path.parts)
        or any(part == ".." for part in windows_path.parts)
    ):
        raise ValueError("recorded file path must be relative to its route root")
    return path


def _inspect_file(path: Path) -> tuple[int | None, str | None, str | None]:
    """Return the data-row count and SHA-256 digest, or the reason the file could not be read."""
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        return None, None, "missing_file"
    except OSError:
        return None, None, "file_stat_error"

    if not stat.S_ISREG(metadata.st_mode):
        return None, None, "not_regular_file"

    digest = hashlib.sha256()
    line_count = 0
    try:
        with path.open("rb") as handle:
            if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
                return None, None, "not_regular_file"
            while line := handle.readline():
                digest.update(line)
                line_count += 1
    except OSError:
        return None, None, "file_read_error"

    # The first line is the header, so an empty file has zero data rows and
    # fails the row-count check.
    return max(line_count - 1, 0), digest.hexdigest(), None


def verify_oracles(
    route_b_v2: Path | str,
    route_b_v3: Path | str,
    *,
    checks: Sequence[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Check each recorded file in manifest order and stop at the first mismatch."""
    roots = {
        "route_b_v2": Path(route_b_v2),
        "route_b_v3": Path(route_b_v3),
    }
    oracle_checks = list(_load_oracles() if checks is None else checks)
    report: dict[str, Any] = {
        "schema_version": 1,
        "scope": "FROZEN_ORACLE_INTEGRITY_ONLY",
        "method_equivalence": "NOT_EVALUATED",
        "status": "PASS",
        "checks": [],
    }

    for index, check in enumerate(oracle_checks, start=1):
        route_root = check["route_root"]
        relative_file = _relative_oracle_path(check["file"])
        expected_rows = check["data_rows"]
        expected_sha256 = check["sha256"]
        evidence: dict[str, Any] = {
            "index": index,
            "layer": check["layer"],
            "dataset": check["dataset"],
            "route_root": route_root,
            "file": relative_file.as_posix(),
            "expected_data_rows": expected_rows,
            "actual_data_rows": None,
            "expected_sha256": expected_sha256,
            "actual_sha256": None,
            "status": "PASS",
        }

        if route_root not in roots:
            reason = "unknown_route_root"
            actual_rows, actual_sha256 = None, None
        else:
            actual_rows, actual_sha256, reason = _inspect_file(
                roots[route_root] / relative_file
            )
        evidence["actual_data_rows"] = actual_rows
        evidence["actual_sha256"] = actual_sha256

        if reason is None and actual_rows != expected_rows:
            reason = "data_row_count_mismatch"
        if reason is None and actual_sha256 != expected_sha256:
            reason = "sha256_mismatch"
        if reason is not None:
            evidence["status"] = "FAIL"
            evidence["reason"] = reason
            report["status"] = "FAIL"
            report["checks"].append(evidence)
            break

        report["checks"].append(evidence)

    return report


def _render_report(report: dict[str, Any]) -> str:
    return json.dumps(report, indent=2, sort_keys=True) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route-b-v2", required=True, type=Path)
    parser.add_argument("--route-b-v3", required=True, type=Path)
    parser.add_argument("--output", type=Path, help="also save the JSON report here")
    args = parser.parse_args(argv)

    try:
        report = verify_oracles(args.route_b_v2, args.route_b_v3)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        parser.error(f"could not load the manifest of recorded files: {exc}")
    rendered = _render_report(report)

    if args.output is not None:
        try:
            args.output.write_text(rendered, encoding="utf-8", newline="\n")
        except OSError as exc:
            print(f"could not write report: {exc}", file=sys.stderr)
            return 2
    sys.stdout.write(rendered)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
