"""Workspace-local test fixtures for managed Windows environments."""

from __future__ import annotations

import os
import re
import shutil
import stat
from pathlib import Path
from typing import Any, Callable

import pytest

_TEMP_ROOT = Path.cwd() / ".evitea-test-tmp"


def _remove_tree(path: Path) -> None:
    if not path.exists():
        return

    def make_writable_and_retry(
        function: Callable[..., Any], target: str, _: BaseException
    ) -> None:
        os.chmod(target, stat.S_IWRITE)
        function(target)

    shutil.rmtree(path, onexc=make_writable_and_retry)


@pytest.fixture
def tmp_path(request: pytest.FixtureRequest) -> Path:
    """Provide a writable test directory without pytest's mode-0700 temp root."""

    _TEMP_ROOT.mkdir(mode=0o777, exist_ok=True)
    safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", request.node.nodeid)
    case_path = _TEMP_ROOT / safe_name
    _remove_tree(case_path)
    case_path.mkdir(mode=0o777)
    try:
        yield case_path
    finally:
        _remove_tree(case_path)


def pytest_sessionfinish() -> None:
    """Remove the empty workspace-local test root after every run."""

    _remove_tree(_TEMP_ROOT)
