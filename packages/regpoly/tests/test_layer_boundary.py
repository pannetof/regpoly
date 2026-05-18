# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Architectural guard: tests in this directory must not depend on the
web layer.

The workspace has a strict one-way dependency arrow
`regpoly_web -> regpoly -> regpoly_cpp` (enforced for `src/` code by
import-linter contracts in the workspace `pyproject.toml`). Test files
live outside the import-linter root packages, so this test plugs the
gap on the `regpoly` side: a test that imports `regpoly_web.*` should
live under `packages/regpoly-web/tests/`, not here.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_SELF = Path(__file__).resolve()

_FORBIDDEN_ROOTS = ("regpoly_web",)


def _sibling_python_files():
    return sorted(p for p in _HERE.glob("*.py") if p.resolve() != _SELF)


def _imported_top_levels(tree: ast.AST):
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            yield node.module.split(".", 1)[0], node
        elif isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name.split(".", 1)[0], node


@pytest.mark.parametrize("path", _sibling_python_files(), ids=lambda p: p.name)
def test_file_does_not_import_regpoly_web(path: Path) -> None:
    tree = ast.parse(path.read_text(), filename=str(path))
    bad = [
        (top, node.lineno)
        for top, node in _imported_top_levels(tree)
        if top in _FORBIDDEN_ROOTS
    ]
    assert not bad, (
        f"{path.name} imports forbidden top-level module(s) {bad!r}. "
        f"This directory is reserved for tests of the `regpoly` Python "
        f"wrapper layer. Move this file to packages/regpoly-web/tests/ "
        f"to respect the one-way dependency arrow "
        f"regpoly_web -> regpoly -> regpoly_cpp."
    )
