# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

from __future__ import annotations

import json

import pytest

from pathlib import Path

from regpoly_legacy.reader import parse_transformation_specs

_THIS = Path(__file__).resolve()
GOLDEN_DIR = _THIS.parent / "data" / "golden"


def fixtures_dir() -> Path:
    return _THIS.parents[2] / "shared" / "legacy_parameters"


def _norm(o):
    if isinstance(o, tuple):
        return [_norm(x) for x in o]
    if isinstance(o, list):
        return [_norm(x) for x in o]
    if isinstance(o, dict):
        return {k: _norm(v) for k, v in o.items()}
    return o


def _trans_pairs():
    if not GOLDEN_DIR.is_dir():
        return []
    out = []
    fdir = fixtures_dir()
    for g in sorted(GOLDEN_DIR.iterdir()):
        if not g.name.endswith(".trans.json"):
            continue
        stem = g.name[: -len(".trans.json")]
        fixture = fdir / stem
        if not fixture.is_file():
            continue
        out.append((fixture, g))
    return out


@pytest.mark.parametrize(
    "fixture,golden",
    _trans_pairs(),
    ids=lambda p: p.name if hasattr(p, "name") else str(p),
)
def test_matches_golden(fixture, golden):
    expected = json.loads(golden.read_text())
    actual = _norm(parse_transformation_specs(str(fixture)))
    assert actual == expected
