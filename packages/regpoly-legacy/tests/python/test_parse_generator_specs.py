# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Parametrised tests asserting that the pure-Python parser produces the
frozen golden output captured from the C++ parser before its removal.

If this test ever fails, either the parser drifted or the golden file
needs to be regenerated against a corrected expectation.
"""

from __future__ import annotations

import json

import pytest

from pathlib import Path

from regpoly_legacy.reader import parse_generator_specs

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


def _gen_golden_pairs():
    if not GOLDEN_DIR.is_dir():
        return []
    out = []
    fdir = fixtures_dir()
    for g in sorted(GOLDEN_DIR.iterdir()):
        if not g.name.endswith(".gen.json"):
            continue
        stem = g.name[: -len(".gen.json")]
        fixture = fdir / stem
        if not fixture.is_file():
            # marsaxorshift_subcodes_21_25 lives next to the goldens, not in fixtures/
            alt = GOLDEN_DIR.parent / stem
            if alt.is_file():
                fixture = alt
            else:
                continue
        out.append((fixture, g))
    return out


@pytest.mark.parametrize(
    "fixture,golden",
    _gen_golden_pairs(),
    ids=lambda p: p.name if hasattr(p, "name") else str(p),
)
def test_matches_golden(fixture, golden):
    expected = json.loads(golden.read_text())
    actual = _norm(parse_generator_specs(str(fixture), 32))
    assert actual == expected
