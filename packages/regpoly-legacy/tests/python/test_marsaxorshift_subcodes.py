# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Coverage gain: exercise marsaxorshift sub-codes 21..25 which the
deleted C++ test suite did not cover. Fixture is hand-authored and
captured against the C++ parser before its removal.
"""

from __future__ import annotations

import json

from pathlib import Path

from regpoly_legacy.reader import parse_generator_specs

_THIS = Path(__file__).resolve()
EXTRA_DATA_DIR = _THIS.parent / "data"
GOLDEN_DIR = _THIS.parent / "data" / "golden"


def _norm(o):
    if isinstance(o, tuple):
        return [_norm(x) for x in o]
    if isinstance(o, list):
        return [_norm(x) for x in o]
    if isinstance(o, dict):
        return {k: _norm(v) for k, v in o.items()}
    return o


def test_subcodes_21_25_expand_correctly():
    fixture = EXTRA_DATA_DIR / "marsaxorshift_subcodes_21_25.dat"
    golden = GOLDEN_DIR / "marsaxorshift_subcodes_21_25.dat.gen.json"
    expected = json.loads(golden.read_text())
    actual = _norm(parse_generator_specs(str(fixture), 32))
    assert actual == expected


def test_subcodes_21_25_variant_counts():
    """Each sub-code expands to a specific number of variants:
    21→3, 22→3, 23→2, 24→2, 25→4. Total 14."""
    fixture = EXTRA_DATA_DIR / "marsaxorshift_subcodes_21_25.dat"
    specs = parse_generator_specs(str(fixture), 32)
    assert len(specs) == 3 + 3 + 2 + 2 + 4  # = 14
    for family, _params in specs:
        assert family == "MarsaXorshiftGen"
