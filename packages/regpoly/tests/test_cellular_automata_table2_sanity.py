# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Sanity check the Cattell-Zhang Table 2 transcription.

The indexing convention is resolved by Bhuvaneswari & Bhattacharjee 2026
Example 1 (p.16): Table 2 positions are 1-indexed; subtract 1 for
0-indexed positions stored in cellular_automata.json.

This test asserts that the C++ CA built with the transcribed positions
has a primitive characteristic polynomial — i.e. period 2^k - 1.
Failure means a transcription bug, not a code bug.

Covers every k in cellular_automata.json["cattell_zhang_1995"]
(k = 29..128) plus the Adak-Das CA(90')/CA(150') families.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from regpoly.core.generator import Generator


def _load_table():
    p = (Path(__file__).resolve().parents[3]
         / "packages" / "regpoly" / "src" / "regpoly" / "data" / "cellular_automata.json")
    with open(p) as f:
        return json.load(f)


_DATA = _load_table()
_CATTELL = _DATA["cattell_zhang_1995"]
_AD90 = _DATA["adak_das_2021_CA90prime_k_values"]
_AD150 = _DATA["adak_das_2021_CA150prime_k_values"]


def _ca_L(k: int) -> int:
    """Default L for a single CellularAutomataGen: min(k, 64)."""
    return min(k, 64)


@pytest.mark.parametrize(
    "k,positions",
    sorted([(int(k), p) for k, p in _CATTELL.items()]),
    ids=lambda v: f"k={v}" if isinstance(v, int) else str(v),
)
def test_cattell_zhang_primitive(k, positions):
    gen = Generator.create(
        "CellularAutomataGen", L=_ca_L(k),
        k=k, rule150_positions=positions, s=1,
    )
    assert gen._cpp_gen.is_full_period(), (
        f"k={k} positions={positions} not primitive"
    )


@pytest.mark.parametrize("k", _AD90)
def test_adak_das_CA90prime_primitive(k):
    """CA(90') at size k: rule-150 only at cell 0, rule-90 elsewhere."""
    gen = Generator.create(
        "CellularAutomataGen", L=_ca_L(k),
        k=k, rule150_positions=[0], s=1,
    )
    assert gen._cpp_gen.is_full_period()


@pytest.mark.parametrize("k", _AD150)
def test_adak_das_CA150prime_primitive(k):
    """CA(150') at size k: rule-150 everywhere except cell 0."""
    gen = Generator.create(
        "CellularAutomataGen", L=_ca_L(k),
        k=k, rule150_positions=list(range(1, k)), s=1,
    )
    assert gen._cpp_gen.is_full_period()
