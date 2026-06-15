# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Audit single-component cellular-automaton generators against the
claims in Bhuvaneswari & Bhattacharjee 2026:

- R1 (k=32, weak 16-rule150 config) — paper: NOT primitive, NOT ME.
- R3 (k=35, CA(150')) — paper: primitive, NOT ME.
- R4 (k=64, rule150 at [2,4]) — paper: NOT primitive, NOT ME.
- Table 2 / Cattell-Zhang 1995 single CAs (100 entries, k=29..128) —
  paper: all primitive.
- Adak-Das 2021 CA(90')/CA(150') — paper: primitive at listed k.
- Table 1 R1 per-t equi verdicts.

R2 (k=1409) is skipped here; ``test_cellular_automata_r2_slow.py``
covers it via the same C++ primitivity entry point with the precomputed
factor table.
"""

from __future__ import annotations

import importlib.resources
import json
from pathlib import Path

import pytest

from regpoly import make_combined
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.generator import Generator

_ROOT = Path(__file__).resolve().parents[3]
_CA_DATA_PATH = importlib.resources.files("regpoly") / "data" / "cellular_automata.json"


def _load_ca_data():
    with open(_CA_DATA_PATH) as f:
        return json.load(f)


_CA_DATA = _load_ca_data()
_CZ = _CA_DATA["cattell_zhang_1995"]
_ADAK_90P = _CA_DATA["adak_das_2021_CA90prime_k_values"]
_ADAK_150P = _CA_DATA["adak_das_2021_CA150prime_k_values"]


def _make_single(k, positions, s=1):
    return Generator.create(
        "CellularAutomataGen", L=min(k, 64),
        k=k, rule150_positions=positions, s=s,
    )


def _is_me_single(k, positions, s=1):
    g = _make_single(k, positions, s)
    comb = make_combined(g, Lmax=64)
    test = EquidistributionTest(
        L=min(k, 64),
        delta=[10**9] * (min(k, 64) + 1),
        mse=10**9, method=None,
    )
    res = test.run(comb)
    return res.is_me(), tuple(res.ecart)


# Each row: (label, k, positions, kernel_primitive, kernel_me).
# The kernel verdicts are what regpoly_cpp returns today; they are the
# regression target. The paper's stated claims sometimes disagree with
# the kernel — see docs/generators/CellularAutomataGen.md "Two
# conventions for the equidistribution matrix B" for the convention
# investigation. Disagreements are flagged below in comments.
_R_QUARTET = [
    pytest.param(
        # Paper: NOT primitive, NOT ME. Kernel: primitive, NOT ME
        # (paper-vs-kernel disagreement on primitivity).
        "R1", 32,
        [1, 5, 6, 12, 15, 16, 18, 19, 20, 22, 23, 24, 25, 27, 29, 31],
        True, False, id="R1_k32_weak16",
    ),
    pytest.param(
        # Paper: primitive, NOT ME. Kernel agrees on both.
        "R3", 35, list(range(1, 35)), True, False, id="R3_k35_CA150prime",
    ),
    pytest.param(
        # Paper: NOT primitive, NOT ME. Kernel: primitive, NOT ME
        # (paper-vs-kernel disagreement on primitivity).
        "R4", 64, [2, 4], True, False, id="R4_k64_rule150_2_4",
    ),
]


@pytest.mark.parametrize(
    ("label", "k", "positions", "kernel_primitive", "kernel_me"),
    _R_QUARTET,
)
def test_r_quartet_kernel_verdicts(label, k, positions, kernel_primitive, kernel_me):
    """R1, R3, R4 primitive/ME verdicts lock in the C++ kernel output."""
    g = _make_single(k, positions)
    prim = g._cpp_gen.is_full_period()
    assert prim is kernel_primitive, (
        f"{label} primitive: kernel returned {prim}, expected {kernel_primitive}"
    )
    me, _ = _is_me_single(k, positions)
    assert me is kernel_me, (
        f"{label} ME: kernel returned {me}, expected {kernel_me}"
    )


@pytest.mark.slow
def test_cattell_zhang_1995_all_primitive():
    """Every Cattell-Zhang 1995 single CA (100 entries, k=29..128) is primitive."""
    failures = []
    for k_str, positions in sorted(_CZ.items(), key=lambda kv: int(kv[0])):
        k = int(k_str)
        g = _make_single(k, positions)
        if not g._cpp_gen.is_full_period():
            failures.append((k, positions))
    assert not failures, (
        f"{len(failures)} Cattell-Zhang entries not primitive: "
        f"{failures[:5]}{'...' if len(failures) > 5 else ''}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("k", sorted(_ADAK_90P))
def test_adak_das_2021_ca90prime_primitive(k):
    """Adak-Das 2021 CA(90') at every listed k is primitive."""
    g = _make_single(k, [0])
    assert g._cpp_gen.is_full_period(), f"CA(90') k={k} not primitive"


@pytest.mark.slow
@pytest.mark.parametrize("k", sorted(_ADAK_150P))
def test_adak_das_2021_ca150prime_primitive(k):
    """Adak-Das 2021 CA(150') at every listed k is primitive."""
    positions = list(range(1, k))
    g = _make_single(k, positions)
    assert g._cpp_gen.is_full_period(), f"CA(150') k={k} not primitive"


# Paper Table 1, R1 per-t verdicts. The kernel and paper disagree on
# numerical rank values; we compare the equi verdict only (see the docs:
# docs/generators/CellularAutomataGen.md "Two conventions for B").
_TABLE_1_ROWS = [
    # (t, ell_star_t, paper_says_equi)
    (2, 16, False),
    (3, 10, False),
    (4,  8, False),
    (5,  6, False),
    (6,  5, False),
    (8,  4, False),
    (10, 3, False),
    (16, 2, False),
    (32, 1, True),
]


def test_table1_r1_per_t_equi_verdicts():
    """Paper Table 1 R1: kernel ME verdict per t agrees with paper."""
    positions = [1, 5, 6, 12, 15, 16, 18, 19, 20, 22, 23, 24, 25, 27, 29, 31]
    g = _make_single(32, positions)
    comb = make_combined(g, Lmax=64)
    test = EquidistributionTest(L=32, delta=[10**9] * 33,
                                mse=10**9, method=None)
    res = test.run(comb)
    lam = res._conv_ecarts(comb)
    for t, _ell_star, paper_equi in _TABLE_1_ROWS:
        lam_t = lam[t] if t < len(lam) else 0
        kernel_equi = (lam_t <= 0)
        assert kernel_equi == paper_equi, (
            f"t={t}: kernel ME={kernel_equi}, paper={paper_equi} "
            f"(Lambda_t={lam_t})"
        )
