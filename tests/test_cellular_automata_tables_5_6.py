# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Verify Bhuvaneswari & Bhattacharjee 2026 Tables 5 and 6 period
claims against the C++ production kernel.

For each (k1, k2) combination in the paper, the kernel's
``is_full_period(s)`` verdict on each component must agree with the
gcd-coprimality predicate ``gcd(s, rho_i) = 1``: the time-spaced CA
T^s is primitive iff T is primitive AND gcd(s, rho_i) = 1.

The paper's ME claims are **not** re-checked here: the matricial-DE
kernel and the paper use different equidistribution conventions and
return divergent verdicts on many of these entries (see
``docs/generators/CellularAutomataGen.md`` "Two conventions for the
equidistribution matrix B"). The paper-vs-kernel ME comparison was the
investigation conducted by the (deleted) ``check_tables_5_6.py``
script; its finding is now documented and is not a stable regression
target.

Table 5 (4 entries) runs in full. Table 6 (~LARGE entries) is sampled
to the first ``TABLE_6_LIMIT`` entries (default 20); set
``TABLE_6_LIMIT=0`` to run the whole table.
"""

from __future__ import annotations

import importlib.resources
import importlib.util
import json
import math
import os
from pathlib import Path

import pytest

from regpoly.core.generator import Generator

_THIS_DIR = Path(__file__).parent
_ROOT = _THIS_DIR.parents[2]
_CA_DATA_PATH = importlib.resources.files("regpoly") / "data" / "cellular_automata.json"


def _load_paper_tables():
    spec = importlib.util.spec_from_file_location(
        "_paper_tables_5_6", _THIS_DIR / "data" / "paper_tables_5_6.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.TABLE_5, mod.TABLE_6


TABLE_5, TABLE_6 = _load_paper_tables()

with open(_CA_DATA_PATH) as f:
    _CA_DATA = json.load(f)
_CZ = _CA_DATA["cattell_zhang_1995"]


def _cz_positions(k):
    return _CZ[str(k)]


def _make_ca(k, s):
    return Generator.create(
        "CellularAutomataGen",
        L=min(k, 64), k=k, rule150_positions=_cz_positions(k), s=s,
    )


def _check_period(k1, k2, s):
    """Returns (period_kernel_ok, gcd_s_rho1, gcd_s_rho2)."""
    rho1 = (1 << k1) - 1
    rho2 = (1 << k2) - 1
    g1 = math.gcd(s, rho1)
    g2 = math.gcd(s, rho2)
    ca1 = _make_ca(k1, s)
    ca2 = _make_ca(k2, s)
    prim1 = ca1._cpp_gen.is_full_period()
    prim2 = ca2._cpp_gen.is_full_period()
    coprime = (g1 == 1) and (g2 == 1)
    return (prim1 and prim2) == coprime, g1, g2


def _verify_rows_period_only(rows):
    """Run period verification across rows; return (ok_count, total_count, failures)."""
    period_ok = period_total = 0
    failures = []
    for (k1, k2, s_period, _s_me, _rho_me) in rows:
        for s in s_period:
            ok, g1, g2 = _check_period(k1, k2, s)
            period_total += 1
            if ok:
                period_ok += 1
            else:
                failures.append((k1, k2, s, g1, g2))
    return period_ok, period_total, failures


@pytest.mark.slow
def test_paper_table_5_period_claim():
    """Every Table 5 row's period claim agrees with the C++ kernel."""
    ok, tot, failures = _verify_rows_period_only(TABLE_5)
    assert ok == tot, f"period: {ok}/{tot} agree; failures={failures[:5]}"


@pytest.mark.slow
def test_paper_table_6_period_claim():
    """Table 6 rows (sampled per TABLE_6_LIMIT, default 20) period claim."""
    limit_env = os.environ.get("TABLE_6_LIMIT", "20")
    try:
        limit_int = int(limit_env)
    except ValueError:
        limit_int = 20
    rows = TABLE_6 if limit_int <= 0 else TABLE_6[:limit_int]
    ok, tot, failures = _verify_rows_period_only(rows)
    assert ok == tot, (
        f"period: {ok}/{tot} agree; failures={failures[:5]}"
        f"{'...' if len(failures) > 5 else ''}"
    )
