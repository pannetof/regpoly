"""Verify paper Tables 5 and 6 (Bhuvaneswari & Bhattacharjee 2026):

For each (k1, k2) combination:
  1. Period claim: for s in s_period, paper says gcd(s, ρ_i) = 1 for each
     component i (so the combined PRNG reaches close-to-maximal period
     LCM(ρ_1, ρ_2)).  We re-check via gcd directly and via
     CellularAutomataGen.is_full_period(s=s).
  2. ME claim: for s in s_ME, paper says the combined PRNG is maximally
     equidistributed.  We re-test via the C++ production kernel
     (EquidistributionTest with the matricial method) and compare
     verdicts.

Component rule150 positions are taken from
packages/regpoly/src/regpoly/data/cellular_automata.json
(Cattell-Zhang 1995, 0-indexed).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "packages/regpoly/tests"))

from data.paper_tables_5_6 import TABLE_5, TABLE_6   # noqa: E402

from regpoly.analyses.equidistribution_test import (
    EquidistributionTest, METHOD_MATRICIAL, METHOD_NOTPRIMITIVE,
)
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator

# Override via env var: METHOD=matricial (default) or METHOD=notprimitive
import os as _os
_METHOD = {
    "matricial":    METHOD_MATRICIAL,
    "notprimitive": METHOD_NOTPRIMITIVE,
    "default":      None,
}[_os.environ.get("ME_METHOD", "default")]


with open(ROOT / "packages/regpoly/src/regpoly/data/cellular_automata.json") as f:
    CA_DATA = json.load(f)
CZ = CA_DATA["cattell_zhang_1995"]


def cz_positions(k):
    return CZ[str(k)]


def check_period(k1, k2, s):
    """Returns (period_kernel_ok, gcd_s_rho1, gcd_s_rho2).
    The paper's selection criterion: gcd(s, ρ_i) = 1 for each component
    (so the time-spaced CA T^s is still primitive, i.e. reaches the full
    component period 2^k - 1)."""
    rho1 = (1 << k1) - 1
    rho2 = (1 << k2) - 1
    g1 = gcd_int = math.gcd(s, rho1)
    g2 = math.gcd(s, rho2)

    # Kernel cross-check: each component as CA(s) primitive iff
    # underlying T is primitive AND gcd(s, ρ_i) = 1.
    p1 = cz_positions(k1)
    p2 = cz_positions(k2)
    ca1 = Generator.create("CellularAutomataGen", L=min(k1, 64),
                           k=k1, rule150_positions=p1, s=s)
    ca2 = Generator.create("CellularAutomataGen", L=min(k2, 64),
                           k=k2, rule150_positions=p2, s=s)
    prim1 = ca1._cpp_gen.is_full_period()
    prim2 = ca2._cpp_gen.is_full_period()
    coprime = (g1 == 1) and (g2 == 1)
    return (prim1 and prim2) == coprime, g1, g2


def check_me(k1, k2, s):
    """Returns the kernel's ME verdict (bool) for the combined PRNG at s."""
    p1 = cz_positions(k1)
    p2 = cz_positions(k2)
    g1 = Generator.create("CellularAutomataGen", L=min(k1, 64),
                          k=k1, rule150_positions=p1, s=s)
    g2 = Generator.create("CellularAutomataGen", L=min(k2, 64),
                          k=k2, rule150_positions=p2, s=s)
    comb = Combination.CreateFromFiles([[g1], [g2]], Lmax=64,
                                       temperings=[[], []])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=_METHOD)
    res = test.run(comb)
    return res.is_me()


def run_table(name, rows, limit=None):
    print(f"\n{'=' * 78}")
    print(f"Paper {name} — {len(rows)} (k1,k2) combinations"
          + (f"  [sampling first {limit}]" if limit else ""))
    print(f"{'=' * 78}")
    period_ok = period_total = 0
    me_ok = me_total = me_paper_me = me_paper_notme = 0

    rows_to_check = rows[:limit] if limit else rows

    for (k1, k2, s_period, s_me, rho_me) in rows_to_check:
        # Period: for each s in s_period, the kernel-vs-coprime expectation
        # should agree.
        for s in s_period:
            ok, g1, g2 = check_period(k1, k2, s)
            period_total += 1
            period_ok += int(ok)
            if not ok:
                print(f"  PERIOD MISMATCH (k1={k1}, k2={k2}, s={s}): "
                      f"gcd(s,ρ1)={g1} gcd(s,ρ2)={g2}")

        # ME: paper claims ME for s in s_me.
        for s in range(2, 11):
            paper_says_me = (s in s_me)
            kernel_says_me = check_me(k1, k2, s)
            me_total += 1
            if paper_says_me:
                me_paper_me += 1
                if kernel_says_me:
                    me_ok += 1
            else:
                me_paper_notme += 1
                if not kernel_says_me:
                    me_ok += 1

    print(f"\n  Period verification: {period_ok}/{period_total} agree")
    print(f"  ME verification:     {me_ok}/{me_total} agree")
    print(f"    (paper ME=true:    {me_paper_me};  paper ME=false: {me_paper_notme})")
    return period_ok, period_total, me_ok, me_total


if __name__ == "__main__":
    # Table 5 (4 entries) — full
    p_ok_5, p_tot_5, m_ok_5, m_tot_5 = run_table("Table 5", TABLE_5)

    # Table 6 — sample (first 20 entries, ~120 ME calls).  Set TABLE_6_LIMIT=0
    # in env to run all (slow).
    import os
    limit = int(os.environ.get("TABLE_6_LIMIT", "20") or 0) or None
    p_ok_6, p_tot_6, m_ok_6, m_tot_6 = run_table("Table 6", TABLE_6, limit=limit)

    print(f"\n{'=' * 78}")
    print("AGGREGATE")
    print(f"{'=' * 78}")
    print(f"  Period: {p_ok_5 + p_ok_6}/{p_tot_5 + p_tot_6} agree")
    print(f"  ME:     {m_ok_5 + m_ok_6}/{m_tot_5 + m_tot_6} agree")
