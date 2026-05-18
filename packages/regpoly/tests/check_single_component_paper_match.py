"""Audit every single-component cellular-automaton generator the paper
makes a claim about, and tabulate kernel agreement with the paper.

The paper's single-component claims:

  1. Section 3.1, R1: k=32 single CA with rule-vector
     <90,150,90,90,90,150,150,90,90,90,90,90,150,90,90,150,
      150,90,150,150,150,90,150,150,150,150,90,150,90,150,90,150>
     → 0-indexed rule-150 positions [1,5,6,12,15,16,18,19,20,22,23,24,25,
                                       27,29,31].
     Paper claim: NOT primitive, NOT ME (Table 1 gives per-t numbers).

  2. Section 3.1, R2: k=1409 CA(150′) (rule 150 at every cell except 0).
     Paper claim: primitive, NOT ME (no equi table; too costly).

  3. Section 3.1, R3: k=35 CA(150′).  Paper claim: primitive, NOT ME.

  4. Section 3.1, R4: k=64 single CA, rule 150 only at cells {3,5} (1-idx)
     = 0-idx positions [2,4].  Paper claim: NOT primitive, NOT ME.

  5. Table 2 / Cattell-Zhang 1995, k=29..128: 100 entries, each rule-150
     position list certified primitive by Cattell-Zhang.  Paper restates
     this in Table 2.

  6. Adak-Das 2021 CA(90′)/CA(150′) for k ∈ {26, 29, 35, 39, 65, 69, 105,
     113, 119}: each claimed primitive (paper's Section 6.2 / Table 4).

This script runs is_full_period(.) on every single component and the
matricial-DE ME test on R1/R3/R4 (we skip the multi-minute R2).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


with open(ROOT / "packages/regpoly/src/regpoly/data/cellular_automata.json") as f:
    CA_DATA = json.load(f)
CZ        = CA_DATA["cattell_zhang_1995"]
ADAK_90P  = CA_DATA["adak_das_2021_CA90prime_k_values"]
ADAK_150P = CA_DATA["adak_das_2021_CA150prime_k_values"]


def make_single(k, positions, s=1):
    return Generator.create("CellularAutomataGen", L=min(k, 64),
                            k=k, rule150_positions=positions, s=s)


def is_me_single(k, positions, s=1):
    """Run the matricial equidistribution test on a single-component CA."""
    g = make_single(k, positions, s)
    comb = Combination.CreateFromFiles([[g]], Lmax=64, temperings=[[]])
    next(iter(comb))
    test = EquidistributionTest(L=min(k, 64),
                                delta=[10**9] * (min(k, 64) + 1),
                                mse=10**9, method=None)
    res = test.run(comb)
    return res.is_me(), tuple(res.ecart)


def heading(s):
    print(f"\n{'=' * 78}")
    print(s)
    print(f"{'=' * 78}")


def check_r_quartet():
    """R1, R3, R4 — paper claims NOT ME (with R1 explicitly NOT primitive).
    R2 (k=1409) is skipped here; covered by test_cellular_automata_r2_slow.py."""
    heading("Section 3.1 weak generators: R1, R3, R4")
    cases = [
        ("R1 (k=32, 16 rule150 cells)", 32,
            [1, 5, 6, 12, 15, 16, 18, 19, 20, 22, 23, 24, 25, 27, 29, 31],
            "paper: NOT primitive, NOT ME"),
        ("R3 (k=35, CA(150′))", 35, list(range(1, 35)),
            "paper: primitive, NOT ME"),
        ("R4 (k=64, rule150 at [2,4])", 64, [2, 4],
            "paper: NOT primitive, NOT ME"),
    ]
    rows = []
    for name, k, pos, claim in cases:
        g = make_single(k, pos)
        prim = g._cpp_gen.is_full_period()
        me, ecart = is_me_single(k, pos)
        ecart_sum = sum(e for e in ecart if e > 0)
        rows.append((name, claim, prim, me, ecart_sum))
        print(f"  {name:32s} {claim}")
        print(f"    kernel: primitive={prim}, ME={me}, ecart_sum={ecart_sum}")
    return rows


def check_cattell_zhang():
    heading("Table 2 — Cattell-Zhang 1995 single CAs (paper: all primitive)")
    total = len(CZ)
    matches = 0
    failures = []
    for k_str, positions in sorted(CZ.items(), key=lambda kv: int(kv[0])):
        k = int(k_str)
        g = make_single(k, positions)
        prim = g._cpp_gen.is_full_period()
        if prim:
            matches += 1
        else:
            failures.append((k, positions))
    print(f"  Primitivity: {matches}/{total} match paper claim")
    if failures:
        print(f"  Failures (paper says primitive, kernel says not):")
        for k, p in failures:
            print(f"    k={k:4d}, positions={p}")
    return matches, total


def check_adak_das():
    heading("Adak-Das 2021 CA(90′) and CA(150′) — paper: primitive at listed k")
    matches_90 = matches_150 = 0
    total_90 = len(ADAK_90P)
    total_150 = len(ADAK_150P)
    for k in ADAK_90P:
        g = make_single(k, [0])           # CA(90′): rule150 only at cell 0
        if g._cpp_gen.is_full_period():
            matches_90 += 1
        else:
            print(f"  CA(90′) k={k}: kernel says NOT primitive (paper says primitive)")
    for k in ADAK_150P:
        positions = list(range(1, k))     # CA(150′): rule150 at every cell except 0
        g = make_single(k, positions)
        if g._cpp_gen.is_full_period():
            matches_150 += 1
        else:
            print(f"  CA(150′) k={k}: kernel says NOT primitive (paper says primitive)")
    print(f"  CA(90′)  primitivity: {matches_90}/{total_90}")
    print(f"  CA(150′) primitivity: {matches_150}/{total_150}")
    return matches_90, total_90, matches_150, total_150


def check_table1_r1_per_t():
    """Reproduce paper's Table 1 per-t breakdown for R1.

    Paper Table 1 lists (t, l*_t, l_t, rank, verdict) for R1.  The verdicts
    are all 'NOT' except t=32 which is 'equi' (which is just t=L → trivial)."""
    heading("Table 1 — R1 per-t equidistribution numbers")
    print("  (Paper rank values are NOT what our matricial kernel returns; the")
    print("   kernel returns Λ_t = ℓ*_t - ℓ_t, not a B-matrix rank.  We compare")
    print("   the verdict columns only.)")
    paper_rows = [
        # (t, l_star_t, l_t_paper, rank_paper, equi_paper)
        (2, 16, 1, 18, False),
        (3, 10, 1, 13, False),
        (4,  8, 1, 12, False),
        (5,  6, 1, 11, False),
        (6,  5, 1, 11, False),
        (8,  4, 1, 12, False),
        (10, 3, 1, 13, False),
        (16, 2, 1, 17, False),
        (32, 1, 1, 32, True),
    ]
    pos = [1, 5, 6, 12, 15, 16, 18, 19, 20, 22, 23, 24, 25, 27, 29, 31]
    g = make_single(32, pos)
    comb = Combination.CreateFromFiles([[g]], Lmax=64, temperings=[[]])
    next(iter(comb))
    test = EquidistributionTest(L=32, delta=[10**9] * 33, mse=10**9, method=None)
    res = test.run(comb)
    # res.ecart is Λ_l per resolution; convert to Λ_t for paper-compatible t.
    lam = res._conv_ecarts(comb)
    print(f"\n  {'t':>3s}  {'ℓ*_t':>5s}  {'kernel ℓ_t':>10s}  {'paper ℓ_t':>9s}  "
          f"{'kernel equi':>12s}  {'paper equi':>11s}  match?")
    matches = 0
    for (t, l_star, l_paper, rank_paper, equi_paper) in paper_rows:
        lam_t = lam[t] if t < len(lam) else 0
        l_kernel = l_star - max(lam_t, 0)
        equi_kernel = (lam_t <= 0)
        ok = (equi_kernel == equi_paper)
        matches += int(ok)
        print(f"  {t:3d}  {l_star:5d}  {l_kernel:10d}  {l_paper:9d}  "
              f"{str(equi_kernel):>12s}  {str(equi_paper):>11s}  {'OK' if ok else 'mismatch'}")
    print(f"\n  Verdict match: {matches}/{len(paper_rows)}")
    return matches, len(paper_rows)


def main():
    r_rows = check_r_quartet()
    cz_matches, cz_total = check_cattell_zhang()
    ad_90, t_90, ad_150, t_150 = check_adak_das()
    t1_match, t1_total = check_table1_r1_per_t()

    heading("AGGREGATE — single-component agreement with paper")
    # R1, R3, R4 paper-vs-kernel: paper claims NOT ME, kernel should agree.
    r_me_match = sum(1 for (_, _, _, me, _) in r_rows if me is False)
    print(f"  R1/R3/R4 ME verdict (paper: all NOT ME):    {r_me_match}/3 match")
    print(f"  R1/R3/R4 primitivity (mixed paper claims):")
    for (name, claim, prim, me, _) in r_rows:
        print(f"    {name:35s}  paper {claim}; kernel primitive={prim}")
    print(f"  Cattell-Zhang 1995 single CAs primitive:    {cz_matches}/{cz_total}")
    print(f"  Adak-Das CA(90′) primitive:                 {ad_90}/{t_90}")
    print(f"  Adak-Das CA(150′) primitive:                {ad_150}/{t_150}")
    print(f"  Table 1 R1 per-t equi-verdict match:        {t1_match}/{t1_total}")


if __name__ == "__main__":
    main()
