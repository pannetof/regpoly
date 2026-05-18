"""Confirm what equidistribution method the kernel actually used for
the paper-reproduction tests, and whether switching method changes the
verdict on the (31, 32, s=7) and (31, 32, s=8) cases (Tables 9, 10).

We compare:
  - 'matricial'         (current default for CellularAutomataGen, and the
                         method that propagates up through CombinedGenerator
                         via the unanimous-component path)
  - 'notprimitive'      (the spec-compliant method for reducible combined
                         characteristic polynomials — does BM → factor χ_f
                         → invariant subspace → matricial DE core → guard
                         → verify, per docs/theory/equidistribution-spec.md)
  - 'simd_notprimitive' (SIMD variant of 'notprimitive')
"""

from __future__ import annotations

from regpoly.analyses.equidistribution_test import (
    EquidistributionTest,
    METHOD_MATRICIAL, METHOD_NOTPRIMITIVE, METHOD_SIMD_NOTPRIMITIVE,
)
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator

_METHOD_BY_NAME = {
    "matricial":         METHOD_MATRICIAL,
    "notprimitive":      METHOD_NOTPRIMITIVE,
    "simd_notprimitive": METHOD_SIMD_NOTPRIMITIVE,
}


def make_comb(k1, p1, k2, p2, s):
    g1 = Generator.create("CellularAutomataGen", L=min(k1, 64),
                          k=k1, rule150_positions=p1, s=s)
    g2 = Generator.create("CellularAutomataGen", L=min(k2, 64),
                          k=k2, rule150_positions=p2, s=s)
    comb = Combination.CreateFromFiles([[g1], [g2]], Lmax=64,
                                       temperings=[[], []])
    next(iter(comb))
    return comb


def run_one(comb, method):
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=method)
    res = test.run(comb)
    return res.is_me(), tuple(res.ecart), res.se


def heading(s):
    print(f"\n{'=' * 78}")
    print(s)
    print(f"{'=' * 78}")


def report_pair(label, k1, p1, k2, p2, s):
    heading(f"{label}: combined (k1={k1}, r150={p1}) ⊕ (k2={k2}, r150={p2}) at s={s}")
    comb = make_comb(k1, p1, k2, p2, s)
    print(f"  Combined L = {comb.L}, k_g = {comb.k_g}")
    print(f"  Default test method (per kernel) = {comb[0]._cpp_gen.default_test_method('equidistribution')!r}  (comp 1)")

    results = {}
    for name, code in _METHOD_BY_NAME.items():
        comb = make_comb(k1, p1, k2, p2, s)
        try:
            me, ecart, se = run_one(comb, code)
            results[name] = (me, sum(e for e in ecart if e > 0), se)
            print(f"  method={name:18s}  ME={me}  ecart_sum={sum(e for e in ecart if e > 0):4d}  se={se}")
        except Exception as e:
            results[name] = e
            print(f"  method={name:18s}  ERROR: {e}")
    return results


def main():
    cases = [
        ("Table 9 (paper claim: ME)",  31, [10], 32, [0, 14], 7),
        ("Table 10 (paper claim: ME)", 31, [10], 32, [0, 14], 8),
        ("Table 12 row (31,40,s=8) — paper claim: ME", 31, [10], 40, [7], 8),
        ("Table 3 (paper claim: NOT ME)", 31, [10], 32, [0, 14], 1),
        ("Table 5 row (37,42,s=8) — paper claim: ME", 37, [8], 42, [18], 8),
        ("Table 6 row (31,33,s=8) — paper claim: ME", 31, [10], 33, [0], 8),
    ]
    for args in cases:
        report_pair(*args)


if __name__ == "__main__":
    main()
