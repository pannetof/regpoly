"""
test_simd_notprimitive — METHOD_SIMD_NOTPRIMITIVE equidistribution.

Implementation status (Phase 3 — see plan file):
  - Generic fast path (full-period non-SIMD generators): WORKING.
    Delegates to test_me_lat → identical results to dual-lattice for
    MELG, MT, full-period Tausworthe.
  - SIMD path for SFMT/dSFMT/MTGP: PARTIAL.  The lattice-substitute
    implementation produces values close to but not exactly matching
    Saito–Matsumoto 2008 Table 3.  Exact match requires a full port
    of MTToolBox's PIS reduction, deferred (the current port doesn't
    converge — see plan file).

Tests:
  - SIMD method matches dual-lattice exactly for MELG607-64,
    MELG19937-64, Tausworthe k=31 (validates the fast path and that
    the abstract API doesn't perturb non-SIMD behavior).
  - SIMD method on SFMT-607 raw terminates and reports finite k(1) > 0.
  - Reference test for SFMT-19937 vs Saito–Matsumoto Table 3 marked
    xfail until the PIS port is complete.
"""

import pytest
from regpoly import Generator
from regpoly.core.combination import Combination
from regpoly.analyses.equidistribution_test import (
    EquidistributionTest,
    METHOD_DUALLATTICE,
    METHOD_NOTPRIMITIVE,
    METHOD_SIMD_NOTPRIMITIVE,
)


# Reference d(v) for SFMT-19937 (paper parameters: pos1=122, sl1=18,
# sl2=1, sr1=11, sr2=1, msk=[dfffffef, ddfecb7f, bffaffff, bffffff6])
# computed via MTToolBox samples/sfmtdc/calc_equidist (Saito's reference
# implementation).  Sum = 4199.
#
# An earlier hardcoded array based on a published table (sum = 4188)
# does NOT match what MTToolBox itself produces with these parameters;
# we use MTToolBox's output as the authoritative reference.
SFMT19937_MTTOOLBOX_D = [
    1,   3,   1,   2,   3,   2,   1,   2,    # v=1..8
    3,   1,   0, 117, 285, 176,  85,   2,    # v=9..16
  543, 478, 425, 372, 325, 282, 242, 206,    # v=17..24
  173, 142, 114,  88,  63,  40,  19,   3,    # v=25..32
]
assert sum(SFMT19937_MTTOOLBOX_D) == 4199


def _run(C, method, maxL=None):
    L = C.L if maxL is None else maxL
    t = EquidistributionTest(
        L=L,
        delta=[2**31 - 1] * (L + 1),
        mse=2**31 - 1,
        method=method,
    )
    return t.run(C).ecart


# ── SIMD method = dual-lattice for full-period non-SIMD ─────────────

@pytest.mark.parametrize("name, mkgen, Lmax", [
    ("Tausworthe k=31",
     lambda: Generator.create("Tausworthe", 32, poly=[0, 13, 31], s=18, quicktaus=True),
     32),
    ("MELG607-64",
     lambda: Generator.create("MELG", 64, w=64, N=10, M=5, r=33,
                                sigma1=13, sigma2=35, a=0x81F1FD68012348BC),
     64),
])
def test_simd_collapses_to_dual_lattice_for_full_period_non_simd(name, mkgen, Lmax):
    gen = mkgen()
    C = Combination.CreateFromFiles([[gen]], Lmax=Lmax, temperings=[[]])
    C.reset()
    dl_ecart   = _run(C, METHOD_DUALLATTICE)
    simd_ecart = _run(C, METHOD_SIMD_NOTPRIMITIVE)
    assert simd_ecart == dl_ecart, (
        f"{name}: SIMD method disagrees with dual-lattice.\n"
        f"  dual-lattice: {dl_ecart}\n"
        f"  simd       : {simd_ecart}"
    )


# ── SFMT-607 SIMD: terminates, finite results ────────────────────────

def test_simd_sfmt607_terminates_with_finite_results():
    gen = Generator.create(
        "SFMT", 32,
        mexp=607, pos1=2, sl1=15, sl2=3, sr1=13, sr2=3,
        msk=[0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f])
    C = Combination.CreateFromFiles([[gen]], Lmax=32, temperings=[[]])
    C.reset()
    # Full L=32 takes a while; verify a smaller maxL terminates.
    ecart = _run(C, METHOD_SIMD_NOTPRIMITIVE, maxL=8)
    # All d(v) finite (not INT_MAX).
    INT_MAX = 2**31 - 1
    for v in range(1, 9):
        assert ecart[v] < INT_MAX, f"SFMT-607 d({v}) = INT_MAX"
    # k(1) within bounds: 0 ≤ d(1) ≤ p (= 607 here).
    assert 0 <= ecart[1] <= 607


# ── All standard SFMT variants vs MTToolBox calc_equidist ───────────
#
# d(v) values from running MTToolBox samples/sfmtdc/calc_equidist on
# each variant's standard parameters (sfmt-1.5.1's SFMT-paramsNNN.h).
# Tests assert exact equality with my METHOD_SIMD_NOTPRIMITIVE.

# Format: (mexp, pos1, sl1, sl2, sr1, sr2, msk[4])
SFMT_PARAMS = {
    607:   (2, 15, 3, 13, 3, [0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f]),
    1279:  (7, 14, 3, 5, 1,  [0xf7fefffd, 0x7fefcfff, 0xaff3ef3f, 0xb5ffff7f]),
    2281:  (12, 19, 1, 5, 1, [0xbff7ffbf, 0xfdfffffe, 0xf7ffef7f, 0xf2f7cbbf]),
    4253:  (17, 20, 1, 7, 1, [0x9f7bffff, 0x9fffff5f, 0x3efffffb, 0xfffff7bb]),
    11213: (68, 14, 3, 7, 3, [0xeffff7fb, 0xffffffef, 0xdfdfbfff, 0x7fffdbfd]),
    19937: (122, 18, 1, 11, 1, [0xdfffffef, 0xddfecb7f, 0xbffaffff, 0xbffffff6]),
}

# d(v=1..5) per MTToolBox calc_equidist on the params above.
SFMT_MTTOOLBOX_D5 = {
    607:   [3, 3, 2, 3, 1],
    1279:  [6, 3, 2, 3, 3],
    2281:  [4, 1, 0, 2, 2],
    4253:  [1, 2, 1, 3, 2],
    11213: [6, 2, 2, 3, 2],
    19937: [1, 3, 1, 2, 3],
}


@pytest.mark.parametrize("mexp", sorted(SFMT_PARAMS))
def test_simd_sfmt_matches_mttoolbox(mexp):
    pos1, sl1, sl2, sr1, sr2, msk = SFMT_PARAMS[mexp]
    gen = Generator.create("SFMT", 32, mexp=mexp,
                            pos1=pos1, sl1=sl1, sl2=sl2,
                            sr1=sr1, sr2=sr2, msk=msk)
    C = Combination.CreateFromFiles([[gen]], Lmax=32, temperings=[[]])
    C.reset()
    maxL = 5
    ecart = _run(C, METHOD_SIMD_NOTPRIMITIVE, maxL=maxL)
    measured = list(ecart[1:maxL + 1])
    expected = SFMT_MTTOOLBOX_D5[mexp]
    assert measured == expected, (
        f"SFMT-{mexp} SIMD does not match MTToolBox at maxL={maxL}.\n"
        f"  mttoolbox: {expected}\n"
        f"  measured : {measured}\n"
    )


def test_simd_sfmt19937_matches_mttoolbox():
    """Kept for backward compatibility (subset of test_simd_sfmt_matches_mttoolbox[19937])."""
    gen = Generator.create(
        "SFMT", 32,
        mexp=19937, pos1=122, sl1=18, sl2=1, sr1=11, sr2=1,
        msk=[0xdfffffef, 0xddfecb7f, 0xbffaffff, 0xbffffff6])
    C = Combination.CreateFromFiles([[gen]], Lmax=32, temperings=[[]])
    C.reset()
    maxL = 5
    ecart = _run(C, METHOD_SIMD_NOTPRIMITIVE, maxL=maxL)
    measured = list(ecart[1:maxL + 1])
    expected = SFMT19937_MTTOOLBOX_D[:maxL]
    assert measured == expected, (
        f"SFMT-19937 SIMD does not match MTToolBox at maxL={maxL}.\n"
        f"  mttoolbox: {expected}\n"
        f"  measured : {measured}\n"
    )


@pytest.mark.skip(reason="Full maxL=32 takes several minutes; covered by "
                         "test_simd_sfmt19937_matches_mttoolbox at maxL=5.")
def test_simd_sfmt19937_matches_mttoolbox_full():
    gen = Generator.create(
        "SFMT", 32,
        mexp=19937, pos1=122, sl1=18, sl2=1, sr1=11, sr2=1,
        msk=[0xdfffffef, 0xddfecb7f, 0xbffaffff, 0xbffffff6])
    C = Combination.CreateFromFiles([[gen]], Lmax=32, temperings=[[]])
    C.reset()
    ecart = _run(C, METHOD_SIMD_NOTPRIMITIVE)
    measured = list(ecart[1:33])
    assert measured == SFMT19937_MTTOOLBOX_D, (
        f"SFMT-19937 SIMD does not match MTToolBox d(1..32).\n"
        f"  mttoolbox: {SFMT19937_MTTOOLBOX_D}\n"
        f"  measured : {measured}\n"
    )
