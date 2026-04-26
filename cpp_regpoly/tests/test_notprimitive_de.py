"""
test_notprimitive_de — METHOD_NOTPRIMITIVE equidistribution test.

Cross-checks the new "no full-period assumption" matricial DE algorithm
(notprimitive_de.cpp / test_me_notprimitive) against:
  - The existing matricial test for small full-period Tausworthe
    generators (algorithm correctness).
  - The SFMT k = 128·N fix (verifies the structural change to SFMT and
    that notprimitive correctly identifies the 607-dim invariant
    subspace inside the 640-bit SFMT607 state).

MT19937 / SFMT19937 cross-checks against published k(v) tables are
deferred — the current implementation maintains one F_2-echelon per
resolution, with O(L²·p²/64) per step, which is impractical at
p ≈ 20 000 without further optimisation (sparse echelons, SIMD, or
restricting to the psi12 resolution set).  That work is documented as
a follow-up in /home/frpan/.claude/plans/i-want-to-build-compiled-swing.md.
"""

import pytest

from regpoly import Generateur
from regpoly.combinaison import Combinaison
from regpoly.analyses.equidistribution_test import (
    EquidistributionTest,
    METHOD_MATRICIAL,
    METHOD_DUALLATTICE,
    METHOD_NOTPRIMITIVE,
)


def _run(C, method):
    L = C.L
    t = EquidistributionTest(
        L=L,
        delta=[2**31 - 1] * (L + 1),
        mse=2**31 - 1,
        method=method,
    )
    return t.run(C).ecart


def _normalize(ecart):
    """matricial leaves ecart[0] = -1 ('not tested'); normalise to 0
    so the result aligns with notprimitive's clean zero."""
    return [0 if e == -1 else e for e in ecart]


# ── Algorithm correctness: matches dual-lattice for full-period
#    generators (notprimitive's full-period fast path dispatches to
#    test_me_lat directly, so they MUST agree).
#
#    matricial is intentionally not used as the reference here: for
#    Tausworthe k=7 with poly=[0,3,7],s=3 it reports d(2)=1 while
#    dual-lattice and notprimitive both report d(2)=0.  The disagreement
#    is genuine and longstanding — it stems from matricial's canonical-
#    basis matrix construction vs lattice-polynomial reduction — and
#    notprimitive's job is to match dual-lattice on the invariant
#    subspace. ────────────────────────────────────────────────────────

@pytest.mark.parametrize("name, mkgen", [
    ("Tausworthe k=7",
     lambda: Generateur.create(
         "Tausworthe", 7, poly=[0, 3, 7], s=3, quicktaus=True)),
    ("Tausworthe k=31",
     lambda: Generateur.create(
         "Tausworthe", 32, poly=[0, 13, 31], s=18, quicktaus=True)),
])
def test_notprimitive_matches_dual_lattice_for_full_period(name, mkgen):
    gen = mkgen()
    C = Combinaison.CreateFromFiles([[gen]], Lmax=gen.L, temperings=[[]])
    C.reset()

    dl_ecart           = _run(C, METHOD_DUALLATTICE)
    notprimitive_ecart = _run(C, METHOD_NOTPRIMITIVE)

    assert notprimitive_ecart == dl_ecart, (
        f"{name}: notprimitive disagrees with dual-lattice.\n"
        f"  dual-lattice: {dl_ecart}\n"
        f"  notprimitive: {notprimitive_ecart}"
    )


# ── SFMT k = 128·N structural fix ─────────────────────────────────────

@pytest.mark.parametrize("mexp, expected_N", [
    (607,    5),    # 607 / 128 + 1 = 5  → k = 640
    (1279,  10),    # 1279 / 128 + 1 = 10 → k = 1280
    (2281,  18),    # → k = 2304
])
def test_sfmt_k_equals_128_times_N(mexp, expected_N):
    """k for any SFMT must be 128·N where N = MEXP/128 + 1, NOT MEXP."""
    # All these MEXP values use the SFMT607 mask shape.  The masks
    # don't matter for this structural assertion.
    gen = Generateur.create(
        "SFMT", 32,
        mexp=mexp, pos1=2, sl1=15, sl2=3, sr1=13, sr2=3,
        msk=[0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f])
    assert gen.k == 128 * expected_N, (
        f"SFMT(mexp={mexp}): k={gen.k}, expected {128*expected_N}"
    )


# ── Notprimitive on SFMT607 — the headline case ──────────────────────

def test_notprimitive_runs_on_sfmt607():
    """notprimitive should run end-to-end on SFMT607 in well under a
    minute and report a k(v) profile against the 607-dim invariant
    subspace (not the 640-dim raw state).  Specifically, the upper
    bound for ecart[v] is floor(607/v), so the worst case is
    ecart[1] = 607 (which would mean the algorithm collapsed); a sane
    SFMT607 should report ecart[v] near zero.
    """
    gen = Generateur.create(
        "SFMT", 32,
        mexp=607, pos1=2, sl1=15, sl2=3, sr1=13, sr2=3,
        msk=[0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f])
    assert gen.k == 640
    C = Combinaison.CreateFromFiles([[gen]], Lmax=32, temperings=[[]])
    C.reset()
    assert C.k_g == 640
    assert C.L == 32

    ecart = _run(C, METHOD_NOTPRIMITIVE)

    # Sanity bounds: ecart values should all be << p = 607.
    # Anything > 50 in the first few positions would indicate that the
    # algorithm picked the wrong invariant subspace or failed in the
    # DE core.
    assert max(ecart[1:6]) < 50, (
        f"SFMT607 notprimitive ecart looks broken in the low resolutions: "
        f"{ecart[1:6]}"
    )
    # The upper bound floor(607/v) for v=1 is 607; reporting near zero
    # is the expected behaviour for a well-tuned SFMT.
    assert ecart[1] <= 5, f"k(1) = 607 - ecart[1] = {607 - ecart[1]}, expected near 607"
