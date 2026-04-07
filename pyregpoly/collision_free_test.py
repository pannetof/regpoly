"""
collision_free_test.py — CF collision-free test.

Must be run after EquidistributionTest because the CF property requires ME.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from test_base import AbstractTest
from collision_free_results import CollisionFreeResults

if TYPE_CHECKING:
    from combinaison import Combinaison
    from equidistribution_results import EquidistributionResults


class CollisionFreeTest(AbstractTest):
    """
    CF test configuration and algorithm.

    Holds the test parameters; calling run(C, me_results) executes the
    test and returns a CollisionFreeResults object.

    Attributes
    ----------
    msecf   : int  — maximum allowed sum of CF rank deficits (quasi-CF threshold)
    cfverif : bool — whether the test is enabled
    """

    def __init__(self, msecf: int, cfverif: bool = True) -> None:
        self.msecf   = msecf
        self.cfverif = cfverif

    # -- AbstractTest interface -------------------------------------------

    def run(
        self,
        C: "Combinaison",
        *args,
        me_results: "EquidistributionResults | None" = None,
        **kwargs,
    ) -> CollisionFreeResults:
        """
        TestCF: compute the CF rank deficit for each t in Phi_4.

        Builds a full k_g × k_g generator matrix (using _prepare_mat from
        AbstractTest), then for each dimension t in Phi_4 computes the rank
        of the (k_g × t) sub-matrix using l+1 bits per entry (where
        l = k_g // t).  The deficit is k_g minus that rank.

        me_results is accepted as a keyword argument for callers that supply
        it; it is not used during computation (only needed for is_cf()).
        """
        if not self.cfverif:
            return CollisionFreeResults(
                ecart_cf=[], secf=0, verified=False, msecf=self.msecf
            )

        # L needed for Phi_4: use me_results.L if available, else C.L.
        L = me_results.L if me_results is not None else C.L

        phi4     = self._compute_phi4(C, L)
        ecart_cf = [0] * (C.k_g + 1)
        secf     = 0

        mat      = self._prepare_mat(C, C.k_g)
        raw_full = np.asarray(mat._mat, dtype=np.int64)    # (k_g, k_g * L)

        for t in range(C.k_g, 1, -1):
            if phi4[t]:
                l       = C.k_g // t
                raw     = raw_full[:C.k_g, :t * C.L].copy()
                rank    = self._rang_cf(raw, C.k_g, t, l + 1, C.L)
                gap     = C.k_g - rank
                ecart_cf[l]  = gap
                secf        += gap

        return CollisionFreeResults(
            ecart_cf=ecart_cf, secf=secf, verified=True, msecf=self.msecf
        )

    # -- Private helpers --------------------------------------------------

    @staticmethod
    def _compute_phi4(C: "Combinaison", L: int) -> list[bool]:
        """
        SetPhi4: compute Phi_4 — dimensions t whose CF rank must be checked.

        t is in Phi_4 iff k_g % t != 0, l = k_g // t < L,
        k_g % (l+1) != 0, and k_g // (t-1) > l.
        """
        phi4 = [False] * (C.k_g + 1)
        for t in range(2, C.k_g + 1):
            if C.k_g % t:
                lt = C.k_g // t
                if lt < L and C.k_g % (lt + 1) and C.k_g // (t - 1) > lt:
                    phi4[t] = True
        return phi4

    @staticmethod
    def _rang_cf(
        raw: np.ndarray,
        kg: int,
        t: int,
        l: int,
        L: int,
    ) -> int:
        """
        RangCF: Gaussian elimination for the CF rank check.

        raw  : (kg, t*L) int64 array — modified in place; caller must pass a copy
        t    : number of column groups
        l    : number of bits per column group to examine (= k_g // t + 1)
        L    : full bit width per column group in raw
        Returns the rank found.  Stops early when full rank kg is reached
        or no pivot is found in a column.
        """
        rang = 0
        for j in range(t):
            for cl in range(l):
                col = j * L + cl        # bit cl (MSB-first) of column group j
                i   = rang
                while i < kg and raw[i, col] == 0:
                    i += 1
                if i < kg:              # pivot found
                    if i != rang:
                        raw[[rang, i]] = raw[[i, rang]]
                    lo = j * L
                    hi = t * L
                    for ii in range(rang + 1, kg):
                        if raw[ii, col]:
                            raw[ii, lo:hi] ^= raw[rang, lo:hi]
                    rang += 1
                    if rang == kg:
                        return rang     # full rank reached
                else:
                    return rang         # no pivot: return current rank
        return rang
