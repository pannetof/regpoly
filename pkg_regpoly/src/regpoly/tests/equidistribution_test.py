"""
equidistribution_test.py — ME equidistribution test (METHOD_MATRICIAL).

METHOD_DUALLATTICE raises NotImplementedError if requested.
"""

from __future__ import annotations

import math
import sys
from typing import TYPE_CHECKING

from regpoly.tests.test_base import AbstractTest
from regpoly.tests.equidistribution_results import EquidistributionResults

if TYPE_CHECKING:
    from regpoly.combinaison import Combinaison

METHOD_MATRICIAL   = 0
METHOD_DUALLATTICE = 1
METHOD_NOTHING     = 2


class EquidistributionTest(AbstractTest):
    """
    ME test configuration and algorithm.

    Holds the test parameters; calling run(C) executes the test and
    returns an EquidistributionResults object.

    Attributes
    ----------
    L       : int       — maximum resolution to test
    delta   : list[int] — per-resolution gap bound (delta[l]); delta[0] unused
    mse     : int       — maximum allowed sum of gaps (quasi-ME threshold)
    meverif : bool      — whether the test is enabled
    method  : int       — METHOD_MATRICIAL or METHOD_NOTHING
    """

    def __init__(
        self,
        L: int,
        delta: list[int],
        mse: int,
        meverif: bool = True,
        method: int = METHOD_MATRICIAL,
    ) -> None:
        if method == METHOD_DUALLATTICE:
            raise NotImplementedError("METHOD_DUALLATTICE is not implemented")
        if method not in (METHOD_MATRICIAL, METHOD_NOTHING):
            raise ValueError(f"Unknown method: {method}")
        self.L       = L
        self.delta   = list(delta)      # indexed 0..L; delta[0] unused
        self.mse     = mse
        self.meverif = meverif
        self.method  = method

    # -- AbstractTest interface -------------------------------------------

    def run(self, C: "Combinaison", *args, **kwargs) -> EquidistributionResults:
        """
        TestME (METHOD_MATRICIAL): compute dimension gaps Delta_l for all l
        in Psi_12.
        """
        if not self.meverif or self.method == METHOD_NOTHING:
            return EquidistributionResults(
                L=self.L, ecart=[0] * (self.L + 1),
                psi12=[False] * (self.L + 1), se=0,
                verified=False, mse=self.mse,
                meverif=self.meverif, delta=self.delta,
            )

        if self.L < C.L:
            raise ValueError(
                f"TestME: EquidistributionTest.L ({self.L}) < generator L ({C.L})"
            )

        ecart = [-1] * (self.L + 1)
        se    = 0
        psi12 = self._compute_psi12(C)

        indice_max = min(C.k_g, C.smax)
        mat_full = self._prepare_mat_packed(C, indice_max)

        verif = False
        maxl  = self.L
        l     = 1
        while l <= self.L:
            if ecart[l] == -1 and (psi12[l] or verif):
                t   = min(C.k_g // l, C.smax)
                mat = mat_full.copy()
                t_l = self._dimension_equid(mat, C.k_g, l, C.L, C.smax)
                ecart[l] = t - t_l
                se       += ecart[l]

                if ecart[l] > self.delta[l] or se > self.mse:
                    maxl = l
                    break

                if ecart[l] != 0:
                    verif = True
                    if l != 1:
                        l -= 2
                else:
                    verif = False
            l += 1

        se = 0
        for l in range(1, maxl + 1):
            if ecart[l] == -1:
                ecart[l] = 0
            se += ecart[l]
        for l in range(maxl + 1, self.L + 1):
            if ecart[l] == -1:
                ecart[l] = sys.maxsize

        return EquidistributionResults(
            L=self.L, ecart=ecart, psi12=psi12, se=se,
            verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    # -- Private helpers --------------------------------------------------

    def _compute_psi12(self, C: "Combinaison") -> list[bool]:
        """SetPsi12: compute Psi_12 — resolutions l whose gap must be tested."""
        psi = [False] * (self.L + 1)
        r = int(math.isqrt(C.k_g))
        for l in range(1, min(r, self.L) + 1):
            psi[l] = True
        r2 = int(math.isqrt(C.k_g - 1))
        m  = C.k_g // self.L
        if m < 2:
            m = 2
        for t in range(m, r2 + 1):
            psi[min(C.k_g // t, self.L)] = True
        return psi

    @staticmethod
    def _dimension_equid(mat, kg: int, l: int, L: int, smax: int) -> int:
        """
        DimensionEquid: Gaussian elimination to compute t_l.

        mat : GaussMatrix (PyIntGaussMatrix or PackedGaussMatrix)
        """
        t    = min(kg // l, smax)
        rang = 0
        for j in range(t):
            for cl in range(l):
                col = j * L + cl
                i = mat.find_pivot(col, rang)
                if i >= 0:
                    if i != rang:
                        mat.swap_rows(rang, i)
                    mat.eliminate_column_masked(rang, col, rang + 1,
                                               j * L, t * L)
                    rang += 1
                else:
                    return j
        return t

    @staticmethod
    def _resolution_equid(mat, kg: int, t: int, L: int, indices: list) -> int:
        """
        ResolutionEquid: Gaussian elimination to compute l_t.

        mat : GaussMatrix — modified in place across calls.
        """
        l    = min(kg // t, L)
        rang = 0
        for cl in range(l):
            for j in range(t):
                col = indices[j] * L + cl
                i = mat.find_pivot(col, rang)
                if i >= 0:
                    if i != rang:
                        mat.swap_rows(rang, i)
                    mat.eliminate_column(rang, col, rang + 1)
                    rang += 1
                else:
                    return cl
        return l
