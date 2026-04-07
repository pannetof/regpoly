"""
collision_free_test.py — CF collision-free test.

Must be run after EquidistributionTest because the CF property requires ME.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from regpoly.tests.test_base import AbstractTest
from regpoly.tests.collision_free_results import CollisionFreeResults

if TYPE_CHECKING:
    from regpoly.combinaison import Combinaison
    from regpoly.tests.equidistribution_results import EquidistributionResults


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

    @classmethod
    def _from_params(cls, params: dict, Lmax: int) -> "CollisionFreeTest":
        """
        YAML params:
            max_gap_sum: int — msecf threshold (default: 0)
        """
        return cls(msecf=params.get("max_gap_sum", 0), cfverif=True)

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
        """
        if not self.cfverif:
            return CollisionFreeResults(
                ecart_cf=[], secf=0, verified=False, msecf=self.msecf
            )

        L = me_results.L if me_results is not None else C.L

        phi4     = self._compute_phi4(C, L)
        ecart_cf = [0] * (C.k_g + 1)
        secf     = 0

        mat_full = self._prepare_mat(C, C.k_g)

        for t in range(C.k_g, 1, -1):
            if phi4[t]:
                l       = C.k_g // t
                mat     = mat_full.copy()
                rank    = self._rang_cf(mat, C.k_g, t, l + 1, C.L)
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
        """
        phi4 = [False] * (C.k_g + 1)
        for t in range(2, C.k_g + 1):
            if C.k_g % t:
                lt = C.k_g // t
                if lt < L and C.k_g % (lt + 1) and C.k_g // (t - 1) > lt:
                    phi4[t] = True
        return phi4

    @staticmethod
    def _rang_cf(mat, kg: int, t: int, l: int, L: int) -> int:
        return mat.rang_cf(kg, t, l, L)
