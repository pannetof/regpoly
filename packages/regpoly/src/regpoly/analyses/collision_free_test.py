# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
collision_free_test.py — CF collision-free test.

Must be run after EquidistributionTest because the CF property requires ME.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from regpoly.analyses.abstract_test import AbstractTest
from regpoly.analyses.collision_free_results import CollisionFreeResults

if TYPE_CHECKING:
    from regpoly.analyses.equidistribution_results import EquidistributionResults
    from regpoly.core.combination import Combination


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
        C: "Combination",
        *args,
        me_results: "EquidistributionResults | None" = None,
        **kwargs,
    ) -> CollisionFreeResults:
        """
        TestCF: compute the CF rank deficit for each t in Phi_4. Phase 2.3:
        the inner loop now runs in C++.
        """
        if not self.cfverif:
            return CollisionFreeResults(
                ecart_cf=[], secf=0, verified=False, msecf=self.msecf
            )

        import regpoly._regpoly_cpp as _cpp

        L_for_phi4 = me_results.L if me_results is not None else C.L

        gens = [C[j]._cpp_gen for j in range(C.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in C.components
        ]
        combined = _cpp.CombinedGenerator(gens, trans, C.L)

        result = _cpp.run_collision_free(combined, C.k_g, C.L, L_for_phi4)
        return CollisionFreeResults(
            ecart_cf=list(result['ecart_cf']), secf=result['secf'],
            verified=result['verified'], msecf=self.msecf,
        )

    # -- Private helpers --------------------------------------------------

    @staticmethod
    def _compute_phi4(C: "Combination", L: int) -> list[bool]:
        """SetPhi4: dimensions whose CF rank must be checked. Delegates to C++."""
        import regpoly._regpoly_cpp as _cpp
        return list(_cpp.compute_phi4(C.k_g, L))

    @staticmethod
    def _rang_cf(mat, kg: int, t: int, l: int, L: int) -> int:
        return mat.rang_cf(kg, t, l, L)
