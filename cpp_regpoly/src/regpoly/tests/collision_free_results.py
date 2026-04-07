"""
collision_free_results.py — Results of the CF collision-free test.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from regpoly.tests.test_results_base import AbstractTestResults

if TYPE_CHECKING:
    from regpoly.tests.equidistribution_results import EquidistributionResults


class CollisionFreeResults(AbstractTestResults):
    """
    Results of one run of CollisionFreeTest.

    Attributes
    ----------
    ecart_cf  : list[int] — CF rank deficit at resolution l (indexed 0..k_g)
    secf      : int       — sum of all CF rank deficits
    msecf     : int       — quasi-CF threshold (from test params)
    _verified : bool      — True iff the test actually ran
    """

    def __init__(
        self,
        ecart_cf: list[int],
        secf: int,
        verified: bool,
        msecf: int,
    ) -> None:
        self.ecart_cf  = ecart_cf
        self.secf      = secf
        self._verified = verified
        self.msecf     = msecf

    # -- AbstractTestResults interface ------------------------------------

    @property
    def verified(self) -> bool:
        return self._verified

    def display(self) -> None:
        """DispCF: print CF status."""
        if self.is_cf_unconditional():
            print("\n ===> CF GENERATOR")

    # -- Status predicates ------------------------------------------------

    def is_cf(self, me_results: "EquidistributionResults") -> bool:
        """True iff ME and all CF rank deficits are zero."""
        return self._verified and me_results.is_me() and self.secf == 0

    def is_quasi_cf(self, me_results: "EquidistributionResults") -> bool:
        """True iff ME and the sum of CF rank deficits does not exceed msecf."""
        return self._verified and me_results.is_me() and self.secf <= self.msecf

    def is_cf_unconditional(self) -> bool:
        """True iff all CF rank deficits are zero (no ME check)."""
        return self._verified and self.secf == 0
