"""
equidistribution_results.py — Results of the ME equidistribution test.
"""

from __future__ import annotations

import math
import sys
from typing import TYPE_CHECKING

from regpoly.tests.test_results_base import AbstractTestResults

if TYPE_CHECKING:
    from regpoly.combinaison import Combinaison


class EquidistributionResults(AbstractTestResults):
    """
    Results of one run of EquidistributionTest.

    Attributes
    ----------
    L        : int        — maximum resolution that was tested
    ecart    : list[int]  — dimension gap Delta_l for l in 1..L
    psi12    : list[bool] — resolutions l that were in Psi_12 (tested)
    se       : int        — sum of dimension gaps over Psi_12
    mse      : int        — quasi-ME threshold (from test params)
    meverif  : bool       — whether the test was enabled (from test params)
    delta    : list[int]  — per-resolution bounds (from test params); delta[0] unused
    _verified: bool       — True iff the test actually ran
    """

    def __init__(
        self,
        L: int,
        ecart: list[int],
        psi12: list[bool],
        se: int,
        verified: bool,
        mse: int,
        meverif: bool,
        delta: list[int],
    ) -> None:
        self.L        = L
        self.ecart    = ecart
        self.psi12    = psi12
        self.se       = se
        self._verified = verified
        self.mse      = mse
        self.meverif  = meverif
        self.delta    = delta
        self._phi12: list[bool] = []    # computed lazily by display_table

    # -- AbstractTestResults interface ------------------------------------

    @property
    def verified(self) -> bool:
        return self._verified

    def display(self) -> None:
        """DispME: print ME status."""
        if self.is_me():
            print("\n ===> ME GENERATOR")

    # -- Status predicates ------------------------------------------------

    def is_me(self) -> bool:
        """True iff all dimension gaps are zero."""
        return self._verified and self.se == 0

    def is_quasi_me(self) -> bool:
        """True iff the sum of gaps does not exceed mse."""
        return self._verified and self.se <= self.mse

    def is_presque_me(self) -> bool:
        """True iff quasi-ME and every individual gap is within delta[l]."""
        if not self.meverif:
            return True
        if not self._verified:
            return False
        for l in range(1, self.L + 1):
            if self.psi12[l] and self.ecart[l] > self.delta[l]:
                return False
        return self.se <= self.mse

    # -- Display table ----------------------------------------------------

    def display_table(self, C: "Combinaison", by: str = 'l') -> int:
        """
        DispTable: print the gap table and return the total gap sum.

        by='l' — rows are resolutions l, values are dimension gaps Delta_l
        by='t' — rows are dimensions t, values are resolution gaps Lambda_t
        """
        if not self._verified:
            return -1

        smax = min(C.k_g, C.smax)

        if by == 'l':
            max_i      = min(C.L, self.L)
            table      = self.ecart
            row_label  = "RESOL  "
            gap_label  = "ECART  "
            dual_label = "DIM    "
            def dual(i: int) -> int:
                return min(smax, C.k_g // i) - table[i]
        elif by == 't':
            self._set_phi12(C)
            lambda_    = self._conv_ecarts(C)
            table      = [0] + lambda_      # make 1-indexed
            max_i      = smax
            row_label  = "DIM    "
            gap_label  = "ECART  "
            dual_label = "RESOL  "
            def dual(i: int) -> int:
                return min(self.L, C.k_g // i) - table[i]
        else:
            raise ValueError(f"display_table: unknown type '{by}'")

        nblocks = (max_i - 1) // 16 + 1
        for block in range(nblocks):
            start = block * 16 + 1
            end   = min((block + 1) * 16, max_i)
            cols  = list(range(start, end + 1))
            w     = len(cols)
            eqbase = "=======" + "+=====" * w
            mline  = "-------" + "+-----" * w + "|"
            # Opening =====: trailing +; newline before first block only
            if block == 0:
                print("\n" + eqbase + "+")
            else:
                print(eqbase + "+")
            print(row_label  + "".join(f"|{i:5d}" for i in cols) + "|")
            print(mline)
            print(gap_label  + "".join(
                "|     " if table[i] == 0 else f"|{table[i]:5d}"
                for i in cols
            ) + "|")
            print(mline)
            print(dual_label + "".join(f"|{dual(i):5d}" for i in cols) + "|")
            # Closing =====: no trailing + (next block or final "+\n" adds it)
            if block < nblocks - 1:
                print(eqbase)
            else:
                print(eqbase + "+")

        if by == 'l':
            somme = self.se
            print(f"--------------------------->DIMENSION GAPS SUM (Psi_12) = {somme}")
        else:
            somme = sum(
                table[i] for i in range(1, smax + 1) if self._phi12[i]
            )
            print(f"--------------------------->RESOLUTION GAPS SUM (Phi_12) = {somme}")

        return somme

    # -- Private helpers --------------------------------------------------

    def _set_phi12(self, C: "Combinaison") -> None:
        """SetPhi12: compute Phi_12 — dimensions t to include in the 't' table."""
        smax = min(C.k_g, C.smax)
        phi  = [False] * (smax + 1)
        r    = int(math.isqrt(C.k_g))
        m    = C.k_g // self.L
        if m < 2:
            m = 2
        if m > smax:
            m = smax
        if r > smax:
            r = smax
        for t in range(m, r + 1):
            phi[t] = True
        r2 = int(math.isqrt(C.k_g - 1))
        for l in range(1, r2 + 1):
            phi[min(smax, C.k_g // l)] = True
        self._phi12 = phi

    def _conv_ecarts(self, C: "Combinaison") -> list[int]:
        """ConvEcarts: convert dimension gaps Delta_l to resolution gaps Lambda_t."""
        smax = min(C.k_g, C.smax)
        lam  = [-1] * (smax + 1)
        for t in range(1, smax + 1):
            l = min(C.k_g // t, self.L)
            for i in range(1, l + 1):
                t_i = min(C.k_g // i, smax) - self.ecart[i]
                if t <= t_i:
                    lam[t] = l - i
        return lam
