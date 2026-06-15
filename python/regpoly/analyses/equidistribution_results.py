# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
equidistribution_results.py — Results of the ME equidistribution test.
"""

from __future__ import annotations

import math

from regpoly.analyses.abstract_test import AbstractTestResults


def _kg_L(C) -> tuple[int, int]:
    """Extract `(kg, L)` from a ComboEnumerator-shaped or Generator-shaped input.

    `_cpp.CombinedF2LinearSource` exposes `.k()` / `.L()` as methods (no
    `.k_g`). `_cpp.ComboEnumerator` and the legacy Python `ComboEnumerator`
    expose `.k_g` / `.L` as properties/attributes. This helper unifies
    the shapes so the display path doesn't care which kind it got.
    """
    kg = getattr(C, "k_g", None)
    if kg is None:
        kg = C.k()
    L_attr = getattr(C, "L", None)
    L = L_attr() if callable(L_attr) else L_attr
    return int(kg), int(L)


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

    def display(self) -> str:
        """DispME: return ME status string."""
        if self.is_me():
            return "\n ===> ME GENERATOR"
        return ""

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

    def display_table(self, C, by: str = 'l') -> tuple[str, int]:
        """
        DispTable: return (table_string, total_gap_sum).

        Accepts either a `_cpp.CombinedF2LinearSource` (`.k()` / `.L()` are
        methods) or a `_cpp.ComboEnumerator` search iterator (`.k_g` / `.L`
        are read-only properties). Both shapes are handled via the
        `_kg_L` normalization helper.

        by='l' — rows are resolutions l, values are dimension gaps Delta_l
        by='t' — rows are dimensions t, values are resolution gaps Lambda_t
        """
        if not self._verified:
            return ("", -1)

        kg, L = _kg_L(C)

        if by == 'l':
            max_i      = min(L, self.L)
            table      = self.ecart
            row_label  = "RESOL  "
            gap_label  = "ECART  "
            dual_label = "DIM    "
            def dual(i: int) -> int:
                return min(kg, kg // i) - table[i]
        elif by == 't':
            self._set_phi12(C)
            lambda_    = self._conv_ecarts(C)
            table      = [0] + lambda_      # make 1-indexed
            max_i      = kg
            row_label  = "DIM    "
            gap_label  = "ECART  "
            dual_label = "RESOL  "
            def dual(i: int) -> int:
                return min(self.L, kg // i) - table[i]
        else:
            raise ValueError(f"display_table: unknown type '{by}'")

        lines = []
        nblocks = (max_i - 1) // 16 + 1
        for block in range(nblocks):
            start = block * 16 + 1
            end   = min((block + 1) * 16, max_i)
            cols  = list(range(start, end + 1))
            w     = len(cols)
            eqbase = "=======" + "+=====" * w
            mline  = "-------" + "+-----" * w + "|"
            if block == 0:
                lines.append("\n" + eqbase + "+")
            else:
                lines.append(eqbase + "+")
            lines.append(row_label  + "".join(f"|{i:5d}" for i in cols) + "|")
            lines.append(mline)
            lines.append(gap_label  + "".join(
                "|     " if table[i] == 0 else f"|{table[i]:5d}"
                for i in cols
            ) + "|")
            lines.append(mline)
            lines.append(dual_label + "".join(f"|{dual(i):5d}" for i in cols) + "|")
            if block < nblocks - 1:
                lines.append(eqbase)
            else:
                lines.append(eqbase + "+")

        if by == 'l':
            somme = self.se
            lines.append(f"--------------------------->DIMENSION GAPS SUM (Psi_12) = {somme}")
        else:
            somme = sum(
                table[i] for i in range(1, kg + 1) if self._phi12[i]
            )
            lines.append(f"--------------------------->RESOLUTION GAPS SUM (Phi_12) = {somme}")

        return ("\n".join(lines), somme)

    # -- Private helpers --------------------------------------------------

    def _set_phi12(self, C) -> None:
        """SetPhi12: compute Phi_12 — dimensions t to include in the 't' table."""
        kg, _ = _kg_L(C)
        phi  = [False] * (kg + 1)
        r    = int(math.isqrt(kg))
        m    = kg // self.L
        if m < 2:
            m = 2
        if m > kg:
            m = kg
        if r > kg:
            r = kg
        for t in range(m, r + 1):
            phi[t] = True
        r2 = int(math.isqrt(kg - 1))
        for l in range(1, r2 + 1):
            phi[min(kg, kg // l)] = True
        self._phi12 = phi

    def _conv_ecarts(self, C) -> list[int]:
        """ConvEcarts: convert dimension gaps Delta_l to resolution gaps Lambda_t."""
        kg, _ = _kg_L(C)
        lam  = [-1] * (kg + 1)
        for t in range(1, kg + 1):
            l = min(kg // t, self.L)
            for i in range(1, l + 1):
                t_i = min(kg // i, kg) - self.ecart[i]
                if t <= t_i:
                    lam[t] = l - i
        return lam
