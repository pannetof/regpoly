"""
tuplets_results.py — Results of the Tuplets Delta(t_1,...,t_d) criterion.
"""

from __future__ import annotations

import sys

from regpoly.analyses.test_results_base import AbstractTestResults

_MAX_TYPE = 1
_SUM_TYPE = 0


class TupletsResults(AbstractTestResults):
    """
    Results of one run of TupletsTest.

    Attributes
    ----------
    tupletsverif      : bool      — whether the test was enabled
    tupd              : int       — depth d of the criterion
    tuph              : list[int] — dimensions t[1..d], 1-indexed
    gap               : list[float] — successive resolution gaps, 1-indexed
    DELTA             : list[float] — max gap per non-successive group, 1-indexed
    pourcentage       : list[float] — fraction with zero gap, 1-indexed
    firstpart_max/sum : float     — statistics over successive dimensions
    secondpart_max/sum: float     — statistics over non-successive groups
    treshold          : float     — acceptance threshold
    testtype          : int       — _MAX_TYPE or _SUM_TYPE
    _verified         : bool
    """

    def __init__(
        self,
        tupletsverif: bool,
        tupd: int,
        tuph: list,
        gap: list,
        DELTA: list,
        pourcentage: list,
        firstpart_max: float,
        firstpart_sum: float,
        secondpart_max: float,
        secondpart_sum: float,
        treshold: float,
        testtype: int,
        verified: bool,
    ) -> None:
        self.tupletsverif   = tupletsverif
        self.tupd           = tupd
        self.tuph           = tuph
        self.gap            = gap
        self.DELTA          = DELTA
        self.pourcentage    = pourcentage
        self.firstpart_max  = firstpart_max
        self.firstpart_sum  = firstpart_sum
        self.secondpart_max = secondpart_max
        self.secondpart_sum = secondpart_sum
        self.treshold       = treshold
        self.testtype       = testtype
        self._verified      = verified

    # -- AbstractTestResults interface ------------------------------------

    @property
    def verified(self) -> bool:
        return self._verified

    def display(self) -> str:
        """DispTuplets: return the gap tables as a string."""
        if not self._verified:
            return ""

        lines = []

        # First part: successive dimensions
        lines.append("\n  Tables of gaps obtained (successive dimensions)")
        T    = 1
        maxD = self.tuph[1]
        while maxD > 0:
            cols = min(maxD, 16)
            lines.append("======" + "+=====" * cols + "+")
            lines.append("DIM   " + "".join(f"|{T + j - 1:5d}" for j in range(1, cols + 1)) + "|")
            lines.append("------" + "+-----" * cols + "+")
            lines.append("GAP   " + "".join(
                "|     " if self.gap[T + j - 1] == 0.0
                else f"|{int(self.gap[T + j - 1]):5d}"
                for j in range(1, cols + 1)
            ) + "|")
            lines.append("======" + "+=====" * cols + "+")
            maxD -= 16
            T    += 16
        lines.append("")

        # Second part: non-successive dimensions (if tupd > 1)
        if self.tupd > 1:
            lines.append("\n\n   Tables of maximal gaps obtained (non-successive dimensions)")
            w = self.tupd - 1
            lines.append("===========" + "+" + "======+" * w)
            lines.append("Dimension  |" + "".join(f"{i:6d}|" for i in range(2, self.tupd + 1)))
            lines.append("-----------+" + "------+" * w)
            lines.append("t_i        |" + "".join(f"{self.tuph[i]:6d}|" for i in range(2, self.tupd + 1)))
            lines.append("-----------+" + "------+" * w)
            lines.append("GAP        |" + "".join(f" {int(self.DELTA[i]):5d}|" for i in range(2, self.tupd + 1)))
            lines.append("-----------+" + "------+" * w)
            lines.append("percentage |" + "".join(f"{self.pourcentage[i] * 100:6.2f}|" for i in range(2, self.tupd + 1)))
            lines.append("===========" + "+" + "======+" * w)

        delta_val = max(self.firstpart_max, self.secondpart_max)
        lines.append(f"Value of DELTA( " +
              ", ".join(str(self.tuph[i]) for i in range(1, self.tupd + 1)) +
              f" ) = {delta_val:5.3f}")
        lines.append(f"Sum of all gaps observed = {self.firstpart_sum + self.secondpart_sum:5.3f}")

        return "\n".join(lines)

    # -- Status predicate -------------------------------------------------

    def is_ok(self) -> bool:
        """True iff the test was disabled or the generator passes the criterion."""
        if not self.tupletsverif:
            return True
        if self.testtype == _MAX_TYPE:
            return self.treshold >= max(self.firstpart_max, self.secondpart_max)
        else:  # SUM_TYPE
            return self.treshold >= (self.firstpart_sum + self.secondpart_sum)
