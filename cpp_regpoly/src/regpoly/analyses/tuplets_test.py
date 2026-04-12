"""
tuplets_test.py — TupletsTest: Delta(t_1,...,t_d) equidistribution criterion.
"""

from __future__ import annotations

import sys
from itertools import combinations as _itertools_combinations

from regpoly.analyses.test_base import AbstractTest
from regpoly.analyses.tuplets_results import TupletsResults, _MAX_TYPE, _SUM_TYPE
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.combinaison import Combinaison


class TupletsTest(AbstractTest):
    """
    Delta(t_1,...,t_d) equidistribution criterion for non-successive outputs.

    Attributes
    ----------
    tupletsverif : bool       — whether to run the test
    d            : int        — depth of the criterion
    s            : list[int]  — t values, 1-indexed (s[1..d])
    mDD          : float      — acceptance threshold for the criterion
    testtype     : int        — _MAX_TYPE (default) or _SUM_TYPE
    """

    def __init__(
        self,
        tupletsverif: bool,
        d: int = 0,
        s: list = None,
        mDD: float = 0.0,
        testtype: int = _MAX_TYPE,
    ) -> None:
        self.tupletsverif = tupletsverif
        self.d            = d
        self.s            = s or []   # 1-indexed; s[0] unused
        self.mDD          = mDD
        self.testtype     = testtype

    @classmethod
    def _from_params(cls, params: dict, Lmax: int) -> "TupletsTest":
        """
        YAML params:
            dimensions: list[int] — t values, e.g. [50, 10, 5]
            threshold: float — mDD acceptance threshold (default: 0.0)
            test_type: str — "max" (default) or "sum"
        """
        dims = params.get("dimensions", [])
        d = len(dims)
        # Build 1-indexed s array: s[0] unused, s[1..d] = dimensions
        s = [0] + dims
        threshold = params.get("threshold", 0.0)
        type_str = params.get("test_type", "max")
        testtype = _MAX_TYPE if type_str == "max" else _SUM_TYPE
        return cls(tupletsverif=True, d=d, s=s, mDD=threshold, testtype=testtype)

    # -- AbstractTest interface -------------------------------------------

    def run(self, C: Combinaison, *args, **kwargs) -> TupletsResults:
        """
        TestTuplets: compute Delta(t_1,...,t_d) for the current generator in C.
        """
        if not self.tupletsverif:
            return TupletsResults(
                tupletsverif=False, tupd=0, tuph=[], gap=[], DELTA=[],
                pourcentage=[], firstpart_max=0.0, firstpart_sum=0.0,
                secondpart_max=0.0, secondpart_sum=0.0,
                treshold=self.mDD, testtype=self.testtype, verified=False,
            )

        kg         = C.k_g
        L          = C.L
        tupd       = self.d
        tuph       = self.s          # 1-indexed
        treshold   = self.mDD
        testtype   = self.testtype
        indice_max = max(tuph[1:tupd + 1])

        # Build and keep one working matrix (modified in place across all calls,
        # matching the C code's behaviour with the same tMat object).
        working = self._prepare_mat(C, indice_max)

        # --- First part: successive dimensions dim = 1..tuph[1] ---
        maxindsucc    = min(tuph[1], kg)
        gap           = [0.0] * (tuph[1] + 1)   # 1-indexed
        firstpart_max = float(-sys.maxsize)
        firstpart_sum = 0.0

        for dim in range(1, maxindsucc + 1):
            indices = list(range(dim))
            l_t = EquidistributionTest._resolution_equid(working, kg, dim, L, indices)
            gap[dim] = float(min(L, kg // dim) - l_t)

            firstpart_sum += gap[dim]
            if gap[dim] > firstpart_max:
                firstpart_max = gap[dim]

            if testtype == _MAX_TYPE and firstpart_max > treshold:
                firstpart_max = float(sys.maxsize)
                break
            if testtype == _SUM_TYPE and firstpart_sum > treshold:
                firstpart_sum = float(sys.maxsize)
                break

        # --- Second part: non-successive dimensions dim = 2..d ---
        DELTA       = [float(-sys.maxsize)] * (tupd + 1)   # 1-indexed
        pourcentage = [0.0]                * (tupd + 1)
        secondpart_max = float(-sys.maxsize)
        secondpart_sum = 0.0

        first_ok = (
            (testtype == _MAX_TYPE and firstpart_max <= treshold) or
            (testtype == _SUM_TYPE and firstpart_sum <= treshold)
        )

        if first_ok and tupd > 1:
            for dim_idx in range(2, tupd + 1):
                DELTA[dim_idx] = float(-sys.maxsize)
                bound    = float(min(L, kg // dim_idx))
                nbposs   = 0
                nbeczero = 0
                stop     = False

                for rest in _itertools_combinations(range(1, tuph[dim_idx]), dim_idx - 1):
                    indices = [0] + list(rest)
                    l_t = EquidistributionTest._resolution_equid(working, kg, dim_idx, L, indices)
                    gap_val = bound - l_t

                    if gap_val > DELTA[dim_idx]:
                        DELTA[dim_idx] = gap_val
                    if gap_val == 0.0:
                        nbeczero += 1
                    else:
                        secondpart_sum += gap_val
                    nbposs += 1

                    if testtype == _MAX_TYPE and DELTA[dim_idx] > treshold:
                        secondpart_max = DELTA[dim_idx] = float(sys.maxsize)
                        stop = True
                        break
                    if testtype == _SUM_TYPE and secondpart_sum > treshold:
                        secondpart_sum = DELTA[dim_idx] = float(sys.maxsize)
                        stop = True
                        break

                if nbposs > 0:
                    pourcentage[dim_idx] = nbeczero / nbposs

                if secondpart_max < DELTA[dim_idx]:
                    secondpart_max = DELTA[dim_idx]

                if stop:
                    break
        else:
            secondpart_max = float(sys.maxsize)
            secondpart_sum = float(sys.maxsize)

        return TupletsResults(
            tupletsverif=True,
            tupd=tupd,
            tuph=tuph,
            gap=gap,
            DELTA=DELTA,
            pourcentage=pourcentage,
            firstpart_max=firstpart_max,
            firstpart_sum=firstpart_sum,
            secondpart_max=secondpart_max,
            secondpart_sum=secondpart_sum,
            treshold=treshold,
            testtype=testtype,
            verified=True,
        )
