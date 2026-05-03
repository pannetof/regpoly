"""
tuplets_test.py — TupletsTest: Delta(t_1,...,t_d) equidistribution criterion.
"""

from __future__ import annotations

from regpoly.analyses.abstract_test import AbstractTest
from regpoly.analyses.tuplets_results import _MAX_TYPE, _SUM_TYPE, TupletsResults
from regpoly.core.combination import Combination


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

    def run(self, C: Combination, *args, **kwargs) -> TupletsResults:
        """
        TestTuplets: compute Delta(t_1,...,t_d). Phase 2.3: orchestration
        loop now lives in C++ as run_tuplets.
        """
        if not self.tupletsverif:
            return TupletsResults(
                tupletsverif=False, tupd=0, tuph=[], gap=[], DELTA=[],
                pourcentage=[], firstpart_max=0.0, firstpart_sum=0.0,
                secondpart_max=0.0, secondpart_sum=0.0,
                treshold=self.mDD, testtype=self.testtype, verified=False,
            )

        import regpoly._regpoly_cpp as _cpp

        gens = [C[j]._cpp_gen for j in range(C.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in C.components
        ]
        combined = _cpp.CombinedGenerator(gens, trans, C.L)

        result = _cpp.run_tuplets(
            combined, C.k_g, C.L, self.d, list(self.s),
            float(self.mDD), int(self.testtype),
        )

        return TupletsResults(
            tupletsverif=True,
            tupd=result['tupd'],
            tuph=list(result['tuph']),
            gap=list(result['gap']),
            DELTA=list(result['DELTA']),
            pourcentage=list(result['pourcentage']),
            firstpart_max=result['firstpart_max'],
            firstpart_sum=result['firstpart_sum'],
            secondpart_max=result['secondpart_max'],
            secondpart_sum=result['secondpart_sum'],
            treshold=self.mDD,
            testtype=self.testtype,
            verified=True,
        )
