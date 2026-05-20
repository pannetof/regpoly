# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
tvalue_test.py — (t,m,s)-net t-value test.

Computes the profile `t(s)` for `s = 2..s_max` for the digital net
implied by a `Combination` (typically wrapping one `DigitalNet`).
Acts as a search-rejection criterion via the per-s `delta[s]` and
aggregate `max_t_sum` thresholds, mirroring EquidistributionTest's
contract.

Methods:
- `"schmid"` (default): primal Schmid-style enumeration.
- `"niederreiter_pirsic"`: registered name; raises
  ``NotImplementedError`` when run.
"""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

import regpoly._regpoly_cpp as _cpp
from regpoly.analyses.abstract_test import AbstractTest
from regpoly.analyses.tvalue_results import TValueResults

if TYPE_CHECKING:
    from regpoly.core.combination import Combination


METHOD_SCHMID = 0
METHOD_NIEDERREITER_PIRSIC = 1

_STR_TO_METHOD = {
    "schmid": METHOD_SCHMID,
    "niederreiter_pirsic": METHOD_NIEDERREITER_PIRSIC,
}
_METHOD_TO_STR = {v: k for k, v in _STR_TO_METHOD.items()}


class TValueTest(AbstractTest):
    """t-value test configuration.

    Parameters
    ----------
    s_max : int
        Largest dimension to query. Profile is reported for
        ``s = 2..s_max``.
    max_t_sum : int
        Aggregate cap on cumulative `se` (sum of t-values). When
        ``se > max_t_sum`` during the walk, the test bails and
        marks the result unverified.
    delta : list[int] | None
        Per-s cap on ``tvals[s]``, indexed ``0..s_max``. ``None``
        (default) is equivalent to ``[INT_MAX] * (s_max + 1)``.
    tverif : bool
        Master enable flag (mirrors EquidistributionTest's
        ``meverif``). When ``False``, ``run`` returns a trivial
        unverified result without invoking the C++ kernel.
    method : int | str | None
        Method selector: ``"schmid"``/``0`` or
        ``"niederreiter_pirsic"``/``1``. ``None`` (default) resolves
        to ``"schmid"``.
    """

    def __init__(
        self,
        s_max: int,
        max_t_sum: int,
        delta: list[int] | None = None,
        tverif: bool = True,
        method: int | str | None = None,
    ) -> None:
        if s_max < 2:
            raise ValueError(f"TValueTest: s_max must be >= 2 (got {s_max})")
        self.s_max = s_max
        self.max_t_sum = max_t_sum
        self.delta = list(delta) if delta is not None else [sys.maxsize] * (s_max + 1)
        if len(self.delta) < s_max + 1:
            self.delta = self.delta + [sys.maxsize] * (s_max + 1 - len(self.delta))
        self.tverif = tverif

        if method is None:
            self.method = METHOD_SCHMID
        elif isinstance(method, str):
            self.method = _STR_TO_METHOD[method]
        else:
            if method not in (METHOD_SCHMID, METHOD_NIEDERREITER_PIRSIC):
                raise ValueError(f"TValueTest: unknown method {method}")
            self.method = method

    @classmethod
    def _from_params(cls, params: dict, Lmax: int) -> "TValueTest":
        """Build from a YAML parameter dict.

        YAML keys:
            s_max:      int (required)
            max_t_sum:  int (default: sys.maxsize)
            method:     "schmid" | "niederreiter_pirsic" (default: "schmid")
            delta:      list of {from, to, max} per-s caps (optional)
        """
        s_max = params["s_max"]
        max_t_sum = params.get("max_t_sum", sys.maxsize)

        raw_method = params.get("method")
        if raw_method is None or raw_method == "":
            method = METHOD_SCHMID
        else:
            method = _STR_TO_METHOD.get(raw_method, METHOD_SCHMID)

        delta = [sys.maxsize] * (s_max + 1)
        for rule in params.get("delta", []):
            lo = rule["from"]
            hi = rule["to"]
            val = rule["max"]
            for s in range(lo, min(hi, s_max) + 1):
                delta[s] = val

        return cls(s_max=s_max, max_t_sum=max_t_sum,
                   delta=delta, tverif=True, method=method)

    def run(self, C: "Combination", *args, **kwargs) -> TValueResults:
        """Execute the t-value test on combination ``C``.

        For ``method = METHOD_NIEDERREITER_PIRSIC`` the underlying
        C++ kernel raises; the exception propagates so callers learn
        immediately that the dual method is not yet implemented.
        """
        method_str = _METHOD_TO_STR[self.method]

        if not self.tverif:
            return TValueResults(
                s_max=self.s_max,
                tvals=[0] * (self.s_max + 1),
                se=0,
                verified=False,
                delta=self.delta,
                max_t_sum=self.max_t_sum,
                method=method_str,
            )

        INT_MAX = 2**31 - 1
        delta_capped = [min(d, INT_MAX) for d in self.delta]
        max_t_sum_capped = min(self.max_t_sum, INT_MAX)

        gens = [C[j]._cpp_gen for j in range(C.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in C.components
        ]
        combined = _cpp.CombinedGenerator(gens, trans, C.L)

        result = _cpp.run_tvalue_profile(
            combined, C.k_g, C.L, self.s_max,
            delta_capped, max_t_sum_capped, method_str,
        )

        return TValueResults(
            s_max=self.s_max,
            tvals=list(result["tvals"]),
            se=result["se"],
            verified=result["verified"],
            delta=self.delta,
            max_t_sum=self.max_t_sum,
            method=method_str,
        )
