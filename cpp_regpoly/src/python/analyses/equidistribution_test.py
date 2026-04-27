"""
equidistribution_test.py — ME equidistribution test (METHOD_MATRICIAL).

METHOD_DUALLATTICE raises NotImplementedError if requested.
"""

from __future__ import annotations

import math
import sys
from typing import TYPE_CHECKING

from regpoly.analyses.abstract_test import AbstractTest
from regpoly.analyses.equidistribution_results import EquidistributionResults

if TYPE_CHECKING:
    from regpoly.core.combination import Combination

METHOD_MATRICIAL    = 0
METHOD_DUALLATTICE  = 1
METHOD_NOTHING      = 2
METHOD_HARASE       = 3
METHOD_NOTPRIMITIVE = 4   # matricial DE on BM-recovered invariant subspace
                          # (no full-period assumption)
METHOD_SIMD_NOTPRIMITIVE = 5  # SIMD-aware variant for SFMT/dSFMT/MTGP — same
                              # subspace selection but uses Saito-style PIS
                              # with lane-interleaved super-words to match
                              # Saito-Matsumoto 2008 published k(v) values.


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
        if method not in (METHOD_MATRICIAL, METHOD_DUALLATTICE, METHOD_NOTHING,
                          METHOD_HARASE, METHOD_NOTPRIMITIVE,
                          METHOD_SIMD_NOTPRIMITIVE):
            raise ValueError(f"Unknown method: {method}")
        self.L       = L
        self.delta   = list(delta)      # indexed 0..L; delta[0] unused
        self.mse     = mse
        self.meverif = meverif
        self.method  = method

    @classmethod
    def _from_params(cls, params: dict, Lmax: int) -> "EquidistributionTest":
        """
        YAML params:
            max_gap_sum: int — mse threshold
            method: str — "matricial" (default) or "nothing"
            delta: list of {from, to, max} — per-resolution bounds (optional)
                   default: all delta[l] = sys.maxsize
        """
        import sys as _sys

        mse = params.get("max_gap_sum", _sys.maxsize)
        method_str = params.get("method", "matricial")
        if method_str == "matricial":
            method = METHOD_MATRICIAL
        elif method_str == "lattice":
            method = METHOD_DUALLATTICE
        elif method_str == "harase":
            method = METHOD_HARASE
        elif method_str == "notprimitive":
            method = METHOD_NOTPRIMITIVE
        elif method_str == "simd_notprimitive":
            method = METHOD_SIMD_NOTPRIMITIVE
        else:
            method = METHOD_NOTHING

        # Build delta array indexed 0..Lmax
        delta = [_sys.maxsize] * (Lmax + 1)
        for rule in params.get("delta", []):
            lo = rule["from"]
            hi = rule["to"]
            val = rule["max"]
            for l in range(lo, min(hi, Lmax) + 1):
                delta[l] = val

        return cls(L=Lmax, delta=delta, mse=mse, meverif=True, method=method)

    # -- AbstractTest interface -------------------------------------------

    def run(self, C: "Combination", *args, **kwargs) -> EquidistributionResults:
        """
        Compute dimension gaps Delta_l. Uses either the matricial method
        (Gaussian elimination) or the dual lattice method (Lenstra's algorithm).
        """
        if not self.meverif or self.method == METHOD_NOTHING:
            return EquidistributionResults(
                L=self.L, ecart=[0] * (self.L + 1),
                psi12=[False] * (self.L + 1), se=0,
                verified=False, mse=self.mse,
                meverif=self.meverif, delta=self.delta,
            )

        if self.method == METHOD_DUALLATTICE:
            return self._run_lattice(C)

        if self.method == METHOD_HARASE:
            return self._run_harase(C)

        if self.method == METHOD_NOTPRIMITIVE:
            return self._run_notprimitive(C)

        if self.method == METHOD_SIMD_NOTPRIMITIVE:
            return self._run_simd_notprimitive(C)

        if self.L < C.L:
            raise ValueError(
                f"TestME: EquidistributionTest.L ({self.L}) < generator L ({C.L})"
            )

        ecart = [-1] * (self.L + 1)
        se    = 0
        psi12 = self._compute_psi12(C)

        indice_max = C.k_g
        mat_full = self._prepare_mat(C, indice_max)

        verif = False
        maxl  = self.L
        l     = 1
        while l <= self.L:
            if ecart[l] == -1 and (psi12[l] or verif):
                t   = C.k_g // l
                mat = mat_full.copy()
                t_l = self._dimension_equid(mat, C.k_g, l, C.L)
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

    def _compute_psi12(self, C: "Combination") -> list[bool]:
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
    def _dimension_equid(mat, kg: int, l: int, L: int) -> int:
        return mat.dimension_equid(kg, l, L)

    @staticmethod
    def _resolution_equid(mat, kg: int, t: int, L: int, indices: list) -> int:
        return mat.resolution_equid(kg, t, L, indices)

    def _run_lattice(self, C: "Combination") -> EquidistributionResults:
        """TestMELat: compute dimension gaps via dual lattice basis reduction."""
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        gens = [C[j]._cpp_gen for j in range(C.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in C.components
        ]

        # Cap values to C++ int range
        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)

        result = _cpp.test_me_lat(
            gens, trans,
            C.k_g, C.L, self.L,
            delta_capped, mse_capped,
        )

        psi12 = self._compute_psi12(C)
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    def _run_harase(self, C: "Combination") -> EquidistributionResults:
        """Harase-Matsumoto-Saito: primal lattice with Mulders-Storjohann."""
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        gens = [C[j]._cpp_gen for j in range(C.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in C.components
        ]

        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)

        result = _cpp.test_me_lat_harase(
            gens, trans,
            C.k_g, C.L, self.L,
            delta_capped, mse_capped,
        )

        psi12 = self._compute_psi12(C)
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    def _run_notprimitive(self, C: "Combination") -> EquidistributionResults:
        """Matricial DE on the BM-recovered invariant subspace.

        Does NOT assume the combined characteristic polynomial is primitive
        (i.e. does not assume the generator is full-period).  Reports
        ecart[v] = floor(p/v) - k(v) where p = dim of the certified
        invariant subspace V.  See docs/C.md and notprimitive_de.h.
        """
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        gens = [C[j]._cpp_gen for j in range(C.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in C.components
        ]

        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)

        result = _cpp.test_me_notprimitive(
            gens, trans,
            C.k_g, C.L, self.L,
            delta_capped, mse_capped,
        )

        psi12 = self._compute_psi12(C)
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    def _run_simd_notprimitive(self, C: "Combination") -> EquidistributionResults:
        """SIMD-aware notprimitive — matches Saito-Matsumoto 2008 published
        k(v) for SFMT-like generators (dispatches to test_me_lat for
        non-SIMD full-period generators)."""
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        gens = [C[j]._cpp_gen for j in range(C.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in C.components
        ]

        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)

        result = _cpp.test_me_notprimitive_simd(
            gens, trans,
            C.k_g, C.L, self.L,
            delta_capped, mse_capped,
        )

        psi12 = self._compute_psi12(C)
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )
