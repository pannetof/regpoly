# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
equidistribution_test.py — ME equidistribution test (METHOD_MATRICIAL).

METHOD_DUALLATTICE raises NotImplementedError if requested.
"""

from __future__ import annotations

from regpoly.analyses.abstract_test import AbstractTest
from regpoly.analyses.equidistribution_results import EquidistributionResults

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

# Canonical mapping between the string names (used in YAML configs and
# returned by the C++ Generator::default_test_method) and the integer
# constants used by run(). Kept here so YAML parsing and default
# resolution share a single source of truth.
_STR_TO_METHOD = {
    "matricial":         METHOD_MATRICIAL,
    "lattice":           METHOD_DUALLATTICE,
    "nothing":           METHOD_NOTHING,
    "harase":            METHOD_HARASE,
    "notprimitive":      METHOD_NOTPRIMITIVE,
    "simd_notprimitive": METHOD_SIMD_NOTPRIMITIVE,
}


class EquidistributionTest(AbstractTest):
    """
    ME test configuration and algorithm.

    Holds the test parameters; calling run(C) executes the test and
    returns an EquidistributionResults object.

    Attributes
    ----------
    L       : int        — maximum resolution to test
    delta   : list[int]  — per-resolution gap bound (delta[l]); delta[0] unused
    mse     : int        — maximum allowed sum of gaps (quasi-ME threshold)
    meverif : bool       — whether the test is enabled
    method  : int | None — METHOD_* constant, or None to defer resolution
                           to run() via the generator's default_test_method
    """

    def __init__(
        self,
        L: int,
        delta: list[int],
        mse: int,
        meverif: bool = True,
        method: int | None = None,
    ) -> None:
        if method is not None and method not in (
            METHOD_MATRICIAL, METHOD_DUALLATTICE, METHOD_NOTHING,
            METHOD_HARASE, METHOD_NOTPRIMITIVE,
            METHOD_SIMD_NOTPRIMITIVE,
        ):
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
            method: str — one of "matricial", "lattice", "harase",
                          "notprimitive", "simd_notprimitive", "nothing".
                          Omit to defer to the generator's default
                          (resolved at run() time).
            delta: list of {from, to, max} — per-resolution bounds (optional)
                   default: all delta[l] = sys.maxsize
        """
        import sys as _sys

        mse = params.get("max_gap_sum", _sys.maxsize)
        # Treat None / missing / empty-string as "no method specified" —
        # run() will resolve via the generator's default_test_method.
        # This avoids a silent METHOD_NOTHING fall-through when the
        # caller passes method='' (e.g. from a form whose default
        # was never populated).
        raw_method = params.get("method")
        if raw_method is None or raw_method == "":
            method = None
        else:
            method = _STR_TO_METHOD.get(raw_method)
            if method is None:
                # Unknown non-empty string: treat as METHOD_NOTHING for
                # backward compatibility with prior YAML behaviour.
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

    def run(self, gen_or_C, *args, **kwargs) -> EquidistributionResults:
        """
        Compute dimension gaps Delta_l. Accepts a Generator (any subclass,
        including `CombinedF2LinearSource` or a `DigitalNet`) or — during the
        Phase-1 migration window — a Python `Combination`. Uses either the
        matricial method (Gaussian elimination) or one of the lattice
        methods (dual / harase / notprimitive / simd_notprimitive).
        """
        cpp_gen = self._to_cpp_gen(gen_or_C)
        method = (self._resolve_method(cpp_gen)
                  if self.method is None else self.method)

        if not self.meverif or method == METHOD_NOTHING:
            return EquidistributionResults(
                L=self.L, ecart=[0] * (self.L + 1),
                psi12=[False] * (self.L + 1), se=0,
                verified=False, mse=self.mse,
                meverif=self.meverif, delta=self.delta,
            )

        if method == METHOD_DUALLATTICE:
            return self._run_lattice(cpp_gen)
        if method == METHOD_HARASE:
            return self._run_harase(cpp_gen)
        if method == METHOD_NOTPRIMITIVE:
            return self._run_notprimitive(cpp_gen)
        if method == METHOD_SIMD_NOTPRIMITIVE:
            return self._run_simd_notprimitive(cpp_gen)

        gen_L = cpp_gen.L()
        if self.L < gen_L:
            raise ValueError(
                f"TestME: EquidistributionTest.L ({self.L}) < generator L ({gen_L})"
            )

        return self._run_matricial(cpp_gen)

    def _resolve_method(self, gen_or_C) -> int:
        """
        Ask the C++ generator for its recommended equidistribution method
        when the YAML config did not specify one. Accepts the same shapes
        as `run()` (Generator, Combination, or raw `_cpp.Recurrence`).

        Primitives (e.g. MTGen) return their family default;
        `_cpp.CombinedF2LinearSource` short-circuits to the primitive's answer
        when it wraps a single component, and returns "notprimitive"
        unconditionally for J >= 2.

        Raises ValueError if the generator returns None.
        """
        cpp_gen = self._to_cpp_gen(gen_or_C)
        resolved = cpp_gen.default_test_method("equidistribution")
        if resolved is None:
            raise ValueError(
                "Generator has no default method for equidistribution; "
                "specify 'method:' in the YAML config"
            )
        method = _STR_TO_METHOD.get(resolved)
        if method is None:
            raise ValueError(
                f"Generator returned an unknown default method: {resolved!r}"
            )
        return method

    def _run_matricial(self, cpp_gen) -> EquidistributionResults:
        """Matricial orchestration loop in C++. Takes any `_cpp.Recurrence`."""
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)
        kg = cpp_gen.k()
        L = cpp_gen.L()

        result = _cpp.test_me_matricial(
            cpp_gen, kg, L, self.L, delta_capped, mse_capped,
        )
        psi12 = list(_cpp.compute_psi12(kg, self.L))
        return EquidistributionResults(
            L=self.L, ecart=list(result['ecart']), psi12=psi12,
            se=result['se'], verified=result['verified'], mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    # -- Private helpers --------------------------------------------------

    @staticmethod
    def _dimension_equid(mat, kg: int, l: int, L: int) -> int:
        return mat.dimension_equid(kg, l, L)

    @staticmethod
    def _resolution_equid(mat, kg: int, t: int, L: int, indices: list) -> int:
        return mat.resolution_equid(kg, t, L, indices)

    def _run_lattice(self, cpp_gen) -> EquidistributionResults:
        """TestMELat: dimension gaps via dual lattice basis reduction."""
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)
        kg = cpp_gen.k()
        L = cpp_gen.L()

        result = _cpp.test_me_lat(
            cpp_gen, kg, L, self.L, delta_capped, mse_capped,
        )
        psi12 = list(_cpp.compute_psi12(kg, self.L))
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    def _run_harase(self, cpp_gen) -> EquidistributionResults:
        """Harase-Matsumoto-Saito: primal lattice with Mulders-Storjohann."""
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)
        kg = cpp_gen.k()
        L = cpp_gen.L()

        result = _cpp.test_me_harase(
            cpp_gen, kg, L, self.L, delta_capped, mse_capped,
        )
        psi12 = list(_cpp.compute_psi12(kg, self.L))
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    def _run_notprimitive(self, cpp_gen) -> EquidistributionResults:
        """Matricial DE on the BM-recovered invariant subspace.

        Does NOT assume the combined characteristic polynomial is
        primitive (i.e. does not assume the generator is full-period).
        Reports ecart[v] = floor(p/v) - k(v) where p = dim of the
        certified invariant subspace V. See docs/C.md and notprimitive_de.h.
        """
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)
        kg = cpp_gen.k()
        L = cpp_gen.L()

        result = _cpp.test_me_notprimitive(
            cpp_gen, kg, L, self.L, delta_capped, mse_capped,
        )
        psi12 = list(_cpp.compute_psi12(kg, self.L))
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    def _run_simd_notprimitive(self, cpp_gen) -> EquidistributionResults:
        """SIMD-aware notprimitive — matches Saito-Matsumoto 2008 published
        k(v) for SFMT-like generators (dispatches to test_me_lat for
        non-SIMD full-period generators)."""
        import regpoly._regpoly_cpp as _cpp

        INT_MAX = 2**31 - 1
        delta_capped = [min(d, INT_MAX) for d in self.delta]
        mse_capped = min(self.mse, INT_MAX)
        kg = cpp_gen.k()
        L = cpp_gen.L()

        result = _cpp.test_me_notprimitive_simd(
            cpp_gen, kg, L, self.L, delta_capped, mse_capped,
        )
        psi12 = list(_cpp.compute_psi12(kg, self.L))
        return EquidistributionResults(
            L=self.L, ecart=result['ecart'], psi12=psi12,
            se=result['se'], verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )
