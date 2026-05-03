"""
tempering_optimizer.py — Tempering bitmask optimization.

Phase 2.4d: the recursive optimize(v) inner loop now lives in C++
(see packages/regpoly-cpp/src/search/tempering_optimizer.cpp). Python
still owns:

  * Safe-mask construction (small structural computation that depends
    on `mu` and bit width — runs once per Combination).
  * Building the (cpp_trans, param_name, width) locator list that the
    C++ driver consumes.
  * Syncing the best-found values back into each Transformation's
    Python `_params` dict after the C++ driver returns.

Modes (controlled by delta, mse, n_restarts):

  delta=None, mse=None, n_restarts=1 (default):
      ME mode.  delta=[0]*L, mse=0.  Single call.

  delta=[...], mse=N, n_restarts=1:
      Find b,c with gaps[v] <= delta[v] and se <= mse.

  n_restarts > 1 (delta and mse ignored):
      Minimize se via iterative delta tightening.
"""

from __future__ import annotations

import sys
import time

import regpoly._regpoly_cpp as _cpp
from regpoly.core.combination import Combination

INT_MAX = 2**31 - 1


class TemperingOptimizer:

    def __init__(
        self,
        max_essais: int = 400,
        delta: list[int] | None = None,
        mse: int | None = None,
        n_restarts: int = 1,
        verbose: bool = True,
    ) -> None:
        self.max_essais = max_essais
        self.delta = delta
        self.mse = mse
        self.n_restarts = n_restarts
        self.verbose = verbose

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, comb: Combination) -> OptResult:
        """
        Optimize the bitmask parameters of the transformations in comb.

        Builds all comb-specific state (safe masks, C++ caches), runs
        the optimization, and returns an OptResult.
        """
        ctx = _RunContext(self, comb)

        if self.n_restarts > 1:
            return ctx.run_minimize()

        delta = self.delta
        mse = self.mse
        if delta is None:
            delta = [0] * (ctx.L + 1)
        if mse is None:
            mse = 0

        return ctx.run_once(delta, mse, self.verbose)


class _RunContext:
    """Per-run state tied to a specific Combination."""

    def __init__(self, opt: TemperingOptimizer, comb: Combination) -> None:
        self.opt = opt
        self.comb = comb
        self.L = comb.L

        # Stable, ordered list of optimizable params:
        #   self._locators = [(j, ti, pn, width, trans), ...]
        # — order matters because the C++ driver indexes safe_masks
        # by position (not by (j, ti, pn) tuple).
        self._locators: list[tuple] = []
        for j, comp in enumerate(comb.components):
            for ti, trans in enumerate(comp.trans):
                for pn, width in trans.optimizable_params():
                    self._locators.append((j, ti, pn, width, trans))

        if not self._locators:
            raise ValueError("No optimizable bitmask parameters found")

        self._safe_masks_per_locator = self._compute_safe_masks()
        self._gens = [comb[j]._cpp_gen for j in range(comb.J)]
        self._trans_cpp = [[t._cpp_trans for t in comp.trans]
                           for comp in comb.components]

    # ------------------------------------------------------------------
    # Safe mask computation (analytical, mirrors temperMK.c)
    # ------------------------------------------------------------------

    def _compute_safe_masks(self) -> list[list[int]]:
        """Compute safe masks shaped [L+1][P] for the C++ driver.

        At resolution v, locator i's safe mask is a bitmap of which
        bits of its parameter may be perturbed without violating
        equidistribution invariants from prior resolutions. mu-aware
        masking (where applicable) excludes bit positions that would
        flip output bits already in the perturbed range.
        """
        L = self.L
        P = len(self._locators)
        safe_masks: list[list[int]] = [
            [0] * P for _ in range(L + 1)
        ]

        # For each locator, look up its mu (only meaningful for the `b`
        # parameter of tempMK; absent on other transformations).
        mu_per_locator: list[int] = []
        for (_j, _ti, pn, _w, trans) in self._locators:
            mu = 0
            if pn == 'b':
                try:
                    mu = trans.get_param('mu')
                except (KeyError, AttributeError):
                    mu = 0
            mu_per_locator.append(mu)

        for v in range(1, L + 1):
            for i, (_j, _ti, _pn, w, _trans) in enumerate(self._locators):
                nbits = w - v + 1
                if nbits <= 0:
                    safe_masks[v][i] = 0
                    continue

                mask = (1 << nbits) - 1
                mu = mu_per_locator[i]
                if mu > 0:
                    for p in range(mu, min(v + mu - 1, w)):
                        py_bit = w - 1 - p
                        if 0 <= py_bit < nbits:
                            mask &= ~(1 << py_bit)

                safe_masks[v][i] = mask

        if self.opt.verbose:
            for v in [1, L // 2, L]:
                total = sum(bin(m).count('1') for m in safe_masks[v])
                print(f"  Safe mask v={v}: {total} bits")

        return safe_masks

    # ------------------------------------------------------------------
    # Single recursive optimization pass
    # ------------------------------------------------------------------

    def run_once(self, delta: list[int], mse: int,
                 verbose: bool) -> OptResult:
        L = self.L
        kg = self.comb.k_g

        if verbose:
            print(f"Optimizing {len(self._locators)} params, "
                  f"k_g={kg}, L={L}, mse={mse}")
            sys.stdout.flush()

        cache = _cpp.TemperOptCache(self._gens, self._trans_cpp, kg, L)

        cfg = _cpp.TemperingOptimizerConfig()
        cfg.max_essais = self.opt.max_essais
        cfg.delta = list(delta)
        cfg.mse = mse
        cfg.n_restarts = 1
        cfg.random_seed = 0

        result, final_values = _cpp.run_tempering_optimizer_once(
            cfg, cache, self._locator_tuples(), self._safe_masks_per_locator)

        self._sync_python_params_from_locators(final_values)

        gaps_dict = {v: result.gaps[v] for v in range(1, L + 1)}
        if verbose:
            print(f"\nOptimization complete: se={result.se}, "
                  f"{result.elapsed_seconds:.1f}s")

        return OptResult(se=result.se, gaps=gaps_dict,
                         elapsed=result.elapsed_seconds,
                         essais=result.essais)

    # ------------------------------------------------------------------
    # Iterative delta tightening (minimize se)
    # ------------------------------------------------------------------

    def run_minimize(self) -> OptResult:
        L = self.L
        n_restarts = self.opt.n_restarts
        t_start = time.time()

        if self.opt.verbose:
            print(f"Minimize se: k_g={self.comb.k_g}, L={L}, "
                  f"n_restarts={n_restarts}")
            sys.stdout.flush()

        cache = _cpp.TemperOptCache(self._gens, self._trans_cpp,
                                      self.comb.k_g, L)

        cfg = _cpp.TemperingOptimizerConfig()
        cfg.max_essais = self.opt.max_essais
        cfg.delta = [INT_MAX] * (L + 1)
        cfg.mse = INT_MAX
        cfg.n_restarts = n_restarts
        cfg.random_seed = 0

        result, final_values = _cpp.run_tempering_optimizer_minimize(
            cfg, cache, self._locator_tuples(), self._safe_masks_per_locator)

        self._sync_python_params_from_locators(final_values)

        gaps_dict = {v: result.gaps[v] for v in range(1, L + 1)}
        elapsed = time.time() - t_start
        if self.opt.verbose:
            print(f"\nMinimize complete: se={result.se}, "
                  f"{result.essais} calls, {elapsed:.1f}s")

        return OptResult(se=result.se, gaps=gaps_dict,
                         elapsed=elapsed, essais=result.essais)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _compute_gaps(self) -> tuple[dict[int, int], int]:
        cache = _cpp.TemperOptCache(
            self._gens, self._trans_cpp, self.comb.k_g, self.L)
        ecart = cache.compute_all()
        gaps = {v: ecart[v] for v in range(1, self.L + 1)}
        se = sum(gaps.values())
        return gaps, se

    def _locator_tuples(self) -> list[tuple]:
        """Build the (cpp_trans, param_name, width, current_value) list
        the C++ driver consumes."""
        return [
            (trans._cpp_trans, pn, w, trans.get_param(pn))
            for (_j, _ti, pn, w, trans) in self._locators
        ]

    def _sync_python_params_from_locators(self,
                                           final_values: list[int]) -> None:
        """After the C++ driver returns, write the best-found values
        back into each Transformation's Python `_params` dict so
        downstream `get_param`/`set_param` calls see consistent state.
        The underlying _cpp_trans has already been updated."""
        for (_j, _ti, pn, _w, trans), val in zip(
                self._locators, final_values):
            trans._params[pn] = val


class OptResult:
    __slots__ = ("se", "gaps", "elapsed", "essais")

    def __init__(self, se: int, gaps: dict[int, int],
                 elapsed: float, essais: int) -> None:
        self.se = se
        self.gaps = gaps
        self.elapsed = elapsed
        self.essais = essais
