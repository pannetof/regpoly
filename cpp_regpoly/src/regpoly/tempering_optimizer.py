"""
tempering_optimizer.py — Tempering bitmask optimization.

Mirrors the C code's OptimizeTemper() from temperMK.c:

  - Safe masks: at resolution v, only perturb bits at positions >= v-1.
  - Random multi-bit perturbation within safe mask.
  - Incremental StackBase via step(v).
  - ONE recursive algorithm controlled by delta[v] and mse.

run(delta, mse, n_restarts):

  delta=None, mse=None, n_restarts=1:
      ME mode.  delta=[0]*L, mse=0.  Single call.

  delta=[...], mse=N, n_restarts=1:
      Find b,c with gaps[v] <= delta[v] and se <= mse.

  n_restarts > 1 (delta and mse ignored):
      Minimize se via iterative delta tightening.
      First pass: delta=INT_MAX (one deep sweep).
      Subsequent passes: delta = best gaps so far, mse = best se so far.
      Stop after n_restarts consecutive failures.
"""

from __future__ import annotations

import random
import sys
import time

import regpoly._regpoly_cpp as _cpp
from regpoly.combinaison import Combinaison

INT_MAX = 2**31 - 1


class TemperingOptimizer:

    def __init__(
        self,
        comb: Combinaison,
        max_essais: int = 400,
        verbose: bool = True,
    ) -> None:
        self.comb = comb
        self.max_essais = max_essais
        self.verbose = verbose
        self.L = comb.L

        # (j, ti, param_name) -> bit width
        self._param_widths: dict[tuple[int, int, str], int] = {}
        for j, comp in enumerate(comb.components):
            for ti, trans in enumerate(comp.trans):
                for pn, width in trans.optimizable_params():
                    self._param_widths[(j, ti, pn)] = width

        if not self._param_widths:
            raise ValueError("No optimizable bitmask parameters found")

        self._safe_masks = self._compute_safe_masks()
        self._gens = [comb[j]._cpp_gen for j in range(comb.J)]
        self._trans_cpp = [[t._cpp_trans for t in comp.trans]
                           for comp in comb.components]

    # ------------------------------------------------------------------
    # Safe mask computation (analytical, mirrors temperMK.c)
    # ------------------------------------------------------------------

    def _compute_safe_masks(self) -> list[dict[tuple, int]]:
        L = self.L
        safe_masks: list[dict[tuple, int]] = [{} for _ in range(L + 1)]

        # Detect mu shift for TemperMK mask_b adjustment
        param_mu: dict[tuple, int] = {}
        for (j, ti, pn), w in self._param_widths.items():
            if pn == 'b':
                trans = self.comb.components[j].trans[ti]
                try:
                    param_mu[(j, ti, pn)] = trans.get_param('mu')
                except (KeyError, AttributeError):
                    pass

        for v in range(1, L + 1):
            for key, w in self._param_widths.items():
                nbits = w - v + 1
                if nbits <= 0:
                    safe_masks[v][key] = 0
                    continue

                mask = (1 << nbits) - 1

                # TemperMK: clear bits whose b-change propagates
                # through the c-step to positions below v-1
                mu = param_mu.get(key, 0)
                if mu > 0:
                    for p in range(mu, min(v + mu - 1, w)):
                        py_bit = w - 1 - p
                        if 0 <= py_bit < nbits:
                            mask &= ~(1 << py_bit)

                safe_masks[v][key] = mask

        if self.verbose:
            for v in [1, L // 2, L]:
                total = sum(bin(m).count('1') for m in safe_masks[v].values())
                print(f"  Safe mask v={v}: {total} bits")

        return safe_masks

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, delta: list[int] | None = None,
            mse: int | None = None,
            n_restarts: int = 1) -> OptResult:

        if n_restarts > 1:
            return self._run_minimize(n_restarts)

        if delta is None:
            delta = [0] * (self.L + 1)
        if mse is None:
            mse = 0

        return self._run_once(delta, mse, self.verbose)

    # ------------------------------------------------------------------
    # Single recursive optimization pass
    # ------------------------------------------------------------------

    def _run_once(self, delta: list[int], mse: int,
                  verbose: bool) -> OptResult:
        comb = self.comb
        L = self.L
        kg = comb.k_g
        t_start = time.time()
        safe_masks = self._safe_masks

        if verbose:
            print(f"Optimizing {len(self._param_widths)} params, "
                  f"k_g={kg}, L={L}, mse={mse}")
            sys.stdout.flush()

        cache = _cpp.LatticeOptCache(self._gens, self._trans_cpp, kg, L)
        cache.reset_step()

        ecart = [-1] * (L + 1)
        best_se = [INT_MAX] * (L + 1)
        best_params = self._save_params()
        essais = 0
        max_v = 0

        def optimize(v):
            nonlocal essais, max_v, best_params

            Lim = 5 if v < L // 2 else 2

            for _ in range(Lim):
                if essais >= self.max_essais:
                    break
                essais += 1

                # Random multi-bit perturbation within safe mask
                for (j, ti, pn), w in self._param_widths.items():
                    mask = safe_masks[v].get((j, ti, pn), 0)
                    if mask == 0:
                        continue
                    trans = comb.components[j].trans[ti]
                    perturbation = random.getrandbits(w) & mask
                    trans.set_param(pn, trans.get_param(pn) ^ perturbation)

                ecart[v] = cache.step(v)
                se = sum(ecart[1:v + 1])

                # Save best if improvement at this depth or deeper
                if se < best_se[v] and v >= max_v:
                    essais = 0
                    best_params = self._save_params()
                    best_se[v] = se
                    max_v = v

                    if verbose:
                        print(f"  v={v:>3d}: se={se} (max_v={max_v})")
                        sys.stdout.flush()

                # Recurse if gap within tolerance
                if ecart[v] <= delta[v] and se <= mse:
                    if v >= L:
                        break
                    else:
                        optimize(v + 1)

                # Early exit if target met at deepest level
                if (ecart[L] >= 0
                        and ecart[L] <= delta[L]
                        and best_se[L] <= mse):
                    break

        optimize(1)
        self._restore_params(best_params)

        # Final verification
        gaps, se = self._compute_gaps()

        elapsed = time.time() - t_start
        if verbose:
            print(f"\nOptimization complete: se={se}, "
                  f"max_v={max_v}, {elapsed:.1f}s")

        return OptResult(se=se, gaps=gaps, elapsed=elapsed,
                         essais=essais)

    # ------------------------------------------------------------------
    # Iterative delta tightening (minimize se)
    # ------------------------------------------------------------------

    def _run_minimize(self, n_restarts: int) -> OptResult:
        L = self.L
        t_start = time.time()

        if self.verbose:
            print(f"Minimize se: k_g={self.comb.k_g}, L={L}, "
                  f"n_restarts={n_restarts}")
            sys.stdout.flush()

        best_result = None
        best_params = None
        failures = 0
        n_calls = 0

        while failures < n_restarts:
            n_calls += 1

            # Randomize starting point
            for (j, ti, pn), w in self._param_widths.items():
                self.comb.components[j].trans[ti].set_param(
                    pn, random.getrandbits(w))

            if best_result is None:
                delta = [INT_MAX] * (L + 1)
                mse_val = INT_MAX
            else:
                delta = [best_result.gaps.get(v, INT_MAX)
                         for v in range(L + 1)]
                mse_val = best_result.se

            result = self._run_once(delta, mse_val, verbose=False)

            if best_result is None or result.se < best_result.se:
                best_result = result
                best_params = self._save_params()
                failures = 0
                if self.verbose:
                    print(f"  call {n_calls:>3d}: se={result.se} "
                          f"({result.elapsed:.1f}s)")
                    sys.stdout.flush()
            else:
                failures += 1

        if best_params is not None:
            self._restore_params(best_params)

        elapsed = time.time() - t_start
        if self.verbose:
            print(f"\nMinimize complete: se={best_result.se}, "
                  f"{n_calls} calls, {elapsed:.1f}s")

        return OptResult(se=best_result.se, gaps=best_result.gaps,
                         elapsed=elapsed, essais=n_calls)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _compute_gaps(self) -> tuple[dict[int, int], int]:
        cache = _cpp.LatticeOptCache(
            self._gens, self._trans_cpp, self.comb.k_g, self.L)
        ecart = cache.compute_all()
        gaps = {v: ecart[v] for v in range(1, self.L + 1)}
        se = sum(gaps.values())
        return gaps, se

    def _save_params(self) -> dict:
        saved = {}
        for j, comp in enumerate(self.comb.components):
            for ti, trans in enumerate(comp.trans):
                for pn, _ in trans.optimizable_params():
                    saved[(j, ti, pn)] = trans.get_param(pn)
        return saved

    def _restore_params(self, saved: dict) -> None:
        for (j, ti, pn), val in saved.items():
            self.comb.components[j].trans[ti].set_param(pn, val)


class OptResult:
    __slots__ = ("se", "gaps", "elapsed", "essais")

    def __init__(self, se: int, gaps: dict[int, int],
                 elapsed: float, essais: int) -> None:
        self.se = se
        self.gaps = gaps
        self.elapsed = elapsed
        self.essais = essais
