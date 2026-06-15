# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
tempering_search.py — Search for optimal transformation parameters.

For each combination of generators, randomly samples transformation
parameters and evaluates a test.  When the test is equidistribution
and the transformations have optimizable bitmask parameters, the
TemperingOptimizer is invoked to refine the bitmasks at each try.

Phase 2.4c: the per-combo / per-try loop now lives in C++
(`run_tempering_search`).  This Python wrapper builds a C++
Combination from the Python one, registers callbacks for the per-try
work (which still touches Python Transformation.randomize_params, the
Python optimizer wrapper, and the Python test object), and assembles
the final list of TemperingSearchResult instances for callers.

Usage::

    from regpoly.search.tempering_optimizer import TemperingOptimizer

    optimizer = TemperingOptimizer(comb, max_essais=400,
                                   delta=delta, mse=mse)

    search = TemperingSearch(
        gen_lists=[[gen1, gen2], [gen3]],
        temperings=[[trans_mk], [trans_lag]],
        test=equidist_test,
        Lmax=32,
        nb_tries=100,
        optimizer=optimizer,
    )
    search.run()
"""

from __future__ import annotations

import sys
import time
from typing import TYPE_CHECKING

from regpoly.analyses.abstract_test import AbstractTest, AbstractTestResults
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.helpers import make_combined
from regpoly.io.combo_builder import build_cpp_enumerator
from regpoly.io.tested_generator import save_tested_generator
from regpoly.search.seek import _current_combo_snapshot
from regpoly_cpp import _regpoly_cpp as _cpp

if TYPE_CHECKING:
    from regpoly.search.tempering_optimizer import TemperingOptimizer


class TemperingSearchResult:
    """Result of one generator combination."""
    __slots__ = ("se", "test_result", "params", "elapsed")

    def __init__(self, se: int, test_result: AbstractTestResults,
                 params: list[list[dict]], elapsed: float) -> None:
        self.se = se
        self.test_result = test_result
        self.params = params
        self.elapsed = elapsed


class TemperingSearch:
    """
    Search for optimal transformation parameters over combined generators.

    Parameters
    ----------
    gen_lists : list[list[Generator]]
        Per-component list of candidate generators (fully defined).
        When two entries are the same object, the components share a
        pool and C(n,k) selection is enforced.
    temperings : list[list[Transformation]]
        Per-component tempering chain.  Non-structural params with a
        rand_type are re-randomized at each try.
    test : AbstractTest
        The test to evaluate (e.g. EquidistributionTest).
    Lmax : int
        Output resolution cap.
    nb_tries : int
        Number of random parameter sets to try per combination.
    optimizer : TemperingOptimizer or None
        A pre-configured TemperingOptimizer whose settings (max_essais,
        delta, mse, n_restarts, verbose) are reused for each combination.
        If None and the test is equidistribution with optimizable params,
        a default optimizer is created (ME mode, max_essais=400).
    output_dir : str or None
        If set, save best results to YAML files.
    verbose : bool
        Print progress.
    """

    def __init__(
        self,
        gen_lists: list,
        temperings: list,
        test: AbstractTest,
        Lmax: int,
        nb_tries: int = 100,
        optimizer: "TemperingOptimizer | None" = None,
        output_dir: str | None = None,
        verbose: bool = True,
    ) -> None:
        self.gen_lists = gen_lists
        self.temperings = temperings
        self.test = test
        self.Lmax = Lmax
        self.nb_tries = nb_tries
        self.optimizer = optimizer
        self.output_dir = output_dir
        self.verbose = verbose

    def run(self) -> list[TemperingSearchResult]:
        """
        Run the search over all generator combinations.

        Returns a list of TemperingSearchResult, one per combination
        that produced a valid result.
        """
        # Build the C++ enumerator directly from the per-slot pools.
        # No Python Combination middleman.
        try:
            cpp_comb = build_cpp_enumerator(
                self.gen_lists, self.temperings, self.Lmax)
        except RuntimeError:
            if self.verbose:
                print("No valid generator combination found.")
            return []

        gen_pools = self.gen_lists
        tempering_pools = self.temperings

        use_optimizer = (
            isinstance(self.test, EquidistributionTest)
            and self._has_optimizable(tempering_pools)
        )

        if use_optimizer and self.optimizer is None:
            from regpoly.search.tempering_optimizer import TemperingOptimizer
            self.optimizer = TemperingOptimizer(verbose=False)

        results: list[TemperingSearchResult] = []
        # Per-combo state shared by callbacks.  best_params is rebuilt
        # eagerly whenever a new best score is observed in on_try, so
        # the snapshot is stable when on_combo_done fires.
        state = {
            "best_se": None,
            "best_result": None,
            "best_params": None,
            "t_combo_start": 0.0,
            "is_first_combo": True,
        }

        def _active_gens():
            """Active Python Generator per slot, looked up via current_gen."""
            return [
                gen_pools[j][cpp_comb.pool(j).current_index()]
                for j in range(cpp_comb.J)
            ]

        def on_combo_start(_cpp_comb, combo_idx):
            # The C++ Combination has already advanced before this fires;
            # no Python-side mirror is needed.
            state["is_first_combo"] = False
            state["best_se"] = None
            state["best_result"] = None
            state["best_params"] = None
            state["t_combo_start"] = time.time()
            if self.verbose:
                active = _active_gens()
                gens_str = " + ".join(
                    f"{g.name()}(k={g.k})" for g in active)
                print(f"\nCombination {combo_idx}: {gens_str}, "
                      f"k_g={cpp_comb.k_g}, L={cpp_comb.L}")
                sys.stdout.flush()

        def on_try(_cpp_comb, _combo_idx, try_idx, _is_first):
            for chain in tempering_pools:
                for trans in chain:
                    trans.randomize_params()

            if use_optimizer:
                try:
                    # Optimizer mutates tempering_pools' Transformations in place.
                    self.optimizer.run(_current_combo_snapshot(
                        cpp_comb, gen_pools, tempering_pools))
                except ValueError:
                    pass  # no optimizable params for this combo

            # Rebuild a fresh CombinedF2LinearSource from the current pool +
            # randomized tempering. Same allocation count as before
            # (the previous _run_* helpers built CombinedF2LinearSource per
            # call too); a future _cpp.CombinedF2LinearSource.rebind_tempering
            # would amortize this further.
            combined = make_combined(
                *_active_gens(), trans=tempering_pools, Lmax=cpp_comb.Lmax)
            test_result = self.test.run(combined)
            score = self._score(test_result)

            outcome = _cpp.TemperingTryResult()
            outcome.got_result = True
            outcome.score = score

            if state["best_se"] is None or score < state["best_se"]:
                state["best_se"] = score
                state["best_result"] = test_result
                state["best_params"] = self._save_params(tempering_pools)
                if self.verbose:
                    print(f"  try {try_idx + 1:>4d}: se={score}"
                          f"{'  ME!' if score == 0 else ''}")
                    sys.stdout.flush()
            return outcome

        def on_combo_done(_cpp_comb, _combo_idx, _best_score, _best_try):
            if state["best_params"] is None:
                return
            self._restore_params(tempering_pools, state["best_params"])
            elapsed = time.time() - state["t_combo_start"]
            if self.verbose:
                print(f"  Best: se={state['best_se']} ({elapsed:.1f}s)")
            if self.output_dir and state["best_result"] is not None:
                test_results = self._build_results_dict(state["best_result"])
                snapshot = _current_combo_snapshot(
                    cpp_comb, gen_pools, tempering_pools)
                path = save_tested_generator(
                    self.output_dir, "equidist", snapshot, test_results)
                if self.verbose:
                    print(f"  Saved: {path}")
            results.append(TemperingSearchResult(
                se=state["best_se"],
                test_result=state["best_result"],
                params=state["best_params"],
                elapsed=elapsed))

        cfg = _cpp.TemperingSearchConfig()
        cfg.nb_tries = self.nb_tries
        cfg.progress_interval = 0

        cpp_result = _cpp.run_tempering_search(
            cpp_comb, cfg,
            on_combo_start=on_combo_start,
            on_try=on_try,
            on_combo_done=on_combo_done,
            on_progress=None)

        if self.verbose:
            print(f"\nSearch complete: {cpp_result.nbgen} combinations, "
                  f"{cpp_result.nb_with_result} with results, "
                  f"{cpp_result.elapsed_seconds:.1f}s")

        return results

    # -- Helpers -----------------------------------------------------------

    @staticmethod
    def _has_optimizable(tempering_pools: list[list]) -> bool:
        """Check if any transformation in the pools has optimizable params."""
        for chain in tempering_pools:
            for trans in chain:
                if trans.optimizable_params():
                    return True
        return False

    @staticmethod
    def _score(result: AbstractTestResults) -> int:
        """Extract a numeric score from a test result (lower is better)."""
        if hasattr(result, 'se'):
            return result.se
        if hasattr(result, 'verified') and result.verified:
            return 0
        return 1

    @staticmethod
    def _save_params(tempering_pools: list[list]) -> list[list[dict]]:
        """Save all transformation params for all components."""
        return [[trans.params for trans in chain]
                for chain in tempering_pools]

    @staticmethod
    def _restore_params(tempering_pools: list[list],
                        saved: list[list[dict]]) -> None:
        """Restore transformation params from a saved snapshot."""
        for chain, comp_params in zip(tempering_pools, saved):
            for trans, params in zip(chain, comp_params):
                for name, value in params.items():
                    if value != trans.get_param(name):
                        trans.set_param(name, value)

    @staticmethod
    def _build_results_dict(result: AbstractTestResults) -> dict:
        """Build a results dict for save_tested_generator."""
        results = {}
        if hasattr(result, 'se'):
            eq = {"se": result.se}
            if hasattr(result, 'is_me') and result.is_me():
                eq["status"] = "ME"
            if hasattr(result, 'ecart') and hasattr(result, 'L'):
                ecart = {l: result.ecart[l]
                         for l in range(1, result.L + 1)
                         if result.ecart[l] != 0}
                if ecart:
                    eq["ecart"] = ecart
            results["equidistribution"] = eq
        return results
