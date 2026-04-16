"""
tempering_search.py — Search for optimal transformation parameters.

For each combination of generators, randomly samples transformation
parameters and evaluates a test.  When the test is equidistribution
and the transformations have optimizable bitmask parameters, the
TemperingOptimizer is invoked to refine the bitmasks at each try.

Usage::

    from regpoly.tempering_optimizer import TemperingOptimizer

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

from regpoly.combinaison import Combinaison
from regpoly.analyses.test_base import AbstractTest
from regpoly.analyses.test_results_base import AbstractTestResults
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.tested_generator import save_tested_generator


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
    gen_lists : list[list[Generateur]]
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
        comb = Combinaison.CreateFromFiles(
            self.gen_lists, self.Lmax, self.temperings)

        if not comb.reset():
            if self.verbose:
                print("No valid generator combination found.")
            return []

        use_optimizer = (
            isinstance(self.test, EquidistributionTest)
            and self._has_optimizable(comb)
        )

        if use_optimizer and self.optimizer is None:
            from regpoly.tempering_optimizer import TemperingOptimizer
            self.optimizer = TemperingOptimizer(verbose=False)

        results = []
        t_total = time.time()
        combo_idx = 0

        while True:
            combo_idx += 1
            result = self._search_one_combo(comb, use_optimizer, combo_idx)
            if result is not None:
                results.append(result)

            try:
                next(comb)
            except StopIteration:
                break

        elapsed_total = time.time() - t_total
        if self.verbose:
            print(f"\nSearch complete: {combo_idx} combinations, "
                  f"{len(results)} with results, {elapsed_total:.1f}s")

        return results

    def _search_one_combo(self, comb: Combinaison,
                          use_optimizer: bool,
                          combo_idx: int) -> TemperingSearchResult | None:
        """Search over nb_tries for one generator combination."""
        t_start = time.time()

        if self.verbose:
            gens_str = " + ".join(
                f"{comb[j].name()}(k={comb[j].k})"
                for j in range(comb.J))
            print(f"\nCombination {combo_idx}: {gens_str}, "
                  f"k_g={comb.k_g}, L={comb.L}")
            sys.stdout.flush()

        best_se = None
        best_result = None
        best_params = None

        for t in range(self.nb_tries):
            # Re-randomize non-structural transformation params in place
            for comp in comb.components:
                for trans in comp.trans:
                    trans.randomize_params()

            # Optimize bitmask params if applicable
            if use_optimizer:
                self._run_optimizer(comb)

            # Run the test
            test_result = self.test.run(comb)

            # Extract score (se for equidistribution, generic for others)
            se = self._score(test_result)

            if best_se is None or se < best_se:
                best_se = se
                best_result = test_result
                best_params = self._save_params(comb)

                if self.verbose:
                    print(f"  try {t + 1:>4d}: se={se}"
                          f"{'  ME!' if se == 0 else ''}")
                    sys.stdout.flush()

                if se == 0:
                    break

        if best_params is None:
            return None

        # Restore best params
        self._restore_params(comb, best_params)

        elapsed = time.time() - t_start

        if self.verbose:
            print(f"  Best: se={best_se} ({elapsed:.1f}s)")

        # Save if output_dir specified
        if self.output_dir and best_result is not None:
            test_results = self._build_results_dict(best_result)
            path = save_tested_generator(
                self.output_dir, "equidist", comb, test_results)
            if self.verbose:
                print(f"  Saved: {path}")

        return TemperingSearchResult(
            se=best_se, test_result=best_result,
            params=best_params, elapsed=elapsed)

    def _run_optimizer(self, comb: Combinaison) -> None:
        """Run the optimizer on this comb."""
        try:
            self.optimizer.run(comb)
        except ValueError:
            pass  # no optimizable params for this combo

    # -- Helpers -----------------------------------------------------------

    @staticmethod
    def _has_optimizable(comb: Combinaison) -> bool:
        """Check if any transformation in the combinaison has optimizable params."""
        for comp in comb.components:
            for trans in comp.trans:
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
    def _save_params(comb: Combinaison) -> list[list[dict]]:
        """Save all transformation params for all components."""
        return [[trans.params for trans in comp.trans]
                for comp in comb.components]

    @staticmethod
    def _restore_params(comb: Combinaison,
                        saved: list[list[dict]]) -> None:
        """Restore transformation params from a saved snapshot."""
        for comp, comp_params in zip(comb.components, saved):
            for trans, params in zip(comp.trans, comp_params):
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
