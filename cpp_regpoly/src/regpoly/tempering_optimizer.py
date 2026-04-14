"""
tempering_optimizer.py — Backtracking optimization of tempering bitmask
parameters to minimize equidistribution dimension gaps.

Based on Harase (2009), generalized to any transformation type.

Two phases:
  1. Candidate filtering — for each resolution v, determine which bits
     of which bitmask parameters can affect the top v output bits.
  2. Backtracking optimization — process resolutions in order, trying
     only candidate bits, with backtracking when stuck.

Supports two equidistribution methods:
  - "matricial": Gaussian elimination (fast for small k, O(k^2 L))
  - "lattice":   dual lattice basis reduction (required for large k)
"""

from __future__ import annotations

import sys
import time

import regpoly._regpoly_cpp as _cpp
from regpoly.bitvect import BitVect
from regpoly.combinaison import Combinaison

METHOD_MATRICIAL = 0
METHOD_LATTICE = 1


class TemperingOptimizer:
    """
    Optimize tempering bitmask parameters to minimize the total
    dimension defect (Delta).

    The transformations in *comb* are modified in place.

    Parameters
    ----------
    comb : Combinaison
        Generator + tempering, already initialized (reset).
    method : str
        "matricial" (default for k < 2000) or "lattice" (for large k).
    max_backtracks : int
        Maximum backtrack attempts before giving up.
    verbose : bool
        Print progress to stdout.

    Usage::

        opt = TemperingOptimizer(comb)
        result = opt.run()
        # result.se == 0 means ME
    """

    def __init__(
        self,
        comb: Combinaison,
        method: str | None = None,
        max_backtracks: int = 500,
        verbose: bool = True,
    ) -> None:
        self.comb = comb
        self.max_backtracks = max_backtracks
        self.verbose = verbose
        self.L = comb.L

        if method is None:
            method = "lattice" if comb.k_g > 2000 else "matricial"
        self._method = METHOD_LATTICE if method == "lattice" else METHOD_MATRICIAL

        # Collect all flippable bits
        self._bits: list[tuple[int, int, str, int]] = []
        for j, comp in enumerate(comb.components):
            for ti, trans in enumerate(comp.trans):
                for param_name, width in trans.optimizable_params():
                    for bit in range(width):
                        self._bits.append((j, ti, param_name, bit))

        if not self._bits:
            raise ValueError("No optimizable bitmask parameters found")

    # ══════════════════════════════════════════════════════════════════════
    # Phase 1: Candidate filtering
    # ══════════════════════════════════════════════════════════════════════

    def _compute_candidates(self) -> dict[int, list[int]]:
        """
        For each resolution v = 1..L, return indices into self._bits
        that can affect the top v output bits.
        """
        comb = self.comb
        L = self.L

        # Build a test state by running the generator one step from a
        # random init
        gen = comb[0]
        k = gen.k
        import random
        gen_copy = gen.copy()
        init_bv = BitVect(k)
        rng = random.getrandbits(k)
        for i in range(k):
            init_bv.put_bit(i, (rng >> i) & 1)
        gen_copy.initialize_state(init_bv)
        next(gen_copy)

        # Get the full output state
        ref_output = gen_copy._cpp_gen.get_output()

        # Build the transformation chain
        chain = []
        for comp in comb.components:
            chain.extend(comp.trans)

        def apply_chain(state_bv):
            s = state_bv.copy()
            for t in chain:
                t._cpp_trans.apply(s)
            return s

        ref_result = apply_chain(ref_output)
        ref_nbits = ref_result.nbits()

        candidates: dict[int, list[int]] = {}
        w = L

        for bi, (j, ti, param_name, bit_pos) in enumerate(self._bits):
            trans = comb.components[j].trans[ti]

            trans.flip_bit(param_name, bit_pos)
            flipped_result = apply_chain(ref_output)
            trans.flip_bit(param_name, bit_pos)

            # Find which of the top w output bits differ
            diff_positions = []
            for b in range(min(w, ref_nbits)):
                if ref_result.get_bit(b) != flipped_result.get_bit(b):
                    diff_positions.append(b)

            if diff_positions:
                min_diff = min(diff_positions)
                for v in range(min_diff + 1, L + 1):
                    candidates.setdefault(v, []).append(bi)

        if self.verbose:
            total = sum(len(c) for c in candidates.values())
            naive = len(self._bits) * L
            pct = 100 * total / max(naive, 1)
            print(f"  Candidate filtering: {total} candidates "
                  f"(vs {naive} naive, {pct:.0f}%)")

        return candidates

    # ══════════════════════════════════════════════════════════════════════
    # Phase 2: Backtracking optimization
    # ══════════════════════════════════════════════════════════════════════

    def run(self) -> OptResult:
        comb = self.comb
        L = self.L
        kg = comb.k_g
        t_start = time.time()
        method_name = "lattice" if self._method == METHOD_LATTICE else "matricial"

        if self.verbose:
            n_trans = sum(len(comp.trans) for comp in comb.components)
            print(f"Optimizing {len(self._bits)} bits across "
                  f"{n_trans} transformation(s)")
            print(f"  k_g={kg}, L={L}, method={method_name}")
            sys.stdout.flush()

        # Phase 1
        candidates = self._compute_candidates()

        # Initial gaps
        gaps = self._compute_all_gaps()
        se = sum(gaps[v] for v in range(1, L + 1))

        if self.verbose:
            print(f"  Initial se = {se}")
            sys.stdout.flush()

        if se == 0:
            return OptResult(se=0, gaps=gaps,
                             elapsed=time.time() - t_start, backtracks=0)

        # Backtracking
        stack: list[tuple[int, int, dict]] = []
        backtracks = 0
        v = 1
        c_start = 0

        while v <= L:
            if gaps.get(v, 0) == 0:
                v += 1
                c_start = 0
                continue

            cands = candidates.get(v, [])
            found_fix = False

            for ci in range(c_start, len(cands)):
                bi = cands[ci]
                j, ti, param_name, bit_pos = self._bits[bi]
                trans = comb.components[j].trans[ti]

                saved = self._save_params()
                trans.flip_bit(param_name, bit_pos)

                new_gaps = self._compute_gaps_up_to(v)
                all_zero = all(new_gaps.get(vp, 0) == 0
                               for vp in range(1, v + 1))

                if all_zero:
                    stack.append((v, ci, saved))
                    for vp in range(1, v + 1):
                        gaps[vp] = new_gaps[vp]
                    found_fix = True

                    if self.verbose:
                        print(f"  v={v:>3d}: flip {param_name}[{bit_pos:>2d}] "
                              f"-> gaps 1..{v} all zero")
                        sys.stdout.flush()

                    v = 1
                    c_start = 0
                    break
                else:
                    self._restore_params(saved)

            if not found_fix:
                if not stack or backtracks >= self.max_backtracks:
                    if self.verbose:
                        reason = ("max backtracks reached"
                                  if backtracks >= self.max_backtracks
                                  else "stack empty")
                        print(f"  v={v}: no fix found, {reason}")
                    break

                backtracks += 1
                prev_v, prev_ci, prev_saved = stack.pop()
                self._restore_params(prev_saved)

                restored_gaps = self._compute_gaps_up_to(prev_v)
                for vp in range(1, prev_v + 1):
                    gaps[vp] = restored_gaps[vp]

                v = prev_v
                c_start = prev_ci + 1

                if self.verbose:
                    print(f"  backtrack to v={v} "
                          f"({backtracks}/{self.max_backtracks})")
                    sys.stdout.flush()

        # Final
        gaps = self._compute_all_gaps()
        se = sum(gaps[v] for v in range(1, L + 1))

        elapsed = time.time() - t_start
        if self.verbose:
            print(f"\nOptimization complete: se={se}, "
                  f"{backtracks} backtracks, {elapsed:.1f}s")

        return OptResult(se=se, gaps=gaps, elapsed=elapsed,
                         backtracks=backtracks)

    # ══════════════════════════════════════════════════════════════════════
    # Parameter save/restore
    # ══════════════════════════════════════════════════════════════════════

    def _save_params(self) -> dict:
        saved = {}
        for j, comp in enumerate(self.comb.components):
            for ti, trans in enumerate(comp.trans):
                for param_name, _ in trans.optimizable_params():
                    saved[(j, ti, param_name)] = trans.get_param(param_name)
        return saved

    def _restore_params(self, saved: dict) -> None:
        for (j, ti, param_name), value in saved.items():
            self.comb.components[j].trans[ti].set_param(param_name, value)

    # ══════════════════════════════════════════════════════════════════════
    # Gap computation — dispatches to matricial or lattice method
    # ══════════════════════════════════════════════════════════════════════

    def _compute_all_gaps(self) -> dict[int, int]:
        if self._method == METHOD_LATTICE:
            return self._gaps_lattice()
        return self._gaps_matricial(self.L)

    def _compute_gaps_up_to(self, v_max: int) -> dict[int, int]:
        if self._method == METHOD_LATTICE:
            # Lattice computes all resolutions at once
            all_gaps = self._gaps_lattice()
            return {v: all_gaps.get(v, 0) for v in range(1, v_max + 1)}
        return self._gaps_matricial(v_max)

    def _gaps_matricial(self, v_max: int) -> dict[int, int]:
        comb = self.comb
        gens = [comb[j]._cpp_gen for j in range(comb.J)]
        gen_k = [comb[j].k for j in range(comb.J)]
        trans = [
            [t._cpp_trans for t in comp.trans]
            for comp in comb.components
        ]
        mat = _cpp.prepare_mat(gens, gen_k, trans, comb.k_g,
                               comb.k_g, comb.L)

        kg = comb.k_g
        L = self.L
        gaps = {}
        for v in range(1, v_max + 1):
            t_star = kg // v
            m = mat.copy()
            t_v = m.dimension_equid(kg, v, L)
            gaps[v] = t_star - t_v
        return gaps

    def _gaps_lattice(self) -> dict[int, int]:
        comb = self.comb
        INT_MAX = 2**31 - 1
        gens = [comb[j]._cpp_gen for j in range(comb.J)]
        trans = [
            [t._cpp_trans for t in comp.trans]
            for comp in comb.components
        ]
        delta = [INT_MAX] * (self.L + 1)

        result = _cpp.test_me_lat(
            gens, trans,
            comb.k_g, comb.L, self.L,
            delta, INT_MAX,
        )

        gaps = {}
        ecart = result["ecart"]
        for v in range(1, self.L + 1):
            gaps[v] = ecart[v] if v < len(ecart) else 0
        return gaps


class OptResult:
    """Result of a tempering optimization run."""

    __slots__ = ("se", "gaps", "elapsed", "backtracks")

    def __init__(self, se: int, gaps: dict[int, int],
                 elapsed: float, backtracks: int) -> None:
        self.se = se
        self.gaps = gaps
        self.elapsed = elapsed
        self.backtracks = backtracks
