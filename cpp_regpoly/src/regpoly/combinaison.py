"""
combinaison.py — XOR combination of J Component instances.

Combinaison is a Python iterator over all valid combinations:

    for _ in comb:          # __iter__ calls reset(), __next__ advances
        process(comb)

Or manually:

    if comb.reset():
        process(comb)
        while True:
            try:
                next(comb)
                process(comb)
            except StopIteration:
                break

A "valid" combination is one where no two active generators share the same
degree k.  comb[j] returns the active Generateur of component j.
"""

from __future__ import annotations

import sys

from regpoly.generateur import Generateur
from regpoly.component import Component


def product_no_repeat_ordered(*lists):
    """
    Cartesian product across lists with two constraints:
    1. No object appears twice in a single tuple (identity-based).
    2. When two lists are identical (same object), indices must be
       strictly increasing — enforces C(n,k) selection over permutations.
    """
    def recurse(depth, current_ids, current_indices, current):
        if depth == len(lists):
            yield tuple(current)
            return
        current_list = lists[depth]

        min_index = 0
        for prev_depth in range(depth):
            if lists[prev_depth] is current_list:
                min_index = max(min_index, current_indices[prev_depth] + 1)

        for i in range(min_index, len(current_list)):
            item = current_list[i]
            if id(item) not in current_ids:
                yield from recurse(
                    depth + 1,
                    current_ids | {id(item)},
                    current_indices + [i],
                    current + [item],
                )

    yield from recurse(0, set(), [], [])


class Combinaison:
    """
    An XOR combination of J Component instances.

    The current combination is identified by the tuple
    (components[j].current_gen for j in range(J)).

    Attributes
    ----------
    components : list[Component] — the J component slots
    J          : int  — number of components
    k_g        : int  — sum of active generator degrees
    L          : int  — minimum L across active generators, capped at Lmax
    Lmax       : int  — maximum allowed L
    """

    @classmethod
    def CreateFromFiles(
        cls,
        gen_lists: list,
        Lmax: int,
        temperings: list,
    ) -> "Combinaison":
        """
        Build a fully-populated Combinaison from pre-built data.

        Parameters
        ----------
        gen_lists  : list[list[Generateur]] — pre-built generator pool per component;
                     when two entries are the same object, that component pair shares
                     a pool and C(n,k) selection is enforced automatically.
                     nb_comp is derived from len(gen_lists).
        Lmax       : maximum bit-width
        temperings : list[list[Transformation]] — one list per component
                     (empty list if no tempering)

        Returns
        -------
        comb : Combinaison — fully initialised and populated
        """
        nb_comp = len(gen_lists)
        comb = cls(J=nb_comp, Lmax=Lmax)

        for j, trans_list in enumerate(temperings):
            for t in trans_list:
                comb.components[j].add_trans(t)

        for j, gen_list in enumerate(gen_lists):
            if j > 0 and gen_list is gen_lists[j - 1]:
                comb.components[j].gens = comb.components[j - 1].gens
            else:
                for gen in gen_list:
                    comb.components[j].add_gen(gen)
        return comb

    def __init__(self, J: int = 0, Lmax: int = 0) -> None:
        self.J: int = J
        self.Lmax: int = Lmax
        self.k_g: int = 0
        self.L: int = 0
        self.components: list[Component] = [Component() for _ in range(J)]
        self._combo_iter = iter([])

    # -- Indexing ---------------------------------------------------------

    def __getitem__(self, j: int) -> Generateur:
        """Return the active generator of component j  (comb[j])."""
        comp = self.components[j]
        return comp.gens[comp.current_gen]

    # -- Iterator protocol ------------------------------------------------

    def _make_combo_iter(self):
        return (
            combo
            for combo in product_no_repeat_ordered(*[comp.gens for comp in self.components])
            if len({gen.k for gen in combo}) == len(combo)
        )

    def reset(self) -> bool:
        """
        Reset to the first valid combination.

        Initialises the internal iterator and advances to the first combo.
        Returns True on success, False if no valid combination exists.
        Use next(comb) to advance to subsequent combinations.
        """
        if not any(len(comp) > 0 for comp in self.components):
            return False
        self._combo_iter = self._make_combo_iter()
        try:
            self._apply_combo(next(self._combo_iter))
            return True
        except StopIteration:
            return False

    def __iter__(self) -> "Combinaison":
        """Set up a fresh iterator; the first next() call yields combo 1."""
        self._combo_iter = self._make_combo_iter()
        return self

    def __next__(self) -> "Combinaison":
        self._apply_combo(next(self._combo_iter))   # StopIteration propagates
        return self

    # -- Internal helpers -------------------------------------------------

    def _apply_combo(self, combo: tuple) -> None:
        """Set current_gen indices from a combo tuple and recompute stats."""
        for j, gen in enumerate(combo):
            self.components[j].current_gen = self.components[j].gens.index(gen)
        self._update_stats()

    def _update_stats(self) -> None:
        """Recompute k_g and L from the current active generators."""
        self.k_g = 0
        min_L = sys.maxsize
        for comp in self.components:
            active = comp.gens[comp.current_gen]
            self.k_g += active.k
            if active.L < min_L:
                min_L = active.L
        self.L = self.Lmax if min_L > self.Lmax else min_L
