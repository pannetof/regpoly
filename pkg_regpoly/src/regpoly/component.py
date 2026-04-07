"""
component.py — Pool of candidate Generateur instances + Transformation chain.
"""

from __future__ import annotations

from regpoly.bitvect import BitVect
from regpoly.generateur import Generateur
from regpoly.transformation import Transformation


class Component:
    """
    One slot of a Combinaison.

    Holds candidate Generateur instances and an associated chain of
    Transformation objects. current_gen selects the active generator.

    len(comp) returns the number of generators in the pool.

    Attributes
    ----------
    gens        : list[Generateur]      — candidate generators
    trans       : list[Transformation]  — tempering chain (applied in order)
    current_gen : int                   — index of the active generator
    """

    def __init__(self) -> None:
        self.gens: list[Generateur] = []
        self.trans: list[Transformation] = []
        self.current_gen: int = 0

    def __len__(self) -> int:
        """Number of generators currently in the pool."""
        return len(self.gens)

    @property
    def nb_gen(self) -> int:
        """Alias for len(self) — used by Combinaison._increment."""
        return len(self.gens)

    # -- Generator pool ---------------------------------------------------

    def add_gen(self, gen: Generateur) -> None:
        """Append a deep copy of gen to the generator pool."""
        self.gens.append(gen.copy())

    # -- Transformation chain ---------------------------------------------

    def add_trans(self, t: Transformation) -> None:
        """Append a deep copy of t to the transformation chain."""
        self.trans.append(t.copy())

    def transform(self, state: BitVect) -> BitVect:
        """Apply all transformations in order to state and return the result."""
        result = state
        for t in self.trans:
            result = t(result)
        return result

    def update_trans(self) -> None:
        """
        Synchronise each transformation's width w with the active generator's
        output resolution L, then recompute its internal parameters.

        If w_original == -1, w inherits L directly; otherwise w = min(w_original, L).
        """
        active_L = self.gens[self.current_gen].L
        for t in self.trans:
            if t.w_original == -1:
                t.w = active_L
            else:
                t.w = min(t.w_original, active_L)
            t.update_params()

    # -- Display ----------------------------------------------------------

    def display(self) -> None:
        """Print all transformations in this component."""
        for t in self.trans:
            t.display()
