"""
component.py — Pool of candidate Generateur instances + Transformation chain.
"""

from __future__ import annotations

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

    def update_trans(self) -> None:
        """
        Synchronise each transformation's width w with the active generator's
        output resolution L, then recompute its internal parameters.
        """
        active_L = self.gens[self.current_gen].L
        for t in self.trans:
            t.update_params(active_L)

    # -- Display ----------------------------------------------------------

    def display(self) -> str:
        """Return a string describing all transformations in this component."""
        return "\n".join(t.display() for t in self.trans)
