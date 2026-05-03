"""
component.py — Pool of candidate Generator instances + Transformation chain.
"""

from __future__ import annotations

from regpoly.core.generator import Generator
from regpoly.core.transformation import Transformation


class Component:
    """
    One slot of a Combination.

    Holds candidate Generator instances and an associated chain of
    Transformation objects. current_gen selects the active generator.

    len(comp) returns the number of generators in the pool.

    Attributes
    ----------
    gens        : list[Generator]      — candidate generators
    trans       : list[Transformation]  — tempering chain (applied in order)
    current_gen : int                   — index of the active generator
    """

    def __init__(self) -> None:
        self.gens: list[Generator] = []
        self.trans: list[Transformation] = []
        self.current_gen: int = 0

    def __len__(self) -> int:
        """Number of generators currently in the pool."""
        return len(self.gens)

    @property
    def nb_gen(self) -> int:
        """Alias for len(self) — used by Combination._increment."""
        return len(self.gens)

    # -- Generator pool ---------------------------------------------------

    def add_gen(self, gen: Generator) -> None:
        """Append a deep copy of gen to the generator pool."""
        self.gens.append(gen.copy())

    # -- Transformation chain ---------------------------------------------

    def add_trans(self, t: Transformation) -> None:
        """Append a deep copy of t to the transformation chain."""
        self.trans.append(t.copy())

    # -- Display ----------------------------------------------------------

    def display(self) -> str:
        """Return a string describing all transformations in this component."""
        return "\n".join(t.display() for t in self.trans)
