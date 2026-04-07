"""
regpoly — Analysis and search for combined pseudo-random number generators
based on modulo-2 linear recurrences (LFSRs over GF(2)).
"""

__version__ = "2.0.0"

from regpoly.bitvect import BitVect
from regpoly.matrix import Matrix
from regpoly.generateur import Generateur
from regpoly.component import Component
from regpoly.combinaison import Combinaison
from regpoly.transformation import Transformation
from regpoly.seek import Seek

__all__ = [
    "BitVect",
    "Matrix",
    "Generateur",
    "Component",
    "Combinaison",
    "Transformation",
    "Seek",
]
