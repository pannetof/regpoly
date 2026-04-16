"""
regpoly — Analysis and search for combined pseudo-random number generators
based on modulo-2 linear recurrences (LFSRs over GF(2)).
"""

__version__ = "2.0.0"

from regpoly.bitvect import BitVect
from regpoly.matrix import BitMatrix
from regpoly.generateur import Generateur
from regpoly.component import Component
from regpoly.combinaison import Combinaison
from regpoly.transformation import Transformation
from regpoly.legacy_reader import LegacyReader
from regpoly.seek import Seek
from regpoly.search_primitive import PrimitiveSearch
from regpoly.tempering_optimizer import TemperingOptimizer
from regpoly.tempering_search import TemperingSearch

__all__ = [
    "BitVect",
    "BitMatrix",
    "Generateur",
    "Component",
    "Combinaison",
    "Transformation",
    "LegacyReader",
    "Seek",
    "PrimitiveSearch",
    "TemperingOptimizer",
    "TemperingSearch",
]
