"""
regpoly — Analysis and search for combined pseudo-random number generators
based on modulo-2 linear recurrences (LFSRs over GF(2)).
"""

__version__ = "2.0.0"

from regpoly.core.bitvect import BitVect
from regpoly.core.matrix import BitMatrix
from regpoly.core.generator import Generator
from regpoly.core.component import Component
from regpoly.core.combination import Combination
from regpoly.core.transformation import Transformation
from regpoly.io.legacy_reader import LegacyReader
from regpoly.search.seek import Seek
from regpoly.search.search_primitive import PrimitiveSearch
from regpoly.search.tempering_optimizer import TemperingOptimizer
from regpoly.search.tempering_search import TemperingSearch

# Backwards-compat aliases for the renamed French nouns. Same class
# objects, so isinstance/pickle keep working.
Generateur = Generator
Combinaison = Combination

__all__ = [
    "BitVect",
    "BitMatrix",
    "Generator",
    "Component",
    "Combination",
    "Transformation",
    "LegacyReader",
    "Seek",
    "PrimitiveSearch",
    "TemperingOptimizer",
    "TemperingSearch",
    # legacy aliases
    "Generateur",
    "Combinaison",
]
