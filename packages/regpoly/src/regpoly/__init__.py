# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
regpoly — Analysis and search for combined pseudo-random number generators
based on modulo-2 linear recurrences (LFSRs over GF(2)).
"""

__version__ = "2.0.0"

# The C++ extension lives in the regpoly-cpp package as
# `regpoly_cpp._regpoly_cpp`. Historical Python code imports it as
# `regpoly._regpoly_cpp` (this aliasing predates the monorepo split).
# Until the v2.0 redesign Phase 2 reroutes every call site through a
# proper `regpoly.introspection` / wrapper layer, expose the module
# under both names via sys.modules so existing imports keep working.
import sys as _sys

import regpoly_cpp._regpoly_cpp as _regpoly_cpp_module  # noqa: E402

_sys.modules.setdefault("regpoly._regpoly_cpp", _regpoly_cpp_module)
del _sys, _regpoly_cpp_module

from regpoly.core.bitvect import BitVect
from regpoly.core.combination import Combination
from regpoly.core.component import Component
from regpoly.core.generator import Generator
from regpoly.core.matrix import BitMatrix
from regpoly.core.transformation import Transformation
from regpoly.search.search_primitive import PrimitiveSearch
from regpoly.search.seek import Seek
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
    "Seek",
    "PrimitiveSearch",
    "TemperingOptimizer",
    "TemperingSearch",
    # legacy aliases
    "Generateur",
    "Combinaison",
]
