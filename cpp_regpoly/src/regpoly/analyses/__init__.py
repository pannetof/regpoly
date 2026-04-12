"""
regpoly.analyses — Statistical test infrastructure (equidistribution, collision-free, tuplets).
"""

from regpoly.analyses.test_base import AbstractTest
from regpoly.analyses.test_results_base import AbstractTestResults
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.analyses.equidistribution_results import EquidistributionResults
from regpoly.analyses.collision_free_test import CollisionFreeTest
from regpoly.analyses.collision_free_results import CollisionFreeResults
from regpoly.analyses.tuplets_test import TupletsTest
from regpoly.analyses.tuplets_results import TupletsResults

__all__ = [
    "AbstractTest",
    "AbstractTestResults",
    "EquidistributionTest",
    "EquidistributionResults",
    "CollisionFreeTest",
    "CollisionFreeResults",
    "TupletsTest",
    "TupletsResults",
]
