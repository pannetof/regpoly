"""
regpoly.tests — Statistical test infrastructure (equidistribution, collision-free, tuplets).
"""

from regpoly.tests.test_base import AbstractTest
from regpoly.tests.test_results_base import AbstractTestResults
from regpoly.tests.equidistribution_test import EquidistributionTest
from regpoly.tests.equidistribution_results import EquidistributionResults
from regpoly.tests.collision_free_test import CollisionFreeTest
from regpoly.tests.collision_free_results import CollisionFreeResults
from regpoly.tests.tuplets_test import TupletsTest
from regpoly.tests.tuplets_results import TupletsResults

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
