"""
test_results_base.py — Abstract base class for ME/CF test results.
"""

from __future__ import annotations

from abc import ABC, abstractmethod


class AbstractTestResults(ABC):
    """
    Abstract base for all test result containers.

    A results object is produced by a test's run() method.  It is
    immutable after construction: predicates read the stored data,
    display() prints a summary to stdout.

    Attributes
    ----------
    verified : bool — True when the test has been successfully executed
                      and the stored values are meaningful.
    """

    @property
    @abstractmethod
    def verified(self) -> bool:
        """True iff the test has been run and results are available."""

    @abstractmethod
    def display(self) -> str:
        """Return a human-readable summary of the results."""
