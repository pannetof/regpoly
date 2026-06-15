# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
abstract_test.py — Abstract base classes for ME/CF/Tuplets tests
and their result containers.

Provides the shared _prepare_mat() static method
used by EquidistributionTest, CollisionFreeTest, and TupletsTest.
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


class AbstractTest(ABC):
    """
    Abstract base for all statistical test classes.

    A test object holds the configuration parameters.  Calling run()
    executes the test against the supplied generator and returns an
    AbstractTestResults object.
    """

    @abstractmethod
    def run(self, gen, *args, **kwargs) -> AbstractTestResults:
        """Execute the test on ``gen`` and return a results object.

        ``gen`` may be a :class:`~regpoly.core.generator.Generator` wrapper,
        a bare ``_cpp.Recurrence`` / ``_cpp.F2LinearSource``, or a
        ``_cpp.CombinedF2LinearSource`` produced by
        :func:`regpoly.make_combined`. Concrete subclasses normalise it
        via :func:`_to_cpp_gen`.
        """

    @classmethod
    @abstractmethod
    def _from_params(cls, params: dict, Lmax: int) -> "AbstractTest":
        """Construct a test from a parameter dict (YAML)."""

    @classmethod
    def from_yaml(cls, filename: str, Lmax: int) -> list["AbstractTest"]:
        """
        Read tests from a YAML file.

        File format::

            tests:
              - type: equidistribution
                max_gap_sum: 100
                ...
              - type: collision_free
                max_gap_sum: 0
              - type: tuplets
                dimensions: [50, 10, 5]
                ...
              - type: tvalue
                s_max: 8
                max_t_sum: 4
                ...

        Returns a list of AbstractTest instances. New tests are added
        by registering them in
        :mod:`regpoly.analyses.test_registry`; this dispatch site
        consumes the registry rather than carrying its own hard-coded
        table.
        """
        import yaml

        # Late import: the registry module itself imports every test
        # class, so importing it at module top would create a cycle
        # (each test module imports AbstractTest).
        from regpoly.analyses import test_registry

        with open(filename) as f:
            data = yaml.safe_load(f)

        tests = []
        for entry in data["tests"]:
            type_str = entry["type"]
            try:
                spec = test_registry.get(type_str)
            except KeyError as exc:
                raise ValueError(str(exc)) from exc
            tests.append(spec.test_cls._from_params(entry, Lmax))
        return tests

    @staticmethod
    def _to_cpp_gen(gen):
        """Normalize a `run(...)` argument to a `_cpp.CombinedF2LinearSource`.

        Every C++ kernel takes a CombinedF2LinearSource (the lattice
        family requires it for `cs.recurrence_components()`; the matricial
        method accepts it via the `ITestable` base). A bare Generator is
        wrapped in a 1-component CombinedF2LinearSource here.

        Accepts:
        - `regpoly.core.generator.Generator` (Python wrapper) — wrapped
          as a J=1 CombinedF2LinearSource.
        - `_cpp.CombinedF2LinearSource` — returned as-is.
        - `_cpp.F2LinearSource` (any C++ generator instance) — wrapped
          as a J=1 CombinedF2LinearSource.

        Other shapes raise `TypeError`.
        """
        import regpoly._regpoly_cpp as _cpp
        from regpoly.core.generator import Generator
        if isinstance(gen, _cpp.CombinedF2LinearSource):
            return gen
        if isinstance(gen, Generator):
            return _cpp.CombinedF2LinearSource([gen._cpp_gen], gen.L)
        if isinstance(gen, _cpp.F2LinearSource):
            return _cpp.CombinedF2LinearSource([gen], gen.L())
        raise TypeError(
            f"{type(gen).__name__} is neither a regpoly.Generator nor a "
            "_cpp.F2LinearSource handle. Use Generator.create(...) or "
            "make_combined(...) to build a valid input."
        )

    @staticmethod
    def _prepare_mat(gen_or_C, indice_max: int):
        """
        PrepareMat: build the generator matrix for Gaussian elimination.

        Returns a C++ GaussMatrix. Accepts any input shape `_to_cpp_gen`
        understands (Generator, Combination, raw _cpp.Recurrence). The
        unpacking (`components()` / `tempering_chains()`) happens on the
        C++ side via `prepare_mat_from_gen`.
        """
        import regpoly._regpoly_cpp as _cpp
        cpp_gen = AbstractTest._to_cpp_gen(gen_or_C)
        return _cpp.prepare_mat_from_gen(cpp_gen, indice_max)
