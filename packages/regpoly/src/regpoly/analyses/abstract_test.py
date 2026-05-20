# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
test_base.py — Abstract base class for ME/CF tests.

Provides the shared _prepare_mat() static method
used by EquidistributionTest, CollisionFreeTest, and TupletsTest.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import regpoly._regpoly_cpp as _cpp
from regpoly.analyses.abstract_results import AbstractTestResults

if TYPE_CHECKING:
    from regpoly.core.combination import Combination


class AbstractTest(ABC):
    """
    Abstract base for all statistical test classes.

    A test object holds the configuration parameters.  Calling run()
    executes the test against the current combination and returns an
    AbstractTestResults object.
    """

    @abstractmethod
    def run(self, C: "Combination", *args, **kwargs) -> AbstractTestResults:
        """Execute the test on C and return a results object."""

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
    def _prepare_mat(C: "Combination", indice_max: int):
        """
        PrepareMat: build the generator matrix for Gaussian elimination.

        Returns a C++ GaussMatrix. The entire matrix build (including
        transformations) runs in C++.
        """
        gens = []
        gen_k = []
        trans = []
        for j, comp in enumerate(C.components):
            gen = C[j]
            gens.append(gen._cpp_gen)
            gen_k.append(gen.k)
            chain = [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            trans.append(chain)

        return _cpp.prepare_mat(gens, gen_k, trans, C.k_g, indice_max, C.L)
