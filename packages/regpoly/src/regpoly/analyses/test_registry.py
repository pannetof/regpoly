# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
test_registry.py — single source of truth for "what tests exist".

Replaces the hard-coded `_dispatch` dict that previously lived in
`abstract_test.py`. Adding a new statistical test requires exactly one
new entry here:

    >>> from regpoly.analyses.my_new_test import MyNewTest, MyNewResults
    >>> register(TestSpec(
    ...     name="my_new",
    ...     test_cls=MyNewTest,
    ...     results_cls=MyNewResults,
    ...     cpp_seek_kind="MyNew",     # if a C++ search-runner dispatch
    ...                                # case exists; else None
    ... ))

`AbstractTest.from_yaml` looks tests up by `name`. The optional
`cpp_seek_kind` field is what `seek.py` should use to dispatch to the
C++ `SeekTestKind` enum (today the enum is still hard-coded; refactor
follow-up tracked alongside this file).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional, Type

if TYPE_CHECKING:
    from regpoly.analyses.abstract_results import AbstractTestResults
    from regpoly.analyses.abstract_test import AbstractTest


@dataclass(frozen=True)
class TestSpec:
    """One row of the registry.

    Attributes
    ----------
    name
        YAML ``type:`` value used in ``tests:`` blocks.
    test_cls
        The :class:`AbstractTest` subclass to instantiate.
    results_cls
        The :class:`AbstractTestResults` subclass produced by
        ``test_cls.run``.
    cpp_seek_kind
        Name of the corresponding ``_cpp.SeekTestKind`` enum value, or
        ``None`` if the test does not (yet) flow through the C++ search
        runner. Used by ``seek.py`` when it wires Python tests into the
        C++ search loop; ``None`` means "Python-only orchestration".
    """

    name: str
    test_cls: Type["AbstractTest"]
    results_cls: Type["AbstractTestResults"]
    cpp_seek_kind: Optional[str] = None


_REGISTRY: dict[str, TestSpec] = {}


def register(spec: TestSpec) -> None:
    """Add a test to the registry.

    Re-registration with the same ``name`` is allowed (used by tests
    that wipe and re-populate the registry); the new spec replaces the
    old.
    """
    _REGISTRY[spec.name] = spec


def get(name: str) -> TestSpec:
    """Look up a test by YAML name. Raises ``KeyError`` if absent."""
    try:
        return _REGISTRY[name]
    except KeyError:
        raise KeyError(
            f"Unknown test type '{name}'. Registered: {sorted(_REGISTRY)}."
        ) from None


def all_specs() -> list[TestSpec]:
    """Return every registered :class:`TestSpec`, sorted by name."""
    return [_REGISTRY[name] for name in sorted(_REGISTRY)]


def names() -> list[str]:
    """Return the registered YAML test-type names, sorted."""
    return sorted(_REGISTRY)


def _populate_default_registry() -> None:
    """Register the canonical tests shipped with regpoly.

    Called on first import of this module. Idempotent — repeat calls
    are harmless because registration is a dict write.
    """
    # Late imports: each test module pulls in heavy machinery; we avoid
    # paying the cost until something actually queries the registry.
    from regpoly.analyses.collision_free_results import CollisionFreeResults
    from regpoly.analyses.collision_free_test import CollisionFreeTest
    from regpoly.analyses.equidistribution_results import EquidistributionResults
    from regpoly.analyses.equidistribution_test import EquidistributionTest
    from regpoly.analyses.tuplets_results import TupletsResults
    from regpoly.analyses.tuplets_test import TupletsTest
    from regpoly.analyses.tvalue_results import TValueResults
    from regpoly.analyses.tvalue_test import TValueTest

    register(TestSpec(
        name="equidistribution",
        test_cls=EquidistributionTest,
        results_cls=EquidistributionResults,
        cpp_seek_kind="Equidistribution",
    ))
    register(TestSpec(
        name="collision_free",
        test_cls=CollisionFreeTest,
        results_cls=CollisionFreeResults,
        cpp_seek_kind="CollisionFree",
    ))
    register(TestSpec(
        name="tuplets",
        test_cls=TupletsTest,
        results_cls=TupletsResults,
        cpp_seek_kind="Tuplets",
    ))
    register(TestSpec(
        name="tvalue",
        test_cls=TValueTest,
        results_cls=TValueResults,
        cpp_seek_kind=None,   # search-runner integration is a follow-up
    ))


_populate_default_registry()
