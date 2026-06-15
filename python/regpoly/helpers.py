# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
helpers.py — user-facing convenience wrappers.

Replaces the four-line `Combination(J=…, Lmax=…)` + `add_gen` + `reset`
ceremony with a single function call when the goal is to build a
single runtime generator (one-component or J-component XOR-combined)
to feed to an `AbstractTest.run(...)`.

The implementation routes directly to the C++ `CombinedF2LinearSource`
constructor — no Python `Combination` is involved.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import regpoly._regpoly_cpp as _cpp

if TYPE_CHECKING:
    from regpoly.core.generator import Generator
    from regpoly.core.transformation import Transformation


def make_combined(
    *gens: "Generator",
    trans: list["Transformation"] | list[list["Transformation"]] | None = None,
    Lmax: int | None = None,
) -> "_cpp.CombinedF2LinearSource":
    """Build a `_cpp.CombinedF2LinearSource` ready to pass to `*.run(...)`.

    The single-generator (`make_combined(gen)`) shape replaces
    `Combination.single(gen)`. The J-component shape
    (`make_combined(g1, g2, trans=[chain1, chain2])`) replaces
    `Combination.CreateFromFiles([[g1], [g2]], Lmax, [chain1, chain2])`
    when the pools-per-slot are size 1.

    Parameters
    ----------
    *gens
        One or more :class:`regpoly.core.generator.Generator` instances.
        At least one is required.
    trans
        Tempering chains. Two accepted shapes:

        - ``None`` (default): no tempering on any component.
        - A flat list of transformations: applied to the *first*
          generator only. Useful when ``len(gens) == 1``. For multi-
          component cases this raises ``ValueError``.
        - A list-of-lists of transformations: ``trans[j]`` is the chain
          for ``gens[j]``. Length must match ``len(gens)``.

        Every transformation must have a ``_cpp_trans`` attribute; the
        check is *strict* (a transformation lacking ``_cpp_trans``
        raises ``TypeError`` — silent drops are not tolerated, mirroring
        the previous list-comprehension's hidden footgun).
    Lmax
        Output bit-width. Defaults to ``min(g.L for g in gens)`` —
        the standard "combined L is the minimum across active gens"
        rule that the search loop has always used.

    Returns
    -------
    `_cpp.CombinedF2LinearSource`
        A runtime generator handle accepted by every `AbstractTest.run`
        and by every C++ kernel taking a `Generator&`.
    """
    if not gens:
        raise ValueError("make_combined: at least one Generator is required")

    cpp_gens = [g._cpp_gen for g in gens]

    if trans is None:
        cpp_trans: list[list] = [[] for _ in gens]
    elif _is_flat_trans_list(trans):
        if len(gens) != 1:
            raise ValueError(
                "make_combined: a flat `trans` list is only valid when "
                "len(gens) == 1. For J > 1 pass a list-of-lists."
            )
        cpp_trans = [[_unwrap_trans(t) for t in trans]]
    else:
        if len(trans) != len(gens):
            raise ValueError(
                f"make_combined: trans has {len(trans)} chains but "
                f"{len(gens)} generators."
            )
        cpp_trans = [[_unwrap_trans(t) for t in chain] for chain in trans]

    L = Lmax if Lmax is not None else min(g.L for g in gens)

    # `_cpp.CombinedF2LinearSource` accepts any F2LinearSource components
    # (Recurrence PRNGs and DigitalNet sources alike). Kernels that need
    # Recurrence specifically (matricial χ-recovery via `recover_char_poly`,
    # SIMD path) check at call time and throw `std::invalid_argument` if
    # a DigitalNet is present.
    return _cpp.CombinedF2LinearSource(cpp_gens, cpp_trans, L)


def _is_flat_trans_list(trans) -> bool:
    """True iff `trans` is a flat list of Transformation-like objects."""
    if not trans:
        return False
    first = trans[0]
    return hasattr(first, "_cpp_trans") or not hasattr(first, "__iter__")


def _unwrap_trans(t):
    """Return `t._cpp_trans`, raising `TypeError` with context if absent."""
    cpp_t = getattr(t, "_cpp_trans", None)
    if cpp_t is None:
        raise TypeError(
            f"{type(t).__name__} has no `_cpp_trans`; transformation "
            "wrappers used in `make_combined(...)` must have a registered "
            "C++ counterpart."
        )
    return cpp_t
