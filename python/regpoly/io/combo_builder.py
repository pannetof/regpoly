# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
combo_builder.py — direct YAML-to-C++ enumerator builder.

Replaces the previous `seek.py:_build_cpp_comb_from_python` middleman.
`Seek.from_yaml` and `TemperedSearch.from_yaml` now go straight from
per-slot generator/tempering pools to a `_cpp.ComboEnumerator` (the
C++ search-loop iterator that Phase 6 renames to `ComboEnumerator`).

No Python `ComboEnumerator` is constructed along the way.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import regpoly._regpoly_cpp as _cpp

if TYPE_CHECKING:
    from regpoly.core.generator import Generator
    from regpoly.core.transformation import Transformation


def build_cpp_enumerator(
    gen_pools: list[list["Generator"]],
    tempering_pools: list[list["Transformation"]],
    Lmax: int,
) -> "_cpp.ComboEnumerator":
    """Build a `_cpp.ComboEnumerator` search iterator from per-slot pools.

    Parameters
    ----------
    gen_pools
        Per-slot list of candidate :class:`Generator` instances.
        Sharing the same Python list object across two slots replicates
        the legacy `copy_pool_from` semantics (C(n,k) selection).
    tempering_pools
        Per-slot ordered list of :class:`Transformation` instances
        (the tempering chain for that slot). Each transformation must
        carry a `_cpp_trans` attribute; missing it is a hard failure
        (closes the silent-drop footgun the previous list-comprehension
        filter hid).
    Lmax
        Maximum output bit-width for the combined point set.

    Returns
    -------
    `_cpp.ComboEnumerator`
        A fresh search iterator, with `reset()` already called. Pass
        directly to `_cpp.run_seek_search` / `_cpp.run_tempering_search`.

    Raises
    ------
    RuntimeError
        If `reset()` returns False (typically because one pool is empty).
    TypeError
        If a tempering item lacks `_cpp_trans`.
    """
    if len(gen_pools) != len(tempering_pools):
        raise ValueError(
            f"gen_pools has {len(gen_pools)} slots but tempering_pools has "
            f"{len(tempering_pools)}; they must match."
        )

    nb_comp = len(gen_pools)
    cpp_comb = _cpp.ComboEnumerator(nb_comp, Lmax)

    # Identity-based shared-pool detection: two slots that point at the
    # same Python list object share their pool on the C++ side.
    pool_owner: dict[int, int] = {}
    for j, gen_list in enumerate(gen_pools):
        cpp_comp = cpp_comb.pool(j)
        pool_key = id(gen_list)
        if pool_key in pool_owner:
            cpp_comp.copy_pool_from(cpp_comb.pool(pool_owner[pool_key]))
        else:
            pool_owner[pool_key] = j
            for gen in gen_list:
                cpp_comp.add_source(gen._cpp_gen)
        for trans in tempering_pools[j]:
            cpp_t = getattr(trans, "_cpp_trans", None)
            if cpp_t is None:
                raise TypeError(
                    f"{type(trans).__name__} has no `_cpp_trans`; "
                    "transformations passed to build_cpp_enumerator must "
                    "have a registered C++ counterpart."
                )
            cpp_comp.add_trans(cpp_t)

    if not cpp_comb.reset():
        raise RuntimeError(
            "build_cpp_enumerator: search space is empty (one of the "
            "generator pools has no candidates)."
        )
    return cpp_comb
