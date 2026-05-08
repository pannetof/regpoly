# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""WELL-RNG helpers (paper Table I cost model, structured `matrices` map).

A thin Python facade over the C++ free functions
`well_random_matrices` and `well_total_cost`. The C++ side owns the
sampling algorithm (rejection + greedy-budgeted fallback) and the
M-class cost table. Both the CLI search driver and the web worker
consume this module so the two paths stay in lockstep.
"""

from __future__ import annotations

from typing import Mapping

from regpoly_cpp._regpoly_cpp import (
    well_random_matrices as _random_matrices,
    well_total_cost as _total_cost,
)


def random_matrices(w: int, max_cost: int, seed: int = 0) -> dict:
    """Sample a WELL `matrices` map whose total cost is ≤ ``max_cost``.

    Slots ``T0..T7`` each carry an ``M``-class (paper Table I) plus the
    args required for that class. The sum of per-Mi costs (M0=0, M1=1,
    M2=2, M3=3, M4=5, M5=4, M6=8) is bounded by ``max_cost``.

    Algorithm: 64-attempt rejection sampling first; on failure, falls
    back to a greedy-budgeted draw that always succeeds. Returns the
    nested dict shape consumed by ``Generator.create("WELLRNG", ...)``.

    ``seed`` is consumed by a per-call RNG so the function is
    reproducible without touching global state.

    Raises ``ValueError`` (mapped from C++ ``std::invalid_argument``)
    when ``max_cost <= 0`` or ``w != 32``.
    """
    return _random_matrices(w, max_cost, seed)


def total_cost(matrices: Mapping[str, Mapping[str, object]]) -> int:
    """Sum of per-Mi costs across the slots in a `matrices` map."""
    return _total_cost(dict(matrices))


__all__ = ["random_matrices", "total_cost"]
