"""Phase 2.4b-pre: pybind11 binding for Component + Combination.

C++ iterator semantics are exercised in detail by test_combination.cpp.
Here we only check the binding shape and that Python can drive the
iterator end-to-end.
"""

from __future__ import annotations


def _make_tgfsr(a: int):
    from regpoly_cpp._regpoly_cpp import create_generator
    return create_generator(
        "TGFSRGen",
        {"w": 32, "r": 3, "m": 1, "a": a},
        32,
    )


def test_combination_independent_pools_iterates_cartesian() -> None:
    from regpoly_cpp._regpoly_cpp import Combination

    c = Combination(2, 32)
    c.component(0).add_gen(_make_tgfsr(1))
    c.component(0).add_gen(_make_tgfsr(2))
    c.component(1).add_gen(_make_tgfsr(3))
    c.component(1).add_gen(_make_tgfsr(4))

    assert c.reset() is True
    count = 1
    while c.next():
        count += 1
    assert count == 4
    assert c.exhausted() is True


def test_combination_shared_pool_enforces_cnk() -> None:
    from regpoly_cpp._regpoly_cpp import Combination

    c = Combination(2, 32)
    for a in (1, 2, 3, 4):
        c.component(0).add_gen(_make_tgfsr(a))
    c.component(1).share_pool_with(c.component(0))

    assert c.reset() is True
    count = 1
    while c.next():
        count += 1
    # C(4, 2) = 6.
    assert count == 6


def test_combination_kg_l_and_at() -> None:
    from regpoly_cpp._regpoly_cpp import Combination

    c = Combination(1, 32)
    c.component(0).add_gen(_make_tgfsr(0xdeadbeef))
    assert c.reset() is True
    assert c.k_g == c[0].k()
    assert c.L > 0


def test_combination_empty_returns_false() -> None:
    from regpoly_cpp._regpoly_cpp import Combination

    c = Combination(2, 32)
    # No generators added.
    assert c.reset() is False
    assert c.exhausted() is True
