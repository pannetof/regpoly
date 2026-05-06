# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Tests for the per-generator default_test_methods map exposed via
GET /api/library/generators/{gen_id}.

The map is computed server-side by instantiating the catalog entry and
asking the underlying C++ Generator for its recommendation per test
type. Expected values reflect the rules in
src/generators/{sfmt,dsfmt,mtgp,rmt64,combined}.cpp and the base
implementation in src/core/generator.cpp:

  - SFMT family (sfmt607)              -> simd_notprimitive
  - Combined generator (taus88, J>1)   -> notprimitive
  - All entries return None for
    collision_free / tuplets in this change.

Exhaustive per-rule coverage lives in the C++ gtest suite
(test_default_test_method.cpp) and the Python binding tests
(test_default_test_method.py). This file checks the endpoint shape
and that the map propagates through the FastAPI route.
"""

from __future__ import annotations

_KNOWN_TEST_TYPES = ("equidistribution", "collision_free", "tuplets")


def _detail(client, gen_id: str) -> dict:
    r = client.get(f"/api/library/generators/{gen_id}")
    assert r.status_code == 200, r.text
    return r.json()


def test_sfmt_default_is_simd_notprimitive(client) -> None:
    body = _detail(client, "sfmt607")
    assert "default_test_methods" in body
    methods = body["default_test_methods"]
    assert set(methods.keys()) == set(_KNOWN_TEST_TYPES)
    assert methods["equidistribution"] == "simd_notprimitive"
    assert methods["collision_free"] is None
    assert methods["tuplets"] is None


def test_combined_tausworthe_default_is_notprimitive(client) -> None:
    body = _detail(client, "taus88")
    methods = body["default_test_methods"]
    assert methods["equidistribution"] == "notprimitive"
    assert methods["collision_free"] is None
    assert methods["tuplets"] is None


def test_full_period_mt_default_is_harase(client) -> None:
    """mt19937 is a single MT with k=19937 (>> 100), full-period.
    Base rule: full-period & k > 100 -> harase."""
    body = _detail(client, "mt19937")
    methods = body["default_test_methods"]
    assert methods["equidistribution"] == "harase"


def test_unknown_generator_404(client) -> None:
    r = client.get("/api/library/generators/this-does-not-exist")
    assert r.status_code == 404


def test_default_methods_cache_is_warm_on_repeat(client, monkeypatch) -> None:
    """Two GETs on the same generator must call _default_test_methods_for
    at most ONCE — the (gen_id, source_mtime)-keyed cache short-circuits
    the second request.
    """
    from regpoly_web.routes import library as lib_routes

    real_compute = lib_routes._default_test_methods_for
    call_count = {"n": 0}

    def counting(g):
        call_count["n"] += 1
        return real_compute(g)

    monkeypatch.setattr(lib_routes, "_default_test_methods_for", counting)

    # Cold cache: first request triggers compute.
    r1 = client.get("/api/library/generators/sfmt607")
    assert r1.status_code == 200
    assert call_count["n"] == 1

    # Warm cache: subsequent requests must NOT recompute.
    for _ in range(5):
        r = client.get("/api/library/generators/sfmt607")
        assert r.status_code == 200
        # All return the same cached payload.
        assert r.json()["default_test_methods"] == r1.json()["default_test_methods"]
    assert call_count["n"] == 1, (
        f"cache miss: compute called {call_count['n']} times across 6 requests"
    )
