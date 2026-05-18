# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Primitive-search route coverage for ``XoroshiroGen`` and
``XoshiroGen`` (Blackman & Vigna 2022).

Both families take structural params ``w`` and ``r``; the rest of the
parameter space is one or two non-structural shifts (A, B, [C]) in
[1, w-1] that can be fixed, randomised, or exhaustively enumerated.
"""

from __future__ import annotations


# ── /api/families — both families are discoverable ───────────────────────


def test_families_endpoint_lists_xoroshiro_and_xoshiro(client) -> None:
    r = client.get("/api/families")
    assert r.status_code == 200
    names = {f["name"] for f in r.json()}
    assert "XoroshiroGen" in names
    assert "XoshiroGen" in names


def test_family_detail_marks_xoroshiro_enumerable(client) -> None:
    r = client.get("/api/families/XoroshiroGen")
    assert r.status_code == 200
    body = r.json()
    assert body["enumerable"] is True
    names = [p["name"] for p in body["params"]]
    # w and r are structural; A B C are non-structural with rand_type=range
    assert names == ["w", "r", "A", "B", "C"]


def test_family_detail_marks_xoshiro_enumerable(client) -> None:
    r = client.get("/api/families/XoshiroGen")
    assert r.status_code == 200
    body = r.json()
    assert body["enumerable"] is True
    names = [p["name"] for p in body["params"]]
    assert names == ["w", "r", "A", "B"]


# ── Exhaustive-mode estimate: enumerator size + axes ─────────────────────


def _xoroshiro_body(*, w: int, r: int, mode: str, fixed=None) -> dict:
    return {
        "family": "XoroshiroGen",
        "structural_params": {"w": w, "r": r},
        "fixed_params": fixed or {},
        "max_tries": None,
        "max_seconds": None,
        "max_cost": None,
        "max_se": None,
        "search_mode": mode,
        "confirm_huge": False,
    }


def _xoshiro_body(*, w: int, r: int, mode: str, fixed=None) -> dict:
    return {
        "family": "XoshiroGen",
        "structural_params": {"w": w, "r": r},
        "fixed_params": fixed or {},
        "max_tries": None,
        "max_seconds": None,
        "max_cost": None,
        "max_se": None,
        "search_mode": mode,
        "confirm_huge": False,
    }


def test_xoroshiro_exhaustive_estimate(client) -> None:
    """Total = (w-1)^3 — axes A, B, C each over [1, w-1]."""
    r = client.post(
        "/api/primitive-searches/estimate",
        json=_xoroshiro_body(w=16, r=4, mode="exhaustive"),
    )
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["search_mode"] == "exhaustive"
    assert int(data["total"]) == 15 ** 3
    assert [a["name"] for a in data["axes"]] == ["A", "B", "C"]


def test_xoshiro_exhaustive_estimate(client) -> None:
    """Total = (w-1)^2 — axes A, B."""
    r = client.post(
        "/api/primitive-searches/estimate",
        json=_xoshiro_body(w=32, r=4, mode="exhaustive"),
    )
    assert r.status_code == 200, r.text
    data = r.json()
    assert int(data["total"]) == 31 ** 2
    assert [a["name"] for a in data["axes"]] == ["A", "B"]


# ── Create (persist): exhaustive ─────────────────────────────────────────


def test_xoroshiro_exhaustive_create_persists_with_correct_k_and_L(
    client,
) -> None:
    """End-to-end: exhaustive create must build the enumerator, probe
    k from at(0), and store k = w*r with L = w."""
    r = client.post(
        "/api/primitive-searches",
        json=_xoroshiro_body(w=16, r=8, mode="exhaustive"),
    )
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 16 * 8
    assert payload["search_mode"] == "exhaustive"
    assert payload["enum_total"] == str(15 ** 3)
    assert payload["L"] == 16  # one w-bit word per output


def test_xoshiro_exhaustive_create_persists(client) -> None:
    r = client.post(
        "/api/primitive-searches",
        json=_xoshiro_body(w=32, r=4, mode="exhaustive"),
    )
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 32 * 4
    assert payload["enum_total"] == str(31 ** 2)
    assert payload["L"] == 32


# ── Create: random mode (structural-only payload) ────────────────────────


def test_xoroshiro_random_create_with_no_fixed_params(client) -> None:
    """Random mode with only w and r pinned: A, B, C are auto-sampled
    via the range sampler each try."""
    r = client.post(
        "/api/primitive-searches",
        json=_xoroshiro_body(w=32, r=2, mode="random") | {
            "max_tries": 1, "max_seconds": 5,
        },
    )
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 32 * 2
    assert payload["search_mode"] == "random"


def test_xoshiro_random_create_partial_fix(client) -> None:
    """Mixed fix/random: A pinned, B sampled."""
    body = _xoshiro_body(
        w=64, r=4, mode="random", fixed={"A": 17},
    ) | {"max_tries": 1, "max_seconds": 5}
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 64 * 4


def test_form_renders_with_xoroshiro_bv_hiding_helpers(client) -> None:
    """The form template carries the JS helpers needed to hide
    non-structural inputs in exhaustive mode for xoroshiro/xoshiro."""
    r = client.get("/primitive-search?family=XoroshiroGen")
    assert r.status_code == 200
    body = r.text
    for symbol in (
        "isBlackmanVignaFamily",
        "bvParamEnumerated_",
        # the non-structural card is gated on visibleNonStructural().length
        "visibleNonStructural().length",
    ):
        assert symbol in body, f"missing JS symbol: {symbol}"


def test_xoroshiro_random_create_all_fixed_matches_table5(client) -> None:
    """All three shifts pinned to the published Table 5 row: should
    succeed and produce a row with k=64.  The runtime then validates
    that this particular triple is primitive (it is — weight 31)."""
    body = _xoroshiro_body(
        w=32, r=2, mode="random",
        fixed={"A": 26, "B": 9, "C": 13},
    ) | {"max_tries": 1, "max_seconds": 5}
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 64
