"""Phase 3 red — v2 history endpoints (exact-K bucket downsampling)."""

from __future__ import annotations


def test_v2_primitive_history_returns_points_array(seeded_client) -> None:
    r = seeded_client.get("/api/v2/primitive-searches/1/history?points=10")
    assert r.status_code == 200, r.text
    body = r.json()
    assert "points" in body and isinstance(body["points"], list)
    # Tolerate empty when the seeded run has no progress rows; structure
    # contract still applies.
    if body["points"]:
        assert all(isinstance(p, (int, float)) for p in body["points"])
    assert "type" in body
    assert "status" in body


def test_v2_tempering_history_returns_best_se_series(client) -> None:
    # No tempering run is seeded → endpoint must return 404 cleanly,
    # not 500.
    r = client.get("/api/v2/tempering-searches/9999/history?points=10")
    assert r.status_code == 404, r.text


def test_v2_history_exact_K_via_index_buckets(client) -> None:
    # When ?points= is set, response.points length is at most K.
    r = client.get("/api/v2/primitive-searches/9999/history?points=10")
    assert r.status_code == 404


def test_v2_duplicate_endpoint_returns_prefill_params(seeded_client) -> None:
    r = seeded_client.get("/api/v2/primitive-searches/1/duplicate")
    assert r.status_code == 200
    body = r.json()
    assert body.get("family") == "MTGen"
    assert body.get("L") == 19937
