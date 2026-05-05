"""Phase 5 red — /generators/compare?ids=… side-by-side view."""

from __future__ import annotations


def test_compare_route_with_two_ids(seeded_client) -> None:
    r = seeded_client.get("/generators/compare?ids=1,1")
    assert r.status_code == 200, r.text
    body = r.text
    assert "compare" in body.lower()


def test_v2_compare_endpoint_returns_one_row_per_id(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/compare?ids=1,1")
    assert r.status_code == 200, r.text
    body = r.json()
    assert "rows" in body
    assert len(body["rows"]) == 2


def test_v2_compare_endpoint_skips_missing_ids(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/compare?ids=1,99999")
    assert r.status_code == 200
    body = r.json()
    # Missing ids appear as null entries so the UI can render a hole.
    assert len(body["rows"]) == 2
    assert any(row is None for row in body["rows"])
