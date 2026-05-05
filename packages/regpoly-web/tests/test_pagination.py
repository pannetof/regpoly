"""Phase 2 red — server-side pagination on /api/v2/ list endpoints.

Every v2 list endpoint follows the same contract:
  GET /api/v2/{entity}?...&limit&offset
  → { rows: [...], total: int }

Unified across generators, tested-generators, library, library/papers
(per the developer-pillar S10.4 mandate).
"""

from __future__ import annotations


def test_v2_generators_list_limit_offset(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators?limit=10&offset=0")
    assert r.status_code == 200
    body = r.json()
    assert "rows" in body, "v2 list shape uses `rows` (not `items`)"
    assert "total" in body
    assert isinstance(body["rows"], list)
    assert isinstance(body["total"], int)
    assert "items" not in body, "v2 must not carry the v1 `items` key"


def test_v2_pagination_total_count_in_response(seeded_client) -> None:
    """`total` reflects the unfiltered match count, not the page size."""
    body = seeded_client.get("/api/v2/generators?limit=1&offset=0").json()
    # seeded DB has at least one primitive_generator
    assert body["total"] >= 1
    assert len(body["rows"]) <= 1


def test_v2_tested_generators_list_limit_offset(seeded_client) -> None:
    r = seeded_client.get("/api/v2/tested-generators?limit=10&offset=0")
    assert r.status_code == 200
    body = r.json()
    assert "rows" in body and "total" in body
    assert isinstance(body["rows"], list)


def test_v2_library_list_limit_offset(seeded_client) -> None:
    """Per the Developer-pillar D-2 demand: unified pagination across
    every v2 list endpoint, including the library."""
    r = seeded_client.get("/api/v2/library?limit=10&offset=0")
    assert r.status_code == 200
    body = r.json()
    assert "rows" in body and "total" in body


def test_v2_library_papers_list_limit_offset(seeded_client) -> None:
    r = seeded_client.get("/api/v2/library/papers?limit=10&offset=0")
    assert r.status_code == 200
    body = r.json()
    assert "rows" in body and "total" in body


def test_v2_offset_actually_skips_rows(seeded_client) -> None:
    """offset=N skips the first N matching rows."""
    full = seeded_client.get("/api/v2/generators?limit=100&offset=0").json()
    if full["total"] < 2:
        # Only one seeded row — verify offset=1 returns empty rows.
        r = seeded_client.get("/api/v2/generators?limit=100&offset=1").json()
        assert r["rows"] == []
        return
    skipped = seeded_client.get("/api/v2/generators?limit=100&offset=1").json()
    assert len(skipped["rows"]) == len(full["rows"]) - 1
    assert skipped["rows"][0]["id"] != full["rows"][0]["id"]
