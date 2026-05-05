"""Phase 2 regression guard — the existing /api/ namespace contract is
preserved unchanged across the v2 redesign. v2 (renames, wrapper shape
changes) lives under /api/v2/. External CLI / SDK consumers continue
to read {items, total, page, per_page} from /api/{generators,tested-generators}
without any rewrite.

If any of these tests fails, the v2 redesign has leaked a breaking
change into the v1 contract — revert before continuing.
"""

from __future__ import annotations


def test_api_generators_returns_items_total_page_per_page(seeded_client) -> None:
    r = seeded_client.get("/api/generators")
    assert r.status_code == 200
    body = r.json()
    assert "items" in body, "v1 list shape uses `items` (not `rows`)"
    assert "total" in body
    assert "page" in body
    assert "per_page" in body
    # Reject silent wrapper renames.
    assert "rows" not in body, "v1 must NOT carry the v2 `rows` key"
    assert "limit" not in body, "v1 must NOT carry the v2 `limit` key"
    assert "offset" not in body, "v1 must NOT carry the v2 `offset` key"


def test_api_tested_generators_returns_items_total_page_per_page(
    seeded_client,
) -> None:
    r = seeded_client.get("/api/tested-generators")
    assert r.status_code == 200
    body = r.json()
    assert "items" in body
    assert "total" in body
    assert "page" in body
    assert "per_page" in body
    assert "rows" not in body
    assert "limit" not in body


def test_api_generators_pagination_uses_page_per_page(seeded_client) -> None:
    """v1 paginates via ?page=1&per_page=50 (1-indexed). The v2
    endpoint takes ?limit&offset (0-indexed); they must not share a
    namespace."""
    r = seeded_client.get("/api/generators?page=1&per_page=10")
    assert r.status_code == 200
    body = r.json()
    assert body["page"] == 1
    assert body["per_page"] == 10


def test_api_v1_does_not_advertise_v2_routes(seeded_client) -> None:
    """OpenAPI grouping: v1 endpoints have no tag of `v2`. The v2
    namespace is reserved for breaking changes."""
    spec = seeded_client.get("/openapi.json").json()
    for path, ops in spec.get("paths", {}).items():
        if path.startswith("/api/v2/"):
            continue
        for verb, op in ops.items():
            if verb == "parameters":
                continue
            tags = op.get("tags", [])
            assert "v2" not in tags, (
                f"v1 endpoint {verb.upper()} {path} carries v2 tag"
            )
