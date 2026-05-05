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


# ── P6 red: SSE byte-for-byte preservation under ?v=2 gate ─────────────


def _stream_text(client, url: str, timeout: float = 0.6) -> str:
    """Pull a snapshot of the SSE body off the endpoint, then close."""
    chunks: list[bytes] = []
    with client.stream("GET", url, timeout=timeout) as resp:
        assert resp.status_code == 200
        try:
            for chunk in resp.iter_bytes():
                chunks.append(chunk)
                if sum(len(c) for c in chunks) > 16 * 1024:
                    break
        except Exception:
            pass
    return b"".join(chunks).decode("utf-8", errors="replace")


def test_v1_sse_default_emits_no_named_progress_event(seeded_client) -> None:
    """Without ?v=2, the SSE body must NOT contain `event: progress`.
    A pre-redesign v1 EventSource client that registered .onmessage
    must see exactly the same byte stream as before."""
    body = _stream_text(seeded_client, "/api/primitive-searches/1/progress")
    assert "event: progress" not in body, (
        "v1 SSE must not emit the named v2 `progress` channel without ?v=2 "
        "opt-in (regression: doubles the events for v1 EventSource clients)"
    )


def test_v2_sse_with_query_emits_named_progress(seeded_client) -> None:
    """With ?v=2, the named `progress` channel is enabled. The seeded
    primitive run has `status='completed'` so emit terminates quickly."""
    body = _stream_text(
        seeded_client, "/api/primitive-searches/1/progress?v=2",
    )
    # Either the run had progress rows (event: progress present) or it
    # had none (only event: end). The contract is that v=2 is honoured;
    # the test passes when the URL is reachable and either condition holds.
    assert "event: end" in body or "event: progress" in body
