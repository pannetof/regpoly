"""Phase 1+ — OpenAPI fixture + path presence.

Phase 1 only asserts the OpenAPI spec is reachable. The /api/v2/
namespace and the committed fixture diff arrive in P2+.
"""

from __future__ import annotations


def test_openapi_reachable(client) -> None:
    """FastAPI auto-emits /openapi.json. Just confirm it's served."""
    r = client.get("/openapi.json")
    assert r.status_code == 200
    body = r.json()
    assert "paths" in body
    assert "openapi" in body


def test_openapi_v2_namespace_reserved(client) -> None:
    """P1 reserves the /api/v2/ namespace by mounting an empty router.
    Asserting the prefix is present in the spec ensures P2 can add
    endpoints without re-plumbing the app."""
    r = client.get("/openapi.json")
    body = r.json()
    paths = body.get("paths", {})
    # Phase 1 may not have any v2 endpoints yet, but the prefix should
    # be reachable (e.g. via a stub /api/v2/healthz route or merely a
    # mounted router). We accept either: at least one /api/v2/ path
    # registered, OR the app exposes the v2 router on app.state.
    has_v2_path = any(p.startswith("/api/v2/") for p in paths)
    has_v2_router = hasattr(client.app.state, "v2_router_registered")
    assert has_v2_path or has_v2_router, (
        "v2 namespace must be reserved at P1 so P2 endpoints attach cleanly"
    )
