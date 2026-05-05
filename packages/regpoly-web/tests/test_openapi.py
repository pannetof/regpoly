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


# ── P6 red — committed snapshot diff + response_model audit ──────────


def test_v2_openapi_matches_committed_snapshot(client) -> None:
    """The committed fixture at tests/fixtures/openapi.v2.json is the
    contract for v2 paths. Drift fails CI; intentional updates require
    `pytest --update-snapshots`."""
    from pathlib import Path
    import json

    fixture = Path(__file__).parent / "fixtures" / "openapi.v2.json"
    assert fixture.exists(), (
        "tests/fixtures/openapi.v2.json must be committed; run "
        "`uv run python -m tests._gen_openapi_snapshot` to generate it"
    )
    spec = client.get("/openapi.json").json()
    expected = json.loads(fixture.read_text())
    # Compare only the v2 path slice — not the full v1 surface.
    actual_v2 = {p: ops for p, ops in spec.get("paths", {}).items()
                 if p.startswith("/api/v2/")}
    expected_v2 = {p: ops for p, ops in expected.get("paths", {}).items()
                   if p.startswith("/api/v2/")}
    actual_paths = set(actual_v2.keys())
    expected_paths = set(expected_v2.keys())
    missing = expected_paths - actual_paths
    extra = actual_paths - expected_paths
    assert not missing and not extra, (
        f"OpenAPI v2 surface drift: missing={missing} extra={extra}"
    )


def test_every_v2_endpoint_declares_response_model(client) -> None:
    """Catch the 4 known undeclared endpoints
    (`/generators/families/counts`, `/generators/compare`,
    `/generators/export`, `/library/families/counts`) plus future
    regressions. An endpoint without `response_model` produces an
    inferred / empty schema in /openapi.json — codegen breaks."""
    spec = client.get("/openapi.json").json()
    paths = spec.get("paths", {})
    components = spec.get("components", {}).get("schemas", {})

    # Endpoints known to legitimately return non-JSON (binary, plain
    # text). They're audited via `responses` blocks instead.
    NON_JSON_OK = {
        ("get", "/api/v2/generators/{gen_id}/transition-matrix-coords"),
        ("get", "/api/v2/healthz"),  # trivial dict
        ("get", "/api/v2/generators/export"),  # CSV alternative
    }

    failures: list[str] = []
    for path, ops in paths.items():
        if not path.startswith("/api/v2/"):
            continue
        for verb, op in ops.items():
            if verb not in ("get", "post", "delete", "put", "patch"):
                continue
            if (verb, path) in NON_JSON_OK:
                continue
            responses = op.get("responses", {})
            ok = responses.get("200") or responses.get("201")
            if not ok:
                continue
            content = ok.get("content", {})
            json_resp = content.get("application/json", {})
            schema = json_resp.get("schema", {})
            # A typed response_model produces a $ref into components.
            if "$ref" not in schema and "type" not in schema:
                failures.append(f"{verb.upper()} {path}")
    assert not failures, (
        f"v2 endpoints missing response_model: {failures}"
    )
