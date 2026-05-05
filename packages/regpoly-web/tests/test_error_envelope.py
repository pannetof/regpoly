"""Phase 2 red — JSON error envelope contract.

Per Developer persona D-7: standardise on FastAPI's default
`{detail: str | object}` shape for both HTTP errors and validation
errors. SDK consumers can rely on it.
"""

from __future__ import annotations


def test_404_returns_detail_field(seeded_client) -> None:
    """Unknown family endpoint returns 404 with {detail: str}."""
    r = seeded_client.get("/api/families/NoSuchFamily")
    assert r.status_code == 404
    body = r.json()
    assert "detail" in body, "404 must use the {detail} envelope"
    assert isinstance(body["detail"], str)


def test_422_validation_error_returns_detail_field(seeded_client) -> None:
    """A bad query-param value triggers Pydantic validation; the
    response carries `detail` (an array of validator failures by
    default in FastAPI)."""
    # /api/generators expects per_page<=500; 9999 trips the validator.
    r = seeded_client.get("/api/generators?per_page=9999")
    assert r.status_code == 422
    body = r.json()
    assert "detail" in body
    # FastAPI's default 422 detail is a list of error dicts.
    assert isinstance(body["detail"], (list, str))


def test_v2_404_uses_same_envelope(seeded_client) -> None:
    """v2 endpoints honour the same envelope — no RFC 7807."""
    r = seeded_client.get("/api/v2/dashboard/summary/no-such-suffix")
    # 404 (or 405) — either returns `{detail}` per FastAPI default.
    assert r.status_code in (404, 405)
    body = r.json()
    assert "detail" in body
