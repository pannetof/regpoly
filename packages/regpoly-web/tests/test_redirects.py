# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 2 red — /library/papers preserved via 307 redirect.

S7.2 user-story URLs are pinned. The Library merge collapses
`templates/library/papers.html` into `/library?tab=papers`. The
redirect must:
  - Return 307 (preserves method)
  - Land on /library?tab=papers
  - Be declared in OpenAPI (via responses={307: ...}) so SDK
    consumers see the indirection.
"""

from __future__ import annotations


def test_library_papers_returns_307(seeded_client) -> None:
    r = seeded_client.get("/library/papers", follow_redirects=False)
    assert r.status_code == 307, (
        f"/library/papers must 307 → /library?tab=papers; got {r.status_code}"
    )
    assert "/library" in r.headers["location"]
    assert "tab=papers" in r.headers["location"]


def test_library_papers_redirects_to_tab_papers(seeded_client) -> None:
    r = seeded_client.get("/library/papers", follow_redirects=True)
    assert r.status_code == 200
    # Followed redirect ends up on the Library tab content.
    final = str(r.url)
    assert "/library" in final


def test_library_papers_redirect_documented_in_openapi(seeded_client) -> None:
    spec = seeded_client.get("/openapi.json").json()
    paths = spec.get("paths", {})
    op = paths.get("/library/papers", {}).get("get")
    assert op is not None, "/library/papers must remain in the OpenAPI spec"
    responses = op.get("responses", {})
    assert "307" in responses, (
        "307 response must be declared (responses={307: ...}) so "
        "SDK consumers see the indirection"
    )
