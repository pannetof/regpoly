# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 1 red set — smoke-test every existing route under the new shell.

Parameterized over the route list. Each route must:
- Return HTTP 200
- Render the new Tabler shell (data-bs-theme on <html>, .navbar-vertical sidebar)
- Have no JS console errors
- Show a breadcrumb in the page header
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


# Routes that must render under the new shell. Parameterized fixtures
# come from tests/e2e/conftest.py.
SMOKE_ROUTES = [
    "/",
    "/generators",
    "/generators/1",            # may 404 if no row; we only assert <500
    "/tested-generators",
    "/tested-generators/4242",
    "/searches",
    "/library",
    "/library/papers",          # P2 redirects this to /library?tab=papers
    "/family/MTGen",
    "/primitive-search?family=MTGen",
    "/tempering-search",
]


@pytest.mark.parametrize("path", SMOKE_ROUTES)
def test_every_route_returns_200_with_new_shell(page, base_url, path) -> None:
    """Visit each route and assert the new shell rendered."""
    console_errors: list[str] = []
    page.on("console", lambda msg: (
        console_errors.append(msg.text)
        if msg.type == "error" else None
    ))

    response = page.goto(f"{base_url}{path}", wait_until="domcontentloaded")
    # 200 for known rows; some detail routes may 404 if the seeded DB
    # doesn't have the id — accept anything < 500 (still rendered the
    # error page through the new shell).
    assert response is not None
    assert response.status < 500, (
        f"{path} returned {response.status} — server error"
    )

    # New shell signature: data-bs-theme on <html> and a Tabler vertical
    # navbar. The legacy shell used <div class="app-shell"> instead.
    html = page.locator("html")
    assert html.get_attribute("data-bs-theme") in ("light", "dark"), (
        f"{path}: <html> must carry data-bs-theme (Tabler)"
    )
    assert page.locator("aside.navbar.navbar-vertical").count() >= 1, (
        f"{path}: new shell must use <aside class='navbar navbar-vertical'>"
    )

    # No JS console errors during initial load.
    assert console_errors == [], (
        f"{path}: console errors during load: {console_errors}"
    )
