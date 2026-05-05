"""Phase 1 red set — breadcrumb + KaTeX + no-console-errors across the shell."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


BREADCRUMB_ROUTES = [
    "/",
    "/generators",
    "/tested-generators",
    "/searches",
    "/library",
    "/family/MTGen",
    "/tempering-search",
]


@pytest.mark.parametrize("path", BREADCRUMB_ROUTES)
def test_breadcrumb_visible_on_every_route(page, base_url, path) -> None:
    page.goto(f"{base_url}{path}", wait_until="domcontentloaded")
    crumb = page.locator(".page-header .breadcrumb, ol.breadcrumb")
    assert crumb.count() >= 1, (
        f"{path}: breadcrumb missing — page header must include a Tabler breadcrumb"
    )


def test_katex_renders_char_poly_on_generators_detail(page, base_url) -> None:
    """The seeded primitive_generator (id=1) has char_poly=0xdeadbeef.
    The detail page renders char_poly_card via KaTeX. Asserting at least
    one .katex element appears proves renderMathInElement was invoked."""
    page.goto(f"{base_url}/generators/1", wait_until="networkidle")
    katex = page.locator(".katex")
    assert katex.count() >= 1, (
        "char_poly_card must render at least one KaTeX element on generator detail"
    )


def test_no_js_console_errors_on_initial_load(page, base_url) -> None:
    """Loading the dashboard must not produce any console errors."""
    errors: list[str] = []
    page.on("console", lambda msg: (
        errors.append(f"{msg.type}:{msg.text}") if msg.type == "error" else None
    ))
    page.goto(f"{base_url}/", wait_until="networkidle")
    assert errors == [], f"console errors: {errors}"
