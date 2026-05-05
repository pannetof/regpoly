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


def test_katex_runtime_available_on_every_page(page, base_url) -> None:
    """P1 promise: KaTeX is loaded AND `renderMathInElement` is invoked
    on DOMContentLoaded. The full integration with the char_poly_card
    macro lands in P4 (test moves to tests/e2e/test_generator_detail.py
    then). Here we only verify the runtime is available."""
    page.goto(f"{base_url}/", wait_until="networkidle")
    runtime = page.evaluate(
        "() => typeof window.renderMathInElement === 'function'"
    )
    assert runtime, "renderMathInElement must be loaded on every page"


def test_no_js_console_errors_on_initial_load(page, base_url) -> None:
    """Loading the dashboard must not produce any console errors."""
    errors: list[str] = []
    page.on("console", lambda msg: (
        errors.append(f"{msg.type}:{msg.text}") if msg.type == "error" else None
    ))
    page.goto(f"{base_url}/", wait_until="networkidle")
    assert errors == [], f"console errors: {errors}"
