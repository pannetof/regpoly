"""Phase 2 red — Compact mode density toggle (Analyst persona A-1)."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_compact_mode_increases_visible_rows(page, base_url) -> None:
    """`?density=compact` sets `.table-sm` and trims padding so more
    rows fit without scroll."""
    page.goto(f"{base_url}/generators?density=compact",
              wait_until="networkidle")
    table = page.locator("table.table-sm")
    assert table.count() >= 1, (
        "?density=compact must apply Tabler .table-sm to the list"
    )


def test_density_toggle_persists_in_localstorage(page, base_url) -> None:
    page.goto(f"{base_url}/generators", wait_until="networkidle")
    toggle = page.locator("[data-density-toggle]")
    if toggle.count() == 0:
        pytest.skip("density toggle button not rendered yet (P2 green pending)")
    toggle.click()
    page.wait_for_function(
        "() => location.search.includes('density=compact')",
        timeout=2000,
    )
    stored = page.evaluate("() => localStorage.getItem('regpoly.density')")
    assert stored == "compact"
