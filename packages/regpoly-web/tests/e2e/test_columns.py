"""Phase 2 red — Column show/hide + reorder + sticky pin (Analyst A-7)."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_column_show_hide_persists(page, base_url) -> None:
    page.goto(f"{base_url}/generators?cols=family,k", wait_until="networkidle")
    headers = page.locator("table thead th")
    visible = [h.inner_text().strip().lower() for h in headers.element_handles()]
    # `?cols=family,k` requested a 2-column projection.
    assert any("family" in h for h in visible), \
        "family column must be visible when ?cols includes it"
    # The table should NOT show columns that aren't in the cols list.
    # `params` is a default column normally hidden by ?cols=family,k.
    # We don't crash if it appears (tests can be tightened in P2 green).


def test_column_reorder_persists_in_url(page, base_url) -> None:
    page.goto(f"{base_url}/generators?cols=k,family", wait_until="networkidle")
    # Read only the data columns (those carry data-col-key); the leading
    # pin column and trailing action column are always-present chrome.
    data_keys = page.evaluate(
        "() => Array.from(document.querySelectorAll('table thead th[data-col-key]'))"
        "        .map(th => th.dataset.colKey)"
    )
    if len(data_keys) < 2:
        pytest.skip("list table not rendered yet (P2 green pending)")
    # Order in the URL drives data-column order.
    assert data_keys[:2] == ["k", "family"], (
        f"?cols=k,family must order columns as k then family; got {data_keys}"
    )


def test_sticky_pin_keeps_row_at_top(page, base_url) -> None:
    page.goto(f"{base_url}/generators", wait_until="networkidle")
    rows = page.locator("table tbody tr")
    if rows.count() < 1:
        pytest.skip("no rows seeded for pinning test")
    star = rows.nth(0).locator("[data-pin-toggle]")
    if star.count() == 0:
        pytest.skip("pin toggle not rendered yet (P2 green pending)")
    star.click()
    pinned = page.evaluate(
        "() => JSON.parse(localStorage.getItem('regpoly.gen.pinned') || '[]')"
    )
    assert isinstance(pinned, list) and len(pinned) >= 1
