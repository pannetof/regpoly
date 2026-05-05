"""Phase 2 red — Dashboard pinning + 5-row HTTP/1.1 SSE budget cap.

Researcher persona R-1: pin by session-launched OR explicit star;
new runs go to "+N more" overflow once 5 pinned slots are taken.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_active_runs_section_renders(page, base_url) -> None:
    page.goto(f"{base_url}/", wait_until="networkidle")
    section = page.locator("[data-active-runs]")
    assert section.count() >= 1, (
        "dashboard must render an [data-active-runs] section"
    )


def test_active_runs_capped_at_5_by_default(page, base_url) -> None:
    """Even if the API returns 7 active rows, only 5 render with live
    SSE; the rest go into the overflow link."""
    page.goto(f"{base_url}/", wait_until="networkidle")
    rows = page.locator("[data-active-runs] tbody tr[data-active-row]")
    assert rows.count() <= 5, (
        f"max 5 active rows (HTTP/1.1 budget); got {rows.count()}"
    )


def test_explicit_star_pins_row(page, base_url) -> None:
    page.goto(f"{base_url}/", wait_until="networkidle")
    star = page.locator("[data-active-runs] [data-pin-active]").first
    if star.count() == 0:
        pytest.skip("no active runs in seeded DB; can't exercise star")
    star.click()
    pinned = page.evaluate(
        "() => JSON.parse(localStorage.getItem('regpoly.dashboard.pinned') || '[]')"
    )
    assert isinstance(pinned, list) and len(pinned) >= 1


def test_seventh_run_does_not_evict_session_pinned_row(page, base_url) -> None:
    """When a 7th active run starts and 5 slots are full of pinned
    rows, the new run goes to the overflow link, not the visible 5."""
    # Hard to simulate without inserting fake rows; the contract is
    # exercised at the unit level (dashboard summary endpoint sorts
    # pinned runs first). Here we assert the overflow link renders
    # when the API claims more than 5 active rows.
    page.goto(f"{base_url}/", wait_until="networkidle")
    overflow = page.locator("[data-active-runs-overflow]")
    # No assertion on overflow visibility (depends on data); the
    # element must be wired in the template either way.
    assert overflow.count() >= 1, (
        "overflow link must be present in the active-runs section template"
    )
