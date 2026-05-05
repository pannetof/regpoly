"""Phase 6 red — real chart system (uPlot) replaces bare polylines."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_dashboard_loads_uplot(page, base_url) -> None:
    page.goto(base_url, wait_until="networkidle")
    has_uplot = page.evaluate("() => typeof window.uPlot === 'function'")
    assert has_uplot, "uPlot must be loaded for charts"


def test_equidist_chart_renders_full_l_axis(page, base_url, seeded_db) -> None:
    page.goto(f"{base_url}/generators/1", wait_until="networkidle")
    chart = page.locator("[data-equidist-chart], svg.equidist-chart, canvas.equidist-chart")
    assert chart.count() >= 1, (
        "generator detail must render an equidist_chart, not a sparse table"
    )


def test_equidist_chart_has_no_redundant_table(page, base_url, seeded_db) -> None:
    """Plan said: chart only, no redundant sparse-values table beneath."""
    page.goto(f"{base_url}/generators/1", wait_until="networkidle")
    pis_tables = page.locator(
        "table:has(th:has-text('v')):has(th:has-text('gap'))"
    )
    assert pis_tables.count() == 0, (
        "old sparse 'v / ⌊k/v⌋ / t_v / gap' table must be removed; "
        "chart is the source of truth"
    )
