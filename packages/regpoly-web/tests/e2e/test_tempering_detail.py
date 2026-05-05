"""Phase 6 red — tempering-detail page Tabler rewrite (Researcher #1).

After P6 the page must:
  - have a header strip with Pause/Resume/Cancel/Restart/Duplicate buttons
  - show 4 metric cards (Combos / Rate / Best / ETA)
  - render an SVG/canvas best_se sparkline seeded from history on reload
  - use Tabler card vocabulary, NOT Pico (<hgroup>, <article>, .status-pill)
"""

from __future__ import annotations

import json
import sqlite3

import pytest

pytestmark = pytest.mark.e2e


@pytest.fixture
def seeded_tempering_run(seeded_db) -> int:
    """Seed a tempering_search_run we can navigate to for the detail
    page. Returns the run id."""
    conn = sqlite3.connect(seeded_db)
    try:
        cur = conn.execute(
            "INSERT INTO tempering_search_run "
            "(test_type, test_config, Lmax, nb_tries, status, "
            " combos_total, combos_done, best_se, elapsed_seconds) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ("equidistribution", json.dumps({"L": 32}), 32, 1,
             "running", 50, 12, 7, 30.0),
        )
        run_id = cur.lastrowid
        # one component
        conn.execute(
            "INSERT INTO tempering_search_component "
            "(search_run_id, component_index, tempering_config) "
            "VALUES (?, 0, ?)",
            (run_id, json.dumps([])),
        )
        conn.commit()
        return run_id
    finally:
        conn.close()


def test_pause_and_resume_buttons_present(page, base_url, seeded_tempering_run) -> None:
    page.goto(f"{base_url}/tempering-search/{seeded_tempering_run}",
              wait_until="networkidle")
    # Buttons may be `<button>` or `<a>` with data-action attributes.
    pause = page.locator(
        "[data-action='pause'], button:has-text('Pause')"
    )
    resume = page.locator(
        "[data-action='resume'], button:has-text('Resume')"
    )
    assert pause.count() >= 1, "Pause button missing on tempering detail"
    assert resume.count() >= 1, "Resume button missing on tempering detail"


def test_metric_cards_show_combos_rate_best_eta(
    page, base_url, seeded_tempering_run,
) -> None:
    page.goto(f"{base_url}/tempering-search/{seeded_tempering_run}",
              wait_until="networkidle")
    cards = page.locator(".card-sm, .metric-card")
    text = " ".join(c.inner_text() for c in cards.element_handles())
    for label in ("Combos", "Rate", "Best", "ETA"):
        assert label in text, (
            f"tempering detail missing metric card '{label}'"
        )


def test_no_pico_vocabulary(page, base_url, seeded_tempering_run) -> None:
    """Pico CSS holdovers (<hgroup>, <article>, .status-pill, table.compact)
    must be gone after the rewrite."""
    page.goto(f"{base_url}/tempering-search/{seeded_tempering_run}",
              wait_until="networkidle")
    pico_count = page.evaluate(
        "() => document.querySelectorAll("
        "  'hgroup, article, .status-pill, table.compact'"
        ").length"
    )
    assert pico_count == 0, (
        f"tempering detail still uses Pico vocabulary: {pico_count} elements"
    )
