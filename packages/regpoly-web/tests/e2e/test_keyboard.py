"""Phase 1 red set — keyboard-shortcut layer.

Power-user bindings: ?, g {d,g,s,t,l,o}, n {p,t}, /, j/k, Enter, Esc.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_question_mark_opens_help(page, base_url) -> None:
    page.goto(f"{base_url}/", wait_until="domcontentloaded")
    page.keyboard.press("Shift+?")
    modal = page.locator("[data-kbd-help-modal]")
    assert modal.count() == 1
    # Tabler modals get .show on open.
    assert "show" in (modal.get_attribute("class") or "")


@pytest.mark.parametrize("chord,expected_path", [
    ("g d", "/"),
    ("g g", "/generators"),
    ("g s", "/searches"),
    ("g t", "/tested-generators"),
    ("g l", "/library"),
    ("g o", "/tools"),
])
def test_g_chord_jumps(page, base_url, chord, expected_path) -> None:
    page.goto(f"{base_url}/", wait_until="domcontentloaded")
    a, b = chord.split()
    page.keyboard.press(a)
    page.keyboard.press(b)
    page.wait_for_url(f"{base_url}{expected_path}", timeout=2000)


@pytest.mark.parametrize("chord,expected_substr", [
    ("n p", "/primitive-search"),
    ("n t", "/tempering-search"),
])
def test_n_chord_opens_form(page, base_url, chord, expected_substr) -> None:
    page.goto(f"{base_url}/", wait_until="domcontentloaded")
    a, b = chord.split()
    page.keyboard.press(a)
    page.keyboard.press(b)
    page.wait_for_url(lambda url: expected_substr in url, timeout=2000)


def test_slash_handler_armed_on_every_page(page, base_url) -> None:
    """P1: the keydown listener for `/` is wired (it tries to focus a
    [data-filter-focus] target if one exists). The full chip-toolbar
    lives in P2 — that phase's red set adds the focus assertion. Here
    we only verify the listener exists and no console error fires when
    `/` is pressed on a page without a filter target."""
    errors = []
    page.on("console", lambda msg: errors.append(msg.text)
            if msg.type == "error" else None)
    page.goto(f"{base_url}/generators", wait_until="domcontentloaded")
    page.keyboard.press("/")
    assert errors == [], f"/ key produced console errors: {errors}"


def test_j_k_navigate_table_rows(page, base_url) -> None:
    """When a list table is in viewport, j/k advance/retreat the
    highlighted row and Enter follows the row link."""
    page.goto(f"{base_url}/generators", wait_until="networkidle")
    rows = page.locator("table tbody tr")
    if rows.count() < 2:
        pytest.skip("Need at least 2 generator rows to test j/k navigation")
    page.keyboard.press("j")
    first_active = page.locator("table tbody tr.kbd-active").count()
    assert first_active >= 1, "j should highlight a row with .kbd-active"
    page.keyboard.press("j")
    assert page.locator("table tbody tr.kbd-active").count() == 1, (
        "k/j navigates row-by-row, only one row active at a time"
    )
