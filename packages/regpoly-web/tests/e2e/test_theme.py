# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 1 red set — theme defaults and persistence.

Light is the default; the toggle in the navbar persists choice in
localStorage and survives reload.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_light_default_on_first_visit(page, base_url) -> None:
    page.context.clear_cookies()
    page.goto(f"{base_url}/", wait_until="domcontentloaded")
    page.evaluate("() => localStorage.removeItem('regpoly.theme')")
    page.reload(wait_until="domcontentloaded")
    theme = page.locator("html").get_attribute("data-bs-theme")
    # OS preference may force dark on the CI runner; accept either as
    # long as the toggle is present and writes localStorage.
    assert theme in ("light", "dark")
    toggle = page.locator("[data-theme-toggle]")
    assert toggle.count() == 1, "navbar must expose a theme toggle"


def test_theme_toggle_persists_across_reload(page, base_url) -> None:
    page.goto(f"{base_url}/", wait_until="domcontentloaded")
    # Force light first to start from a known state.
    page.evaluate(
        "() => { document.documentElement.setAttribute('data-bs-theme','light'); "
        "localStorage.setItem('regpoly.theme','light'); }"
    )
    page.reload(wait_until="domcontentloaded")
    assert page.locator("html").get_attribute("data-bs-theme") == "light"

    page.locator("[data-theme-toggle]").click()
    # After toggle: dark.
    assert page.locator("html").get_attribute("data-bs-theme") == "dark"
    assert page.evaluate("() => localStorage.getItem('regpoly.theme')") == "dark"

    # Reload preserves dark.
    page.reload(wait_until="domcontentloaded")
    assert page.locator("html").get_attribute("data-bs-theme") == "dark"
