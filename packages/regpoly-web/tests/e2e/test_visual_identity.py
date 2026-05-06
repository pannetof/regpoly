# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — visual identity verification.

After P6.UI lands the dashboard should pass these structural checks
even before pixel-perfect screenshot baselines are taken.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_sidebar_has_svg_logo_mark(page, base_url) -> None:
    page.goto(base_url, wait_until="networkidle")
    # The mark should be an inline <svg> inside the navbar, with class
    # `rp-mark` or carrying the brand fill.
    mark = page.locator("aside.navbar-vertical svg.rp-mark, aside.navbar-vertical [data-logo-mark]")
    assert mark.count() >= 1, (
        "sidebar logo mark missing; nav.html must inline an SVG mark"
    )


def test_brand_color_is_indigo_violet(page, base_url) -> None:
    page.goto(base_url, wait_until="networkidle")
    # The CSS variable --tblr-primary should be overridden to the
    # locked-decision palette (#5B5BD6 indigo-violet).
    primary = page.evaluate(
        "() => getComputedStyle(document.documentElement)"
        ".getPropertyValue('--tblr-primary').trim()"
    )
    assert primary in ("#5B5BD6", "#5b5bd6", "rgb(91, 91, 214)"), (
        f"brand primary not overridden, got {primary!r}"
    )


def test_tabler_icon_sprite_is_injected(page, base_url) -> None:
    page.goto(base_url, wait_until="networkidle")
    sprite = page.locator("svg[data-tabler-sprite] symbol")
    assert sprite.count() >= 5, (
        "Tabler Icons sprite must be injected so <use href='#tabler-…'> resolves"
    )


def test_inter_geist_jetbrains_loaded(page, base_url) -> None:
    page.goto(base_url, wait_until="networkidle")
    families = page.evaluate(
        "() => [...document.fonts.values()].map(f => f.family).join('|')"
    )
    # At least one of each face must be present.
    has_inter = "Inter" in families or "InterVariable" in families
    has_geist = "Geist" in families or "Manrope" in families
    has_mono = "JetBrains" in families or "Fira Code" in families
    assert has_inter and has_geist and has_mono, (
        f"typography stack not loaded: {families}"
    )


def test_no_unicode_glyphs_in_buttons_and_chrome(page, base_url) -> None:
    """The visual sweep replaces ★/☆/▾/×/→/✓/📄/☀/☾ with Tabler icons.
    Allow them to remain inside data-attributes (titles), but ban them
    from rendered button text where SVG icons should live."""
    page.goto(base_url, wait_until="networkidle")
    text = page.evaluate(
        "() => Array.from(document.querySelectorAll('button, a.btn, a.nav-link'))"
        "        .map(b => b.textContent.trim()).join('|')"
    )
    # We tolerate `→` inside the "<count> more →" overflow link; ban
    # the rest from button surfaces.
    banned = ["★", "☆", "▾", "✓", "☀", "☾"]
    found = [g for g in banned if g in text]
    assert not found, (
        f"Unicode glyphs leaked into button chrome: {found} "
        "(should be SVG Tabler icons)"
    )
