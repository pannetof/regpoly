# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 1 red set — macro render tests for the v2 design system.

These tests render Jinja macros against the same environment FastAPI
uses, and assert HTML structure with BeautifulSoup. They fail until
P1 lands the macros under templates/macros/.
"""

from __future__ import annotations

import pytest
from bs4 import BeautifulSoup
from jinja2 import Environment, FileSystemLoader, select_autoescape


def _env() -> Environment:
    """Return a Jinja2 environment rooted at the package's templates dir.

    Mirrors what FastAPI's Jinja2Templates builds for the app.
    """
    from pathlib import Path

    import regpoly_web

    pkg_root = Path(regpoly_web.__file__).parent
    return Environment(
        loader=FileSystemLoader(str(pkg_root / "templates")),
        autoescape=select_autoescape(["html"]),
    )


def _render_macro(template: str, macro: str, *args, **kwargs) -> str:
    env = _env()
    tpl = env.get_template(template)
    return getattr(tpl.module, macro)(*args, **kwargs)


# --- metric_card -------------------------------------------------------


def test_metric_card_renders_value_and_title() -> None:
    html = str(_render_macro("macros/ui.html", "metric_card",
                             title="Period exponent", value="607"))
    soup = BeautifulSoup(html, "html.parser")
    card = soup.find(class_="card-sm") or soup.find(class_="card")
    assert card is not None, "metric_card should emit a Tabler card"
    assert "Period exponent" in card.get_text()
    assert "607" in card.get_text()


# --- status_badge ------------------------------------------------------


@pytest.mark.parametrize("status,expected_class", [
    ("pending", "bg-secondary-lt"),
    ("running", "bg-primary-lt"),
    ("paused", "bg-warning-lt"),
    ("completed", "bg-success-lt"),
    ("cancelled", "bg-orange-lt"),
    ("failed", "bg-danger-lt"),
    ("primitive_yes", "bg-success"),
    ("primitive_no", "bg-danger"),
])
def test_status_badge_for_each_status(status: str, expected_class: str) -> None:
    html = str(_render_macro("macros/ui.html", "status_badge", status))
    soup = BeautifulSoup(html, "html.parser")
    badge = soup.find(class_=lambda c: c and "badge" in c.split())
    assert badge is not None, f"status_badge({status!r}) should emit a .badge"
    classes = badge.get("class", [])
    assert expected_class in classes, (
        f"status_badge({status!r}) should carry {expected_class!r}, got {classes}"
    )


def test_status_badge_cancelled_carries_error_message_tooltip() -> None:
    html = str(_render_macro("macros/ui.html", "status_badge",
                             "cancelled", error_message="reaped at startup"))
    soup = BeautifulSoup(html, "html.parser")
    badge = soup.find(class_=lambda c: c and "badge" in c.split())
    assert badge is not None
    assert badge.get("title") == "reaped at startup", (
        "reaped runs surface error_message as a tooltip"
    )


# --- hex_value (first-8 + middle-12 + last-8 truncation) ---------------


def test_hex_value_short_value_not_truncated_in_table_context() -> None:
    short = "0xabcd1234"
    html = str(_render_macro("macros/ui.html", "hex_value", short, context="table"))
    soup = BeautifulSoup(html, "html.parser")
    span = soup.find(class_="hex-value") or soup.find("span")
    assert span is not None
    assert short in span.get_text()


def test_hex_value_truncates_first_middle_last_in_table_context() -> None:
    # 0x + 64 hex chars = 66 chars total — well above the 32-char threshold.
    long_value = "0x" + "1234abcd5678ef90" * 4  # 0x + 64 = 66 chars
    html = str(_render_macro("macros/ui.html", "hex_value",
                             long_value, context="table"))
    soup = BeautifulSoup(html, "html.parser")
    span = soup.find(class_="hex-value") or soup.find("span")
    assert span is not None
    text = span.get_text()
    # The discriminating middle must remain visible.
    assert "…" in text, "long hex truncates with an ellipsis"
    assert text.startswith("0x"), "first segment kept"
    # First-8 + middle-12 + last-8 total 28 visible hex chars + 0x prefix +
    # two ellipses ≈ 32 chars max.
    visible = text.replace("…", "")
    assert len(visible) <= 32, f"truncated hex too long: {len(visible)} chars"
    # full value preserved for copy
    assert span.get("data-hex") == long_value
    assert span.get("title") == long_value


def test_hex_value_full_in_detail_context() -> None:
    long_value = "0x" + "1234abcd5678ef90" * 4
    html = str(_render_macro("macros/ui.html", "hex_value",
                             long_value, context="detail"))
    soup = BeautifulSoup(html, "html.parser")
    span = soup.find(class_="hex-value") or soup.find("span")
    assert span is not None
    text = span.get_text()
    assert "…" not in text, "detail context renders full value"
    assert long_value in text


def test_hex_value_marks_copy_target() -> None:
    html = str(_render_macro("macros/ui.html", "hex_value",
                             "0xabcd", context="table"))
    soup = BeautifulSoup(html, "html.parser")
    span = soup.find(class_="hex-value") or soup.find("span")
    assert span is not None
    # Click delegate in regpoly.js looks for [data-copy].
    assert span.has_attr("data-copy"), "hex_value should carry data-copy"


# --- char_poly_card ----------------------------------------------------


def test_char_poly_card_emits_three_tabs() -> None:
    html = str(_render_macro("macros/poly.html", "char_poly_card",
                             "0xdeadbeef", 32, 16))
    soup = BeautifulSoup(html, "html.parser")
    # Tabler tabs are <ul class="nav nav-tabs"> with <li class="nav-item">…
    tabs = soup.select(".nav-tabs .nav-item, .nav-tabs .nav-link")
    assert len(tabs) >= 3, (
        "char_poly_card should expose three tabs: Exponents, x-notation, Hex"
    )
    text = soup.get_text().lower()
    assert "exponent" in text
    assert "hex" in text


def test_char_poly_card_emits_katex_delimiters_in_exponents_tab() -> None:
    html = str(_render_macro("macros/poly.html", "char_poly_card",
                             "0xdeadbeef", 32, 16))
    # Exponents and x-notation tabs use \( ... \) so KaTeX renderMathInElement
    # picks them up. Hex tab is plain text.
    assert "\\(" in html or "$$" in html, (
        "char_poly_card must emit KaTeX inline delimiters for math tabs"
    )
