"""Phase 6 red — list pages default to compact + sticky thead."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_generators_default_density_is_compact(page, base_url) -> None:
    page.goto(f"{base_url}/generators", wait_until="networkidle")
    table = page.locator("table.table-sm")
    assert table.count() >= 1, (
        "default density must be compact (table-sm); user prefers density"
    )


def test_generators_thead_is_sticky(page, base_url) -> None:
    page.goto(f"{base_url}/generators", wait_until="networkidle")
    pos = page.evaluate(
        "() => {"
        "  const th = document.querySelector('table thead th');"
        "  if (!th) return null;"
        "  return getComputedStyle(th).position;"
        "}"
    )
    assert pos in ("sticky", "-webkit-sticky"), (
        f"thead must be position: sticky; got {pos!r}"
    )


def test_tested_generators_thead_is_sticky(page, base_url) -> None:
    page.goto(f"{base_url}/tested-generators", wait_until="networkidle")
    pos = page.evaluate(
        "() => {"
        "  const th = document.querySelector('table thead th');"
        "  if (!th) return null;"
        "  return getComputedStyle(th).position;"
        "}"
    )
    assert pos in ("sticky", "-webkit-sticky"), (
        f"thead must be position: sticky; got {pos!r}"
    )
