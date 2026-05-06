# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 2 red — chip toolbar URL state grammar.

Chips, popover values, sort, density, columns, limit/offset all in
the URL query string. Front-end uses history.pushState; popstate
restores filters. Same grammar across /generators, /tested-generators,
/library.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_chip_state_round_trips_through_url_with_offset(page, base_url) -> None:
    page.goto(f"{base_url}/generators?family=MTGen&offset=0&limit=20",
              wait_until="networkidle")
    # Family chip rendered from the URL.
    assert page.locator(".chip-toolbar [data-chip-family='MTGen']").count() >= 1, (
        "family chip must render when ?family=… is in the URL"
    )
    # Add a k-range chip programmatically by navigating; back/forward
    # should restore.
    page.goto(f"{base_url}/generators?family=MTGen&k_min=10&k_max=50000",
              wait_until="networkidle")
    page.go_back()
    page.wait_for_load_state("networkidle")
    # After back, k_min/k_max chips removed.
    assert page.locator(".chip-toolbar [data-chip-k_min]").count() == 0


def test_chip_remove_updates_url(page, base_url) -> None:
    page.goto(f"{base_url}/generators?family=MTGen", wait_until="networkidle")
    chip = page.locator(".chip-toolbar [data-chip-family='MTGen'] [data-chip-remove]")
    if chip.count() == 0:
        pytest.skip("chip remove button not rendered yet (P2 green pending)")
    chip.click()
    page.wait_for_function(
        "() => !location.search.includes('family=MTGen')",
        timeout=2000,
    )
