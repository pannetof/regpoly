# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — command palette (⌘K)."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_meta_k_opens_command_palette(page, base_url) -> None:
    page.goto(base_url, wait_until="networkidle")
    page.keyboard.press("Meta+K")
    palette = page.locator("[data-cmdk-root], [data-command-palette]")
    assert palette.count() >= 1, "⌘K must open the command palette"
    assert palette.first.is_visible(), "command palette must be visible"


def test_ctrl_k_opens_on_non_mac(page, base_url) -> None:
    page.goto(base_url, wait_until="networkidle")
    page.keyboard.press("Control+K")
    palette = page.locator("[data-cmdk-root], [data-command-palette]")
    assert palette.count() >= 1
