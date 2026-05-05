"""Phase 2 red — Quick-actions card on dashboard (drops from spec)."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_quick_actions_card_present(page, base_url) -> None:
    page.goto(f"{base_url}/", wait_until="networkidle")
    card = page.locator("[data-quick-actions]")
    assert card.count() >= 1, (
        "dashboard must include a [data-quick-actions] card "
        "(per regpoly-web design-spec)"
    )
    # Three primary entry points, per the locked decision.
    expected_actions = ["new primitive search",
                        "new tempering search",
                        "browse generators"]
    text = card.inner_text().lower()
    for action in expected_actions:
        assert action in text, f"quick-actions card missing {action!r}"
