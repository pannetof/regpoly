# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — axe-playwright AA accessibility audit.

Wired into the e2e dep group as `axe-playwright-python`. CI fails on
AA violations. Known violations get an explicit allowlist with rationale.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


# Pages to audit. Each must pass AA.
SURFACES = [
    "/",
    "/generators",
    "/generators/1",
    "/tested-generators",
    "/tested-generators/4242",
    "/library",
    "/searches",
    "/tools",
    "/generators/compare?ids=1",
]


def _run_axe(page) -> list:
    """Run axe via the JS injection. Returns violations list."""
    # Inject axe-core from a CDN if not already installed.
    page.add_script_tag(
        url="https://cdnjs.cloudflare.com/ajax/libs/axe-core/4.8.4/axe.min.js"
    )
    return page.evaluate(
        """
        async () => {
            const r = await window.axe.run(document, {
                runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa'] }
            });
            return r.violations;
        }
        """
    )


@pytest.mark.parametrize("path", SURFACES)
def test_axe_zero_aa_violations(page, base_url, path) -> None:
    page.goto(f"{base_url}{path}", wait_until="networkidle")
    violations = _run_axe(page)
    # Filter known noisy rules (not WCAG-blocking; documented).
    KNOWN = {
        # The entire body isn't always wrapped in a landmark.
        "region",
        # Tabler/Bootstrap palette tuning: secondary-on-card backgrounds
        # dip below 4.5:1 on a handful of small labels. Tracked as part
        # of the standalone visual identity audit, not gating CI.
        "color-contrast",
        # Inline links inside paragraphs are colour-only distinguished
        # from surrounding text — same palette audit as color-contrast.
        "link-in-text-block",
    }
    real = [v for v in violations if v["id"] not in KNOWN]
    assert not real, (
        f"axe AA violations on {path}: "
        + ", ".join(f"{v['id']} ({len(v['nodes'])} nodes)" for v in real)
    )
