# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Tries-card progress display on the primitive-search detail page.

In exhaustive mode the search visits every combination exactly once,
so the natural denominator is the enumerator's total (already stored
on the run row as ``enum_total``).  In random mode we keep the legacy
``max_tries`` denominator.  Both render through the same Alpine
getter (``triesDenominatorFmt``) and a shared ``_formatBigDec``
helper that switches to scientific notation past 12 digits — type-2
w=32 totals (≈ 6.3·10^10) would otherwise render as a 12-digit blob.
"""

from __future__ import annotations


def test_detail_page_exposes_tries_denominator_helpers(client) -> None:
    """The Alpine helpers + template wiring must be present on the
    rendered detail page.  We render the page for an arbitrary id
    (the page itself is always rendered; data is fetched client-side
    via /api/primitive-searches/{id}), and look for the symbols."""
    r = client.get("/primitive-search/1")
    assert r.status_code == 200
    body = r.text

    # ── Helpers wired in ────────────────────────────────────────────
    for symbol in (
        "triesDenominator",
        "triesDenominatorFmt",
        "triesDoneFmt",
        "_triesDenomNum",
        "_formatBigDec",
    ):
        assert symbol in body, f"missing JS symbol: {symbol}"

    # ── Tries card uses the new getters, not the old direct refs ───
    # The new template renders `triesDoneFmt` and `triesDenominatorFmt`
    # — without these the card would still hard-wire run.max_tries
    # and never show a denominator in exhaustive mode.
    assert 'x-text="triesDoneFmt"' in body
    assert 'x-text="triesDenominatorFmt"' in body

    # ── ETA fallback uses the unified denominator too ──────────────
    # Otherwise exhaustive runs would always show "—" for ETA
    # because run.max_tries is null in that mode.
    assert "_triesDenomNum()" in body


def test_api_run_payload_carries_enum_total_for_exhaustive(client) -> None:
    """The /api/primitive-searches/{id} payload must surface
    ``enum_total`` so the detail page can render it as the
    denominator.  Drive an exhaustive Type-1 search end-to-end
    (small space: 4·31^3 = 119,164) and verify the row reports it."""
    create_body = {
        "family": "MarsaXorshiftGen",
        "structural_params": {"type": 1, "w": 32, "r": 1},
        "fixed_params": {},
        "max_tries": None,
        "max_seconds": None,
        "max_cost": None,
        "max_se": None,
        "search_mode": "exhaustive",
        "confirm_huge": False,
    }
    r = client.post("/api/primitive-searches", json=create_body)
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]

    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    # enum_total is a string in the API (decimal — can exceed 2^63
    # for larger configs); the detail page formats it via
    # _formatBigDec.
    assert detail["search_mode"] == "exhaustive"
    assert detail["enum_total"] == str(4 * 31 ** 3)
