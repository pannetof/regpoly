# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Smoke test for the MarsaXorshift-specific UI affordances on the
primitive-search form.

The form template lives at
``packages/regpoly-web/src/regpoly_web/templates/primitive_search/form.html``.
It auto-renders structural / non-structural inputs from
``param_specs()`` for every family.  MarsaXorshiftGen additionally:

  * shows a `<select>` dropdown for the `type` structural parameter
    (hard-coded to the 5 valid runtime values), and
  * grays out parameters that the selected `type` doesn't consume
    (per ``MARSA_TYPE_RELEVANT`` in the form's Alpine component).

This test renders the form once and asserts the relevant JS hooks and
HTML strings are present, so a refactor that accidentally drops them
fails loudly.
"""

from __future__ import annotations


def test_form_includes_marsa_type_dropdown_and_disable_logic(client) -> None:
    r = client.get("/primitive-search?family=MarsaXorshiftGen")
    assert r.status_code == 200
    body = r.text

    # ── Alpine helpers wired in ───────────────────────────────────────
    for symbol in (
        "isMarsaXorshiftFamily",
        "MARSA_TYPE_RELEVANT",
        "MARSA_TYPE_ENUMERATED",
        "onMarsaTypeChange",
        "paramDisabledReason",
        "marsaParamIrrelevant_",
        "marsaParamEnumerated_",
    ):
        assert symbol in body, f"missing JS symbol: {symbol}"

    # ── Dropdown placeholder + every runtime type option present ─────
    assert "Choose a type" in body
    for option_label in (
        ">1 — Type I",
        ">2 — Type II",
        ">3 — Marsaglia legacy multi-tap",
        ">4 — Brent four-xorshift",
        ">100 — General multi-tap",
    ):
        assert option_label in body, f"missing dropdown option: {option_label!r}"

    # ── Structural inputs now carry the disabled binding ─────────────
    # The MarsaXorshift change wired :disabled on the structural-card
    # input element — without it, the gray-out is purely cosmetic.
    # We only check the binding is present somewhere; the existing
    # F2w `q` disable test would catch a regression in the dynamic
    # behaviour.
    assert ':disabled="isParamDisabled(p)"' in body

    # ── Enumerator-owned params hidden in exhaustive mode ────────────
    # The visibleNonStructural() filter must reject params produced
    # by MarsaXorshiftEnumerator::at(idx) when mode === 'exhaustive'.
    # Spot-check the per-type enumerator-owned sets are wired in.
    assert "marsaParamEnumerated_(p)" in body
    # Type 100 deliberately keeps nb_taps / shifts_per_tap editable
    # in exhaustive mode (they're enumerator-time pins, not outputs).
    assert "shifts_per_tap" in body
    assert "nb_taps" in body

    # ── Hide-not-gray for type-irrelevant params ─────────────────────
    # The structural and non-structural cards filter on isParamHidden
    # so type-irrelevant params don't render at all.  The previous
    # gray-out hint ("Pinned to 1 for type 1.") should be gone now.
    assert "isParamHidden" in body
    assert "Pinned to 1 for type 1" not in body
    assert "Not used by type" not in body

    # ── Per-element int_vec editor wired in ──────────────────────────
    # intVecLength_ returns the expected slot count per param, and
    # the template renders one input per slot when the length is
    # known.  syncIntVecLengths_ keeps the backing array sized.
    assert "intVecLength_(p)" in body
    assert "syncIntVecLengths_" in body
    # The per-element render uses x-for over the length and writes
    # into inputs[p.name][i - 1].
    assert "inputs[p.name][i - 1]" in body

    # ── canSubmit skips hidden structural params ─────────────────────
    # Without this, the Launch button stays disabled for MarsaXorshift
    # (m is hidden for type 1/3/100, r is hidden for type 1) — the
    # original UX bug ("Launch search button does not become enabled").
    # The presence of this skip in the rendered template is the
    # canonical guard against regression.
    assert "if (this.isParamHidden(p)) continue" in body
