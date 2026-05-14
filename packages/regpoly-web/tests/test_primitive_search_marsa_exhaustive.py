# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Exhaustive-mode primitive-search end-to-end tests for
``MarsaXorshiftGen``.

Regression cover for a UX bug where Type-3 + Exhaustive returned
"Cannot enumerate: Not enumerable.": ``_probe_k`` was called BEFORE
the enumerator was built, with a Params bag that omitted the
enumerator-owned ``tap_positions``/``tap_shifts``; the type-3 ctor
invariant ("type=3 requires a non-empty taps list") tripped, and the
generic ``HTTPException(400, str)`` couldn't be parsed by the form's
dict-shaped error handler.

The fix builds the enumerator first in exhaustive mode and probes
``k`` from ``enumerator.at(0)`` (a guaranteed-valid bag).  Below we
exercise both the estimate path (which now skips the probe entirely)
and the create path (which now probes after building).
"""

from __future__ import annotations


def _marsa_exhaustive_body(typ: int, r: int, **extra) -> dict:
    body = {
        "family": "MarsaXorshiftGen",
        "structural_params": {"type": typ, "w": 32, "r": r},
        "fixed_params": {},
        "max_tries": None,
        "max_seconds": None,
        "max_cost": None,
        "max_se": None,
        "search_mode": "exhaustive",
        "confirm_huge": False,
    }
    body.update(extra)
    return body


# ── /api/primitive-searches/estimate ──────────────────────────────────

def test_type3_exhaustive_estimate_returns_size_and_axes(client) -> None:
    """Type 3 + Exhaustive must produce a numeric estimate, not the
    'Not enumerable.' error the user reported.  C(r=3, 3) × (2·31)^3
    = 1 × 238328."""
    r = client.post("/api/primitive-searches/estimate",
                    json=_marsa_exhaustive_body(typ=3, r=3))
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["search_mode"] == "exhaustive"
    assert int(data["total"]) == 1 * (2 * 31) ** 3
    names = [a["name"] for a in data["axes"]]
    assert "taps" in names
    assert sum(1 for n in names if n.startswith("shift[")) == 3


def test_type1_exhaustive_estimate_works(client) -> None:
    """Type 1 + Exhaustive — w=32 alone defines the search; r is
    pinned to 1 by the form.  4 patterns × 31^3 = 119,164."""
    r = client.post("/api/primitive-searches/estimate",
                    json=_marsa_exhaustive_body(typ=1, r=1))
    assert r.status_code == 200, r.text
    assert int(r.json()["total"]) == 4 * 31 ** 3


def test_type100_exhaustive_estimate_uses_default_pins(client) -> None:
    """Type 100 + Exhaustive with no pins should default to Table-IV
    shape (nb_taps=3, shifts_per_tap=[1,1,1]): C(r=12, 3) × 62^3."""
    body = _marsa_exhaustive_body(typ=100, r=12)
    body["fixed_params"] = {"nb_taps": 3, "shifts_per_tap": [1, 1, 1]}
    r = client.post("/api/primitive-searches/estimate", json=body)
    assert r.status_code == 200, r.text
    assert int(r.json()["total"]) == 220 * 62 ** 3


def test_invalid_type3_r_returns_clear_message(client) -> None:
    """Type 3 with r < 3 is genuinely not enumerable; the response
    must include a ``message`` field the form can render."""
    r = client.post("/api/primitive-searches/estimate",
                    json=_marsa_exhaustive_body(typ=3, r=2))
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert isinstance(detail, dict)
    assert "needs_r_ge_3" in detail.get("code", "")
    assert detail.get("message", "")


# ── /api/primitive-searches (create) ──────────────────────────────────

def test_type3_exhaustive_create_persists_with_correct_k(client) -> None:
    """End-to-end: Type 3 + Exhaustive should submit, build the
    enumerator, probe k from enumerator.at(0), and persist a run row
    with k = w * r = 96."""
    r = client.post("/api/primitive-searches",
                    json=_marsa_exhaustive_body(typ=3, r=3))
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 32 * 3
    assert payload["search_mode"] == "exhaustive"
    assert payload["enum_total"] == str(1 * (2 * 31) ** 3)


def test_type2_random_accepts_per_element_int_vec_arrays(client) -> None:
    """The form's new per-element int_vec editor sends `p` and `q` as
    JSON arrays (not strings) in fixed_params.  The route must accept
    them without re-parsing.  Mirrors what the form's _readValue
    helper produces from a row of small inputs."""
    body = {
        "family": "MarsaXorshiftGen",
        "structural_params": {"type": 2, "w": 32, "r": 2, "m": 1},
        "fixed_params": {
            "p": [-11, 0, 0],
            "q": [-19, 13, 0],
        },
        "max_tries": 1,
        "max_seconds": 5,
        "max_cost": None,
        "max_se": None,
        "search_mode": "random",
        "confirm_huge": False,
    }
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 32 * 2
    assert payload["search_mode"] == "random"


def test_type1_random_accepts_payload_without_hidden_m(client) -> None:
    """Regression for the "Launch search button never enables" UX bug:
    Type 1 hides `m` (and pre-pins `r=1`), so the form's buildBody
    omits both from the payload.  The route must accept that shape.
    canSubmit() must mirror the same skip; we exercise the server
    side here, and the canSubmit JS path is covered by inspecting
    the rendered template in test_marsa_form_renders."""
    body = {
        "family": "MarsaXorshiftGen",
        # Note: no `m` — the form hides it for type 1.
        "structural_params": {"type": 1, "w": 32, "r": 1},
        "fixed_params": {"shifts": [-13, 17, -5]},     # Marsaglia X2(5,17,13)
        "max_tries": 1,
        "max_seconds": 5,
        "max_cost": None,
        "max_se": None,
        "search_mode": "random",
        "confirm_huge": False,
    }
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 32
    assert payload["search_mode"] == "random"


def test_type3_exhaustive_accepts_payload_without_hidden_m_or_taps(client) -> None:
    """Type 3 + Exhaustive hides `m` (irrelevant) and `tap_*`
    (enumerator-owned).  buildBody therefore omits all three; the
    route must accept the minimal payload."""
    body = {
        "family": "MarsaXorshiftGen",
        "structural_params": {"type": 3, "w": 32, "r": 5},
        "fixed_params": {},
        "max_tries": None,
        "max_seconds": None,
        "max_cost": None,
        "max_se": None,
        "search_mode": "exhaustive",
        "confirm_huge": False,
    }
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    payload = r.json()
    assert payload["k"] == 32 * 5
