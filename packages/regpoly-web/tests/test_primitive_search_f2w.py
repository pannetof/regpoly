# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Web-layer tests for the F2w paper-notation primitive search.

A user enters {w, r, nb_terms} (and optionally t, q) in the form; the
server randomizes {t, q, modM, coeff} per iteration via the F2w
param_specs (rand_type ``range`` for t/q, ``irreducible_gf2`` for modM,
``bitmask_vec`` for coeff with scalar length ``nb_terms``).
"""

from __future__ import annotations


def _f2w_body(**overrides) -> dict:
    body = {
        "family": "F2wLFSRGen",
        "L": 32,
        "structural_params": {"w": 32, "r": 8, "nb_terms": 3},
        "fixed_params": {"normal_basis": False, "step": 1},
        "max_tries": 200,
    }
    body.update(overrides)
    return body


def test_f2w_paper_notation_search_persists(client) -> None:
    """POSTing an F2w paper-notation search produces a run row whose
    persisted params match the input."""
    r = client.post("/api/primitive-searches", json=_f2w_body())
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]

    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    assert detail["family"] == "F2wLFSRGen"
    assert detail["L"] == 32
    assert detail["max_tries"] == 200


def test_f2w_search_accepts_fixed_t(client) -> None:
    """Fixing `t` in fixed_params is the paper-form 'fix t' affordance —
    the search loop must skip sampling t when it's already in p."""
    body = _f2w_body(fixed_params={"normal_basis": False, "step": 1, "t": 5})
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text


def test_f2w_search_accepts_nb_terms_equals_2(client) -> None:
    """nb_terms=2 polynomials drop the `q` term; from_params builds nocoeff=[t,0]."""
    body = _f2w_body(structural_params={"w": 32, "r": 8, "nb_terms": 2})
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text


def test_f2w_polylcg_variant_accepted(client) -> None:
    """The same paper-form payload works for F2wPolyLCGGen (twins
    share param_specs)."""
    body = _f2w_body(family="F2wPolyLCGGen")
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
