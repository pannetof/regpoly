# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 5 red — /generators/compare?ids=… side-by-side view."""

from __future__ import annotations


def test_compare_route_with_two_ids(seeded_client) -> None:
    r = seeded_client.get("/generators/compare?ids=1,1")
    assert r.status_code == 200, r.text
    body = r.text
    assert "compare" in body.lower()


def test_v2_compare_endpoint_returns_one_row_per_id(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/compare?ids=1,1")
    assert r.status_code == 200, r.text
    body = r.json()
    assert "rows" in body
    assert len(body["rows"]) == 2


def test_v2_compare_endpoint_skips_missing_ids(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/compare?ids=1,99999")
    assert r.status_code == 200
    body = r.json()
    # Missing ids appear as null entries so the UI can render a hole.
    assert len(body["rows"]) == 2
    assert any(row is None for row in body["rows"])


# ── P6 red — compare endpoint extensions ─────────────────────────────


def test_v2_compare_includes_equidist_data_for_overlay(seeded_client) -> None:
    """The compare endpoint must return equidistribution data so the
    UI can overlay δ profiles on a shared ℓ-axis. Today the response
    has only primitive_generator columns + all_params; chart overlay
    needs `pis_gaps`."""
    r = seeded_client.get("/api/v2/generators/compare?ids=1")
    assert r.status_code == 200
    body = r.json()
    if not body["rows"] or body["rows"][0] is None:
        return  # nothing to assert if the seed row was wiped
    row = body["rows"][0]
    assert "pis_gaps" in row or "pis_se" in row, (
        "compare row must include pis_gaps/pis_se for chart overlay"
    )


def test_v2_tested_generators_compare_endpoint(seeded_client) -> None:
    r = seeded_client.get("/api/v2/tested-generators/compare?ids=4242,4242")
    assert r.status_code == 200, r.text
    body = r.json()
    assert "rows" in body
    assert len(body["rows"]) == 2


def test_v2_library_compare_endpoint(seeded_client) -> None:
    """Library compare endpoint exists. Library entries don't
    necessarily have integer ids; accept whatever id strings the
    catalog uses. Empty list still returns 200 with rows=[]."""
    r = seeded_client.get("/api/v2/library/compare?ids=")
    assert r.status_code == 200
    body = r.json()
    assert "rows" in body
