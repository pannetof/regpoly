# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 5 red — selection-scoped CSV/JSON export endpoints."""

from __future__ import annotations


def test_v2_generators_csv_returns_csv(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/export?ids=1&fmt=csv")
    assert r.status_code == 200
    assert "text/csv" in r.headers.get("content-type", "")
    body = r.text
    assert "id" in body.lower().split("\n")[0], (
        "CSV must include an id column header"
    )


def test_v2_generators_json_returns_json(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/export?ids=1&fmt=json")
    assert r.status_code == 200
    body = r.json()
    assert isinstance(body, list)
    assert len(body) == 1
    assert body[0].get("id") == 1


def test_v2_generators_export_skips_unknown_ids(seeded_client) -> None:
    r = seeded_client.get(
        "/api/v2/generators/export?ids=1,99999&fmt=json"
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body) == 1
