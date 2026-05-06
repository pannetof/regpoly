# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — GET /api/v2/tested-generators/{id} full-detail.

Today an SDK consumer must mix v1 detail with v2 list + v2 publish (3
round-trips, 2 namespaces). The new endpoint returns the full record
including components, per-test results, library_id, with a Pydantic
response_model for OpenAPI typing.
"""

from __future__ import annotations


def test_v2_tested_generator_detail_returns_full_record(seeded_client) -> None:
    r = seeded_client.get("/api/v2/tested-generators/4242")
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["id"] == 4242
    assert "components" in body and isinstance(body["components"], list)
    assert "library_id" in body  # may be null
    # The seeded fixture has J=1 with one component.
    assert len(body["components"]) >= 1


def test_v2_tested_generator_detail_404_on_unknown(client) -> None:
    r = client.get("/api/v2/tested-generators/999999")
    assert r.status_code == 404
