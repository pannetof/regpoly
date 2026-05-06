# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — per-column sort on v2 list endpoints.

Today /api/v2/generators hard-codes ORDER BY id DESC. Add ?sort=<col>
&dir=<asc|desc> with an allowlist of sortable columns, mirrored on
/api/v2/tested-generators and /api/v2/library.
"""

from __future__ import annotations


def test_v2_generators_sort_by_k_asc(seeded_client) -> None:
    r = seeded_client.get(
        "/api/v2/generators?sort=k&dir=asc&limit=200"
    )
    assert r.status_code == 200, r.text
    rows = r.json()["rows"]
    if len(rows) < 2:
        return
    ks = [row["k"] for row in rows]
    assert ks == sorted(ks), f"sort=k&dir=asc not honoured: {ks}"


def test_v2_generators_sort_rejects_unknown_column(seeded_client) -> None:
    r = seeded_client.get(
        "/api/v2/generators?sort=name'); DROP TABLE--"
    )
    # Either a 400 (rejected) or a 200 (silently ignored) is acceptable,
    # but a 500 (interpolated into SQL) is not.
    assert r.status_code in (200, 400, 422), r.text


def test_v2_tested_generators_sort_by_kg_desc(seeded_client) -> None:
    r = seeded_client.get(
        "/api/v2/tested-generators?sort=k_g&dir=desc&limit=200"
    )
    assert r.status_code == 200
    rows = r.json()["rows"]
    if len(rows) < 2:
        return
    kgs = [row.get("k_g") for row in rows if row.get("k_g") is not None]
    assert kgs == sorted(kgs, reverse=True)
