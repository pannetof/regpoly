# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 2 red — chip toolbar filter contract on v2 list endpoints.

Every chip → a query-string facet. AND across facets. URL is the
single source of truth. Same grammar reused by /generators,
/tested-generators, /library.
"""

from __future__ import annotations


def test_v2_generators_filter_by_family(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators?family=MTGen&limit=10")
    assert r.status_code == 200
    body = r.json()
    for row in body["rows"]:
        assert row["family"] == "MTGen"


def test_v2_generators_filter_by_k_range(seeded_client) -> None:
    body = seeded_client.get(
        "/api/v2/generators?k_min=10&k_max=10000&limit=10"
    ).json()
    for row in body["rows"]:
        assert 10 <= row["k"] <= 10000


def test_v2_generators_has_search_results_filter(seeded_client) -> None:
    """`?has_search_results=true` returns only generators whose
    primitive_search_run found at least one result; `=false` is
    the complement."""
    yes = seeded_client.get(
        "/api/v2/generators?has_search_results=true&limit=100"
    ).json()
    no = seeded_client.get(
        "/api/v2/generators?has_search_results=false&limit=100"
    ).json()
    assert isinstance(yes["rows"], list)
    assert isinstance(no["rows"], list)
    # No row should appear in both sides of the partition.
    yes_ids = {r["id"] for r in yes["rows"]}
    no_ids = {r["id"] for r in no["rows"]}
    assert yes_ids.isdisjoint(no_ids)


def test_v2_tested_generators_has_results_filter(seeded_client) -> None:
    """`?has_results=true` joins against the v2 typed result tables
    (equidistribution_result, collision_free_result, tuplets_result).
    The seeded tg id 4242 has one equid row so it must appear here."""
    yes = seeded_client.get(
        "/api/v2/tested-generators?has_results=true&limit=100"
    ).json()
    ids = {r["id"] for r in yes["rows"]}
    assert 4242 in ids, (
        "tested_generator id=4242 has equidistribution_result; "
        "?has_results=true should include it"
    )


def test_v2_tested_generators_chip_facets_kg_test_type_max_sum_delta_is_me(
    seeded_client,
) -> None:
    """Per Analyst persona A-4: enumerate the chip facets explicitly.

    These params must be accepted (200) even when the seeded DB
    contains only one row. Filter effects are exercised in green
    commit; red commit just pins the contract surface."""
    r = seeded_client.get(
        "/api/v2/tested-generators"
        "?family=MTGen&k_g_min=10&k_g_max=100000"
        "&test_type=equidistribution&max_sum_delta=100&is_me=false"
        "&limit=10"
    )
    assert r.status_code == 200, (
        f"all chip facets must be accepted; got {r.status_code}"
    )


def test_v2_generators_families_counts(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/families/counts")
    assert r.status_code == 200
    body = r.json()
    assert isinstance(body, dict)
    assert body.get("MTGen", 0) >= 1, "seeded DB has one MTGen generator"


def test_v2_library_families_counts(seeded_client) -> None:
    """Per Cataloger persona C-10: parity with generators chip filter."""
    r = seeded_client.get("/api/v2/library/families/counts")
    assert r.status_code == 200
    assert isinstance(r.json(), dict)


def test_v2_searches_histories_batch(seeded_client) -> None:
    """Per Researcher persona R-9: avoid N+1 fetches on /searches."""
    r = seeded_client.get(
        "/api/v2/searches/histories?ids=1&type=primitive&points=30"
    )
    assert r.status_code == 200
    body = r.json()
    assert "histories" in body
    assert isinstance(body["histories"], dict)


# ── P6 red — p_<name> filter must work for column-level keys too ─────


def test_p_hamming_weight_filters_against_column_not_json(seeded_client) -> None:
    """Today the v2 endpoint only honours `p_<name>` against
    json_extract(all_params, …). The seeded generator has
    hamming_weight=135 stored as a primitive_generator column, NOT in
    all_params. A user clicking the column popover on Hamming weight
    silently zero-results.

    Fix: widen the popover-filter to also try a top-level column
    match against an allowlist of safe column names (hamming_weight,
    char_poly, library_id, found_at_try) — OR have the backend honour
    these as top-level filters with their own param name.
    """
    r = seeded_client.get(
        "/api/v2/generators?p_hamming_weight=135&limit=10"
    )
    assert r.status_code == 200
    body = r.json()
    assert body["total"] >= 1, (
        "p_hamming_weight=135 silently zero-results because the column "
        "isn't inside all_params JSON"
    )
