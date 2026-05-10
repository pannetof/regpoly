# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 2 red — /api/v2/dashboard/summary contract.

One single fetch on dashboard load:
  {counts: {generators_tested, primitive_searches_total,
            tempering_searches_total, finds_last_24h},
   active: [{id, type, family, k, status, started_at,
             rate, eta_seconds, sparkline}],
   recent: [{id, type, family, status, finished_at, sparkline}],
   recent_papers: [{id, display}]}
"""

from __future__ import annotations


def test_v2_summary_returns_200_and_top_level_shape(seeded_client) -> None:
    r = seeded_client.get("/api/v2/dashboard/summary")
    assert r.status_code == 200
    body = r.json()
    assert set(body.keys()) >= {"counts", "active", "recent", "recent_papers"}


def test_v2_summary_counts_shape(seeded_client) -> None:
    body = seeded_client.get("/api/v2/dashboard/summary").json()
    counts = body["counts"]
    for key in ("generators_tested", "primitive_searches_total",
                "tempering_searches_total", "finds_last_24h"):
        assert key in counts, f"counts.{key} missing"
        assert isinstance(counts[key], int)


def test_v2_summary_active_rows_carry_rate_eta_sparkline(seeded_client) -> None:
    """Each active row must carry the live-update fields the dashboard
    uses to render the row. Empty active list is acceptable in the
    seeded DB (no running runs); the *contract* of the rows is what we
    pin here."""
    body = seeded_client.get("/api/v2/dashboard/summary").json()
    for row in body["active"]:
        for key in ("id", "type", "family", "status",
                    "rate", "eta_seconds", "sparkline"):
            assert key in row, f"active row missing {key}"
        assert row["type"] in ("primitive", "tempering")
        assert isinstance(row["sparkline"], list)


def test_v2_summary_recent_rows_carry_static_sparkline(seeded_client) -> None:
    """Recent (closed) rows expose static sparklines too — the
    dashboard renders them via sparkline_svg, no live SSE."""
    body = seeded_client.get("/api/v2/dashboard/summary").json()
    for row in body["recent"]:
        assert "id" in row and "type" in row and "status" in row
        assert "sparkline" in row
        assert isinstance(row["sparkline"], list)


def test_v2_summary_recent_papers(seeded_client) -> None:
    body = seeded_client.get("/api/v2/dashboard/summary").json()
    papers = body["recent_papers"]
    assert isinstance(papers, list)
    for p in papers:
        assert "id" in p and "display" in p


def test_v2_summary_includes_recent_searches_last_10(seeded_client) -> None:
    """Recent slot caps at 10 closed/cancelled/failed runs."""
    body = seeded_client.get("/api/v2/dashboard/summary").json()
    assert len(body["recent"]) <= 10


# ── P6 red — dashboard live values are populated, not None ─────────────


def test_summary_active_row_rate_and_eta_are_not_always_none(client) -> None:
    """The dashboard's Rate and ETA columns must reflect *live* values,
    not be hard-coded to None. We seed an active primitive run with a
    progress row carrying a known tries_done, then assert that the
    summary endpoint surfaces a non-None `rate` (computed from the
    run's elapsed_seconds + tries_done, or the per-run RollingRate
    snapshot)."""
    import psycopg
    import json
    db_path = client.app.state.settings.db_url
    conn = psycopg.connect(db_path, autocommit=False)
    try:
        conn.execute(
            "INSERT INTO primitive_search_run "
            "(family, L, k, structural_params, fixed_params, status, "
            "tries_done, found_count, elapsed_seconds, started_at) "
            "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s, NOW())",
            ("MTGen", 64, 32, json.dumps({"w": 32}), json.dumps({}),
             "running", 1000, 5, 10.0),
        )
        run_id = conn.execute(
            "SELECT id FROM primitive_search_run ORDER BY id DESC LIMIT 1"
        ).fetchone()[0]
        conn.execute(
            "INSERT INTO search_progress "
            "(search_type, search_run_id, tries_done, found_count, "
            " current_info) VALUES ('primitive', %s, 1000, 5, %s)",
            (run_id, json.dumps({"rate": 100.0})),
        )
        conn.commit()
    finally:
        conn.close()

    body = client.get("/api/v2/dashboard/summary").json()
    active = [r for r in body["active"] if r["id"] == run_id]
    assert len(active) == 1
    row = active[0]
    # Either rate or eta_seconds (or both) must now be non-None for an
    # actively-running seeded run with progress rows.
    assert row["rate"] is not None or row["eta_seconds"] is not None, (
        "dashboard active row stuck at rate=None, eta_seconds=None — "
        "rate plumbing not wired"
    )


def test_summary_sparkline_uses_tries_done_not_found_count(client) -> None:
    """For primitive runs the sparkline metric must be `tries_done` —
    `found_count` is near-zero and useless as a progress indicator."""
    import psycopg
    import json
    db_path = client.app.state.settings.db_url
    conn = psycopg.connect(db_path, autocommit=False)
    try:
        conn.execute(
            "INSERT INTO primitive_search_run "
            "(family, L, k, structural_params, fixed_params, status, "
            "tries_done, found_count, elapsed_seconds) "
            "VALUES (%s,%s,%s,%s,%s,'running',%s,%s,%s)",
            ("MTGen", 64, 32, json.dumps({}), json.dumps({}), 0, 0, 0.0),
        )
        run_id = conn.execute(
            "SELECT id FROM primitive_search_run ORDER BY id DESC LIMIT 1"
        ).fetchone()[0]
        # Insert progress rows where tries_done grows but found_count is 0.
        for i in range(10):
            conn.execute(
                "INSERT INTO search_progress "
                "(search_type, search_run_id, tries_done, found_count) "
                "VALUES ('primitive', %s, %s, 0)",
                (run_id, (i + 1) * 1000),
            )
        conn.commit()
    finally:
        conn.close()

    body = client.get("/api/v2/dashboard/summary").json()
    active = [r for r in body["active"] if r["id"] == run_id]
    assert len(active) == 1
    sparkline = active[0]["sparkline"]
    # If we plotted found_count we'd get [0,0,…,0]. tries_done gives a
    # strictly-increasing sequence.
    assert any(p > 0 for p in sparkline), (
        "sparkline appears to plot found_count (all zero); should be "
        f"tries_done; got {sparkline}"
    )
