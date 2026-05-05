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
