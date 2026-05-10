# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — security guards: path-traversal + CSV-injection."""

from __future__ import annotations

import csv
import io


def test_import_dir_rejects_path_outside_root(client, tmp_path) -> None:
    """`/api/import/generators-dir` must confine the directory to a
    configured `Settings.import_root`. A request with /etc must 400,
    not walk the filesystem."""
    r = client.post(
        "/api/import/generators-dir",
        json={"directory": "/etc"},
    )
    assert r.status_code in (400, 403, 404), (
        f"path-traversal not guarded: status={r.status_code}, body={r.text}"
    )
    body = r.text.lower()
    assert "outside" in body or "not allowed" in body or "forbidden" in body \
        or "not a directory" in body, (
        "rejection must be explicit, not a 500"
    )


def test_csv_export_escapes_formula_injection(seeded_client, tmp_path) -> None:
    """A cell starting with `=`, `+`, `-`, `@`, `\\t`, `\\r` must be
    prefixed with `'` so Excel/LibreOffice don't execute the formula.
    Test by inserting a generator with a malicious-looking param value
    and exporting via /api/v2/generators/export."""
    import psycopg
    import json
    db_url = seeded_client.app.state.settings.db_url
    conn = psycopg.connect(db_url)
    try:
        conn.execute(
            "INSERT INTO primitive_generator "
            "(family, l, k, structural_params, search_params, "
            " all_params) VALUES (%s, %s, %s, %s, %s, %s)",
            ("MTGen", 64, 32,
             json.dumps({"w": 32}),
             json.dumps({"a": "=cmd|' /C calc'!A1"}),
             json.dumps({"w": 32, "a": "=cmd|' /C calc'!A1"})),
        )
        conn.commit()
        gen_id = conn.execute(
            "SELECT id FROM primitive_generator ORDER BY id DESC LIMIT 1"
        ).fetchone()[0]
    finally:
        conn.close()

    r = seeded_client.get(
        f"/api/v2/generators/export?ids={gen_id}&fmt=csv"
    )
    assert r.status_code == 200
    # Parse the CSV body; the 'a' cell must be sanitized (start with `'`).
    reader = csv.reader(io.StringIO(r.text))
    rows = list(reader)
    assert len(rows) >= 2, "CSV must have header + data rows"
    # Find the offending value across all data cells.
    for cell in rows[1]:
        if "cmd" in cell:
            assert cell.startswith("'"), (
                "CSV-injection: cell starting with `=` must be prefixed "
                f"with single quote, got {cell!r}"
            )
            return
    # If we never find the cell, the export missed the row entirely; fail.
    raise AssertionError("malicious cell did not appear in CSV export")


def test_csv_export_escapes_plus_minus_at(seeded_client) -> None:
    """The same rule applies to leading `+`, `-`, `@`."""
    import psycopg
    import json
    db_url = seeded_client.app.state.settings.db_url
    conn = psycopg.connect(db_url)
    try:
        for prefix in ("+evil", "-evil", "@evil"):
            conn.execute(
                "INSERT INTO primitive_generator "
                "(family, l, k, structural_params, search_params, all_params) "
                "VALUES (%s, %s, %s, %s, %s, %s)",
                ("MTGen", 64, 32,
                 json.dumps({"w": 32}),
                 json.dumps({"a": prefix}),
                 json.dumps({"w": 32, "a": prefix})),
            )
        conn.commit()
    finally:
        conn.close()

    r = seeded_client.get("/api/v2/generators/export?ids=1,2,3,4&fmt=csv")
    assert r.status_code == 200
    body = r.text
    for prefix in ("+evil", "-evil", "@evil"):
        # The escaped form prefixes a single quote.
        assert (f"'{prefix}" in body) or (prefix not in body), (
            f"Cell starting with {prefix!r} must be sanitized"
        )
