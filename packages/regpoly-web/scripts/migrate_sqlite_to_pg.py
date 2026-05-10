#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""One-shot SQLite → PostgreSQL data migration (dockerize-plan Phase 1).

Read every row from a legacy ``regpoly.db`` SQLite file and INSERT it
into a fresh PG database. The PG schema must already exist (apply
``regpoly_web.migrations.m001_initial`` and ``m002_library_test_run``
first via ``init_db``, e.g. by running ``regpoly-web`` once or by the
init container in compose).

Type-coercion rules:

- ``datetime('now')`` strings → ``TIMESTAMPTZ``: parsed with
  ``datetime.fromisoformat``; naive timestamps are assumed UTC.
- TEXT-as-JSON columns → ``JSONB``: parsed with ``json.loads``.
  Rows whose JSON fails to parse are skipped and logged (their IDs
  appear in the per-table counter as ``skipped``).
- ``INTEGER`` boolean columns (``verified``, ``is_me``, ``is_cf``)
  → ``BOOLEAN``: ``0/1`` mapped to ``False/True``; NULL preserved.
- ``BIGSERIAL`` sequences are advanced via ``setval`` after the bulk
  insert so subsequent INSERTs from the application don't collide
  with imported PKs.

Usage::

    uv run python packages/regpoly-web/scripts/migrate_sqlite_to_pg.py \\
        var/regpoly.db [--target-url postgresql://...] [--dry-run]

Exit codes:
  0  success
  2  source SQLite file missing or unreadable
  3  runtime error (PG connection, parse failure beyond skip threshold,
     etc.)

Recommended workflow:
  1. Stop the regpoly-web stack (no concurrent writers).
  2. Verify the PG database is empty: row counts in every target table
     must be zero. The script aborts otherwise unless ``--force`` is
     passed.
  3. Run with ``--dry-run`` first; review the per-table counters.
  4. Run for real.
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import psycopg

# Tables to migrate, in dependency order (parents before children).
# Each entry: (table_name, [(sqlite_col, pg_col, coercion), ...]).
# coercion is one of: "scalar", "json", "datetime", "bool".
_TABLES: list[tuple[str, list[tuple[str, str, str]]]] = [
    ("schema_version", [
        ("version", "version", "scalar"),
        ("applied_at", "applied_at", "datetime"),
    ]),
    ("primitive_search_run", [
        ("id", "id", "scalar"),
        ("family", "family", "scalar"),
        ("L", '"L"', "scalar"),
        ("k", "k", "scalar"),
        ("structural_params", "structural_params", "json"),
        ("fixed_params", "fixed_params", "json"),
        ("max_tries", "max_tries", "scalar"),
        ("max_seconds", "max_seconds", "scalar"),
        ("max_cost", "max_cost", "scalar"),
        ("status", "status", "scalar"),
        ("tries_done", "tries_done", "scalar"),
        ("found_count", "found_count", "scalar"),
        ("elapsed_seconds", "elapsed_seconds", "scalar"),
        ("error_message", "error_message", "scalar"),
        ("created_at", "created_at", "datetime"),
        ("started_at", "started_at", "datetime"),
        ("finished_at", "finished_at", "datetime"),
        ("search_mode", "search_mode", "scalar"),
        ("enum_index", "enum_index", "scalar"),
        ("enum_total", "enum_total", "scalar"),
        ("enum_axes", "enum_axes", "json"),
    ]),
    ("primitive_generator", [
        ("id", "id", "scalar"),
        ("search_run_id", "search_run_id", "scalar"),
        ("family", "family", "scalar"),
        ("L", '"L"', "scalar"),
        ("k", "k", "scalar"),
        ("structural_params", "structural_params", "json"),
        ("search_params", "search_params", "json"),
        ("all_params", "all_params", "json"),
        ("found_at_try", "found_at_try", "scalar"),
        ("created_at", "created_at", "datetime"),
        ("char_poly", "char_poly", "scalar"),
        ("hamming_weight", "hamming_weight", "scalar"),
        ("pis_gaps", "pis_gaps", "json"),
        ("pis_se", "pis_se", "scalar"),
        ("pis_elapsed", "pis_elapsed", "scalar"),
        ("pis_computed_at", "pis_computed_at", "datetime"),
        ("pis_error", "pis_error", "scalar"),
        ("library_id", "library_id", "scalar"),
    ]),
    ("tempering_search_run", [
        ("id", "id", "scalar"),
        ("test_type", "test_type", "scalar"),
        ("test_config", "test_config", "json"),
        ("Lmax", '"Lmax"', "scalar"),
        ("nb_tries", "nb_tries", "scalar"),
        ("optimizer_config", "optimizer_config", "json"),
        ("status", "status", "scalar"),
        ("combos_total", "combos_total", "scalar"),
        ("combos_done", "combos_done", "scalar"),
        ("best_se", "best_se", "scalar"),
        ("elapsed_seconds", "elapsed_seconds", "scalar"),
        ("error_message", "error_message", "scalar"),
        ("created_at", "created_at", "datetime"),
        ("started_at", "started_at", "datetime"),
        ("finished_at", "finished_at", "datetime"),
    ]),
    ("tempering_search_component", [
        ("id", "id", "scalar"),
        ("search_run_id", "search_run_id", "scalar"),
        ("component_index", "component_index", "scalar"),
        ("shared_with_component", "shared_with_component", "scalar"),
        ("tempering_config", "tempering_config", "json"),
    ]),
    ("tempering_search_generator", [
        ("id", "id", "scalar"),
        ("component_id", "component_id", "scalar"),
        ("generator_id", "generator_id", "scalar"),
    ]),
    ("tested_generator", [
        ("id", "id", "scalar"),
        ("search_run_id", "search_run_id", "scalar"),
        ("Lmax", '"Lmax"', "scalar"),
        ("k_g", "k_g", "scalar"),
        ("J", '"J"', "scalar"),
        ("created_at", "created_at", "datetime"),
        ("library_id", "library_id", "scalar"),
    ]),
    ("tested_generator_component", [
        ("id", "id", "scalar"),
        ("tested_gen_id", "tested_gen_id", "scalar"),
        ("component_index", "component_index", "scalar"),
        ("generator_id", "generator_id", "scalar"),
        ("family", "family", "scalar"),
        ("L", '"L"', "scalar"),
        ("k", "k", "scalar"),
        ("all_params", "all_params", "json"),
        ("tempering_params", "tempering_params", "json"),
    ]),
    ("test_result", [
        ("id", "id", "scalar"),
        ("tested_gen_id", "tested_gen_id", "scalar"),
        ("test_type", "test_type", "scalar"),
        ("test_config", "test_config", "json"),
        ("se", "se", "scalar"),
        ("is_me", "is_me", "bool"),
        ("secf", "secf", "scalar"),
        ("is_cf", "is_cf", "bool"),
        ("score", "score", "scalar"),
        ("detail", "detail", "json"),
        ("elapsed_seconds", "elapsed_seconds", "scalar"),
        ("created_at", "created_at", "datetime"),
    ]),
    ("equidistribution_result", [
        ("id", "id", "scalar"),
        ("tested_gen_id", "tested_gen_id", "scalar"),
        ("test_config", "test_config", "json"),
        ("kg", "kg", "scalar"),
        ("L", '"L"', "scalar"),
        ("Lmax", '"Lmax"', "scalar"),
        ("ecart_json", "ecart_json", "json"),
        ("se", "se", "scalar"),
        ("verified", "verified", "bool"),
        ("elapsed_seconds", "elapsed_seconds", "scalar"),
        ("created_at", "created_at", "datetime"),
    ]),
    ("collision_free_result", [
        ("id", "id", "scalar"),
        ("tested_gen_id", "tested_gen_id", "scalar"),
        ("test_config", "test_config", "json"),
        ("kg", "kg", "scalar"),
        ("L", '"L"', "scalar"),
        ("ecart_cf_json", "ecart_cf_json", "json"),
        ("secf", "secf", "scalar"),
        ("verified", "verified", "bool"),
        ("elapsed_seconds", "elapsed_seconds", "scalar"),
        ("created_at", "created_at", "datetime"),
    ]),
    ("tuplets_result", [
        ("id", "id", "scalar"),
        ("tested_gen_id", "tested_gen_id", "scalar"),
        ("test_config", "test_config", "json"),
        ("kg", "kg", "scalar"),
        ("L", '"L"', "scalar"),
        ("tupd", "tupd", "scalar"),
        ("testtype", "testtype", "scalar"),
        ("threshold", "threshold", "scalar"),
        ("tuph_json", "tuph_json", "json"),
        ("gap_json", "gap_json", "json"),
        ("DELTA_json", '"DELTA_json"', "json"),
        ("pourcentage_json", "pourcentage_json", "json"),
        ("firstpart_max", "firstpart_max", "scalar"),
        ("firstpart_sum", "firstpart_sum", "scalar"),
        ("secondpart_max", "secondpart_max", "scalar"),
        ("secondpart_sum", "secondpart_sum", "scalar"),
        ("elapsed_seconds", "elapsed_seconds", "scalar"),
        ("created_at", "created_at", "datetime"),
    ]),
    ("search_progress", [
        ("id", "id", "scalar"),
        ("search_type", "search_type", "scalar"),
        ("search_run_id", "search_run_id", "scalar"),
        ("tries_done", "tries_done", "scalar"),
        ("found_count", "found_count", "scalar"),
        ("current_info", "current_info", "scalar"),
        ("message", "message", "scalar"),
        ("updated_at", "updated_at", "datetime"),
    ]),
    ("yaml_import", [
        ("id", "id", "scalar"),
        ("file_path", "file_path", "scalar"),
        ("import_type", "import_type", "scalar"),
        ("row_count", "row_count", "scalar"),
        ("imported_at", "imported_at", "datetime"),
    ]),
]


def _coerce(value, kind: str):
    """Apply per-column coercion rules. Raises ValueError on JSON failure."""
    if value is None:
        return None
    if kind == "scalar":
        return value
    if kind == "datetime":
        # SQLite stores ISO-8601 via datetime('now'). Naive → assume UTC.
        if isinstance(value, datetime):
            return value
        s = str(value).replace(" ", "T")  # tolerate "YYYY-MM-DD HH:MM:SS"
        dt = datetime.fromisoformat(s)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt
    if kind == "json":
        if isinstance(value, (dict, list)):
            return psycopg.types.json.Jsonb(value)
        s = str(value).strip()
        if not s:
            return None
        return psycopg.types.json.Jsonb(json.loads(s))
    if kind == "bool":
        if isinstance(value, bool):
            return value
        return bool(int(value))
    raise ValueError(f"unknown coercion kind: {kind}")


def _ensure_pg_empty(pg: psycopg.Connection, force: bool) -> None:
    """Abort if any target table has rows (unless --force)."""
    if force:
        return
    nonempty: list[tuple[str, int]] = []
    for table, _ in _TABLES:
        if table == "schema_version":
            continue  # init_db populates this; ignore
        with pg.cursor() as cur:
            cur.execute(f"SELECT COUNT(*) FROM {table}")
            row = cur.fetchone()
            n = row[0] if row else 0
        if n > 0:
            nonempty.append((table, n))
    if nonempty:
        details = "; ".join(f"{t}={n}" for t, n in nonempty)
        raise RuntimeError(
            f"target PG database is not empty ({details}). "
            "Pass --force to override (will INSERT alongside existing rows)."
        )


def _migrate_table(
    sqlite: sqlite3.Connection,
    pg: psycopg.Connection,
    table: str,
    cols: list[tuple[str, str, str]],
    dry_run: bool,
) -> tuple[int, int]:
    """Migrate one table; return (inserted, skipped)."""
    src_cols = ", ".join(f'"{c[0]}"' for c in cols)
    dst_cols = ", ".join(c[1] for c in cols)
    placeholders = ", ".join(["%s"] * len(cols))
    insert_sql = (
        f"INSERT INTO {table} ({dst_cols}) VALUES ({placeholders}) "
        "ON CONFLICT DO NOTHING"
    )
    cur = sqlite.execute(f'SELECT {src_cols} FROM "{table}"')
    inserted, skipped = 0, 0
    batch: list[tuple] = []
    for row in cur:
        try:
            coerced = tuple(_coerce(v, kind) for v, (_, _, kind) in zip(row, cols))
        except (ValueError, json.JSONDecodeError) as e:
            print(f"  [skip] {table} row id={row[0]}: {e}", file=sys.stderr)
            skipped += 1
            continue
        batch.append(coerced)
        if len(batch) >= 500:
            if not dry_run:
                with pg.cursor() as pcur:
                    pcur.executemany(insert_sql, batch)
            inserted += len(batch)
            batch = []
    if batch:
        if not dry_run:
            with pg.cursor() as pcur:
                pcur.executemany(insert_sql, batch)
        inserted += len(batch)
    if not dry_run:
        pg.commit()
    return inserted, skipped


def _setval_sequences(pg: psycopg.Connection) -> None:
    """Advance every BIGSERIAL sequence to MAX(id) + 1."""
    for table, cols in _TABLES:
        if not any(c[0] == "id" for c in cols):
            continue
        seq = f"{table}_id_seq"
        with pg.cursor() as cur:
            cur.execute(
                f"SELECT setval('{seq}', GREATEST((SELECT COALESCE(MAX(id), 0) FROM {table}), 1))"
            )
    pg.commit()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="migrate_sqlite_to_pg",
        description="One-shot SQLite → PostgreSQL data migration.",
    )
    parser.add_argument(
        "db_path",
        type=Path,
        help="Path to the source SQLite DB file.",
    )
    parser.add_argument(
        "--target-url",
        default=os.environ.get("REGPOLY_DB_URL"),
        help="Destination PG DSN (default: $REGPOLY_DB_URL).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read source + parse, but do not write to PG.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow INSERT into a non-empty PG database.",
    )
    args = parser.parse_args(argv)

    if not args.db_path.is_file():
        print(f"error: source SQLite file not found: {args.db_path}", file=sys.stderr)
        return 2
    if not args.target_url:
        print(
            "error: --target-url not provided and REGPOLY_DB_URL is unset",
            file=sys.stderr,
        )
        return 3

    print(f"=== sqlite → pg migration ===")
    print(f"source: {args.db_path}")
    print(f"target: {args.target_url}")
    print(f"mode:   {'dry-run' if args.dry_run else 'live'}")
    print()

    try:
        sqlite = sqlite3.connect(str(args.db_path))
        sqlite.row_factory = sqlite3.Row
    except sqlite3.Error as e:
        print(f"error: cannot open SQLite DB: {e}", file=sys.stderr)
        return 2

    try:
        pg = psycopg.connect(args.target_url)
    except psycopg.OperationalError as e:
        print(f"error: cannot connect to PG: {e}", file=sys.stderr)
        return 3

    try:
        _ensure_pg_empty(pg, args.force)
    except RuntimeError as e:
        print(f"error: {e}", file=sys.stderr)
        return 3

    try:
        total_in, total_skip = 0, 0
        for table, cols in _TABLES:
            if table == "schema_version":
                continue  # already populated by init_db
            try:
                ins, skp = _migrate_table(sqlite, pg, table, cols, args.dry_run)
            except sqlite3.OperationalError as e:
                # Source table missing (e.g. older SQLite DB without v2 tables).
                print(f"  [skip table] {table}: {e}", file=sys.stderr)
                continue
            total_in += ins
            total_skip += skp
            print(f"  {table:<32} +{ins:>8}  (skipped {skp})")

        if not args.dry_run:
            _setval_sequences(pg)
            print()
            print("sequences advanced.")

        print()
        print(f"total: inserted={total_in} skipped={total_skip}")
        return 0
    except Exception as e:
        print(f"runtime error: {e}", file=sys.stderr)
        return 3
    finally:
        sqlite.close()
        pg.close()


if __name__ == "__main__":
    raise SystemExit(main())
