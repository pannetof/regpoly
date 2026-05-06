#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Standalone CLI for the WELL mat_types renumbering.

The C++ rename to paper M0..M6 numbering changed the semantic meaning
of every WELL ``mat_types`` integer. JSON-encoded payloads in the
SQLite DB (``primitive_generator``, ``primitive_search_run``,
``tested_generator_component``, ``tempering_search_component``,
``test_result``) are NOT covered by the C++ legacy_reader's read-time
remap — this script does the equivalent for those payloads.

Usage:
    uv run python packages/regpoly-web/scripts/migrate_well_mat_types.py \\
        --db <db_path> [--dry-run] [--force]

The migration is idempotent: re-runs are gated by a marker row in
``yaml_import`` (``file_path = 'well-mat-types-migrated'``). Use
``--force`` to re-run anyway. ``--dry-run`` prints the per-table
counters without committing — recommended on a copy of the DB before
applying to live data.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="migrate_well_mat_types",
        description=(
            "Renumber WELL mat_types integers in JSON columns from the "
            "legacy 0..7 scheme to paper M0..M6 (0..6)."
        ),
    )
    parser.add_argument(
        "--db", type=Path, required=True, help="Path to the SQLite DB file."
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Run the migration but roll back the transaction.",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run even if the idempotency marker is already set.",
    )
    args = parser.parse_args(argv)

    if not args.db.exists():
        print(f"error: DB file does not exist: {args.db}", file=sys.stderr)
        return 2

    from regpoly_web.migrations.well_mat_types import migrate

    conn = sqlite3.connect(str(args.db))
    try:
        conn.execute("PRAGMA foreign_keys=ON")
        try:
            counters = migrate(conn, dry_run=args.dry_run, force=args.force)
        except RuntimeError as exc:
            print(f"error: {exc}", file=sys.stderr)
            conn.rollback()
            return 3

        print("=== well mat_types migration ===")
        print(f"  rows scanned:  {counters.rows_scanned}")
        print(f"  rows touched:  {counters.rows_touched}")
        print(f"  rows skipped:  {counters.rows_skipped}")
        if counters.skipped_old_type6:
            print("  skipped (old type 6 — no paper Mi equivalent):")
            for table, rid in counters.skipped_old_type6:
                print(f"    {table}.id={rid}")
        if counters.skipped_other:
            print("  skipped (other reasons):")
            for table, rid, why in counters.skipped_other:
                print(f"    {table}.id={rid}: {why}")
        if counters.per_table:
            print("  per-table:")
            for table, kv in sorted(counters.per_table.items()):
                bits = ", ".join(f"{k}={v}" for k, v in sorted(kv.items()))
                print(f"    {table}: {bits}")

        if args.dry_run:
            conn.rollback()
            print("--dry-run: rolled back.")
        else:
            conn.commit()
            print("committed.")
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
