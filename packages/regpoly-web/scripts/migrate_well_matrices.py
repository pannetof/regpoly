#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Standalone CLI for the WELL matrices payload migration.

The paper-aligned redesign replaced the flat-triple `mat_types`/
`mat_pi`/`mat_pu` keys in JSON column payloads with a single
structured `matrices: {T0..T7}` map. This CLI walks every JSON
column that may carry a WELL row and rewrites the shape in place.

Supersedes the 45cf403 `migrate_well_mat_types.py` (which only did
integer remap). This CLI handles every state via the marker-row
state machine in `well_matrices.py`:

  - Fresh DB: integer remap + shape conversion + drop legacy keys.
  - 45cf403-marker present: shape conversion only (integer remap was
    already done).
  - This-marker present: skip; nothing to do.

Usage:
    uv run python packages/regpoly-web/scripts/migrate_well_matrices.py \\
        --db <db_path> [--dry-run] [--force]

Recommended workflow:
  1. Stop the regpoly-web server (concurrent writes may interleave).
  2. Back up var/regpoly.db.
  3. Run with --dry-run on a copy of the DB; review the per-table
     counters and the skipped-rows report.
  4. Run for real on the live DB.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="migrate_well_matrices",
        description=(
            "Convert WELL `mat_types`/`mat_pi`/`mat_pu` JSON payloads "
            "to the paper-aligned `matrices: {T0..T7}` form."
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

    from regpoly_web.migrations.well_matrices import migrate

    conn = sqlite3.connect(str(args.db))
    try:
        conn.execute("PRAGMA foreign_keys=ON")
        try:
            counters = migrate(conn, dry_run=args.dry_run, force=args.force)
        except RuntimeError as exc:
            print(f"error: {exc}", file=sys.stderr)
            conn.rollback()
            return 3

        print("=== well matrices migration ===")
        print(f"  starting state: {counters.state}")
        print(f"  rows scanned:   {counters.rows_scanned}")
        print(f"  rows touched:   {counters.rows_touched}")
        print(f"  rows skipped:   {counters.rows_skipped}")
        if counters.skipped_old_type6:
            print("  skipped (legacy old-type-6, no paper Mi equivalent):")
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
