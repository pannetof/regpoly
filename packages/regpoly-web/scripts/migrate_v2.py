#!/usr/bin/env python3
"""Standalone CLI for the v1 → v2 schema migration.

Used to upgrade an existing var/regpoly.db without starting the web
app. The web app's startup path runs the same migration automatically
when it opens an older DB; this script exists for offline rehearsals
(e.g. before deploying a new server build).

Usage:
    uv run python packages/regpoly-web/scripts/migrate_v2.py <db_path>

The migration is idempotent: running it on an already-v2 DB is a no-op.
A copy of the DB file is recommended before running on production data.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="migrate_v2",
        description="Upgrade a regpoly-web SQLite DB from schema v1 to v2.",
    )
    parser.add_argument("db_path", type=Path, help="Path to the SQLite DB file.")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Run the migration but roll back the transaction.",
    )
    args = parser.parse_args(argv)

    if not args.db_path.exists():
        print(f"error: DB file does not exist: {args.db_path}", file=sys.stderr)
        return 2

    from regpoly_web.config import SCHEMA_PATH
    from regpoly_web.database import SCHEMA_VERSION
    from regpoly_web.migrations.v2 import migrate_v2_inplace

    conn = sqlite3.connect(str(args.db_path))
    try:
        conn.execute("PRAGMA foreign_keys=ON")
        # Apply latest schema (CREATE IF NOT EXISTS for v2 tables).
        conn.executescript(Path(SCHEMA_PATH).read_text())
        row = conn.execute("SELECT MAX(version) FROM schema_version").fetchone()
        current = row[0] if row is not None else None
        print(f"current schema_version: {current}")
        if current is not None and current >= SCHEMA_VERSION:
            print(f"already at v{SCHEMA_VERSION} or newer — nothing to do.")
            return 0

        counters = migrate_v2_inplace(conn)
        print("backfill counters:")
        for k, v in counters.items():
            print(f"  {k}: {v}")

        if args.dry_run:
            conn.rollback()
            print("--dry-run: rolled back.")
            return 0

        conn.execute(
            "INSERT INTO schema_version(version) VALUES (?)", (SCHEMA_VERSION,)
        )
        conn.commit()
        print(f"committed; schema_version is now v{SCHEMA_VERSION}.")
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
