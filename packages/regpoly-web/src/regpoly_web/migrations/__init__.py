# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Schema migration package for regpoly-web (PostgreSQL).

The dockerize-and-deploy work (commit ``d200a94…``) switched the web
app from SQLite to PostgreSQL 16. Migrations now follow a numbered
``mNNN_description.py`` naming convention; each module exports:

- ``VERSION`` (int): the schema version this migration produces.
- ``apply(conn: psycopg.Connection) -> None``: run the migration.

``database.init_db`` reads ``MAX(version)`` from the ``schema_version``
table, applies any unapplied modules in ascending VERSION order, and
INSERTs each into ``schema_version`` on success. Each ``apply()`` is
idempotent (every CREATE / ALTER uses ``IF NOT EXISTS``).

Migrations live in this package alongside two legacy helpers kept for
the v1 (SQLite) → v2 (PG) cutover:

- ``v2.py`` — original schema_v1→v2 in-place migration on SQLite.
  Not invoked from ``init_db``; kept for reference and for any operator
  who still needs to migrate an in-place SQLite DB to schema v2 before
  exporting via ``scripts/migrate_sqlite_to_pg.py``.
- ``well_matrices.py`` — opt-in JSON-payload migration for the WELL
  matrices format. Operator-run via
  ``packages/regpoly-web/scripts/migrate_well_matrices.py``; never
  auto-applied.

To add a new migration:

1. Create ``mNNN_short_description.py`` with the next VERSION number.
2. Define ``VERSION`` and ``apply(conn)`` at module level.
3. Test against an ephemeral PG (``pytest-postgresql``); commit.
"""
