# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Initial PG schema (dockerize-plan Phase 1).

PostgreSQL 16 translation of the legacy SQLite ``schema.sql``.

Mapping rules applied:

- ``INTEGER PRIMARY KEY AUTOINCREMENT``  → ``BIGSERIAL PRIMARY KEY``
- ``TEXT NOT NULL DEFAULT (datetime('now'))``
                                          → ``TIMESTAMPTZ NOT NULL DEFAULT NOW()``
- ``TEXT`` for JSON payloads              → ``JSONB``
- ``INTEGER`` for boolean flags
  (``verified``, ``is_me``, ``is_cf``)    → ``BOOLEAN``
- ``REAL``                                → ``DOUBLE PRECISION``

Also adds ``updated_at`` columns to ``primitive_search_run``,
``tempering_search_run`` and ``library_test_run`` so the worker's
``reap_orphans`` stale check has something meaningful to compare.
"""

from __future__ import annotations

import psycopg

VERSION = 1

SQL = r"""
CREATE TABLE IF NOT EXISTS schema_version (
    version     INTEGER NOT NULL,
    applied_at  TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- ================================================================
-- PRIMITIVE GENERATORS
-- ================================================================

CREATE TABLE IF NOT EXISTS primitive_search_run (
    id                BIGSERIAL PRIMARY KEY,
    family            TEXT    NOT NULL,
    "L"               INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    structural_params JSONB   NOT NULL,
    fixed_params      JSONB   NOT NULL,
    max_tries         INTEGER,
    max_seconds       DOUBLE PRECISION,
    max_cost          INTEGER,
    status            TEXT    NOT NULL DEFAULT 'pending',
    tries_done        INTEGER NOT NULL DEFAULT 0,
    found_count       INTEGER NOT NULL DEFAULT 0,
    elapsed_seconds   DOUBLE PRECISION,
    error_message     TEXT,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at        TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    started_at        TIMESTAMPTZ,
    finished_at       TIMESTAMPTZ,
    -- Exhaustive-search mode bookkeeping.
    search_mode       TEXT    NOT NULL DEFAULT 'random',
    enum_index        BIGINT  NOT NULL DEFAULT 0,
    enum_total        TEXT,
    enum_axes         JSONB
);

CREATE INDEX IF NOT EXISTS idx_psr_family ON primitive_search_run(family);
CREATE INDEX IF NOT EXISTS idx_psr_status ON primitive_search_run(status);

CREATE TABLE IF NOT EXISTS primitive_generator (
    id                BIGSERIAL PRIMARY KEY,
    search_run_id     BIGINT REFERENCES primitive_search_run(id) ON DELETE SET NULL,
    family            TEXT    NOT NULL,
    "L"               INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    structural_params JSONB   NOT NULL,
    search_params     JSONB   NOT NULL,
    all_params        JSONB   NOT NULL,
    found_at_try      INTEGER,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    -- Analysis fields populated asynchronously after insert.
    char_poly         TEXT,
    hamming_weight    INTEGER,
    pis_gaps          JSONB,
    pis_se            INTEGER,
    pis_elapsed       DOUBLE PRECISION,
    pis_computed_at   TIMESTAMPTZ,
    pis_error         TEXT,
    library_id        TEXT,
    UNIQUE(family, structural_params, search_params)
);

CREATE INDEX IF NOT EXISTS idx_pg_family ON primitive_generator(family);
CREATE INDEX IF NOT EXISTS idx_pg_search_run ON primitive_generator(search_run_id);
CREATE INDEX IF NOT EXISTS idx_pg_k ON primitive_generator(k);
CREATE INDEX IF NOT EXISTS idx_pg_pending_analysis
    ON primitive_generator(pis_computed_at, pis_error);
CREATE UNIQUE INDEX IF NOT EXISTS idx_pg_library_id
    ON primitive_generator(library_id) WHERE library_id IS NOT NULL;

-- ================================================================
-- TEMPERING SEARCH
-- ================================================================

CREATE TABLE IF NOT EXISTS tempering_search_run (
    id                BIGSERIAL PRIMARY KEY,
    test_type         TEXT    NOT NULL,
    test_config       JSONB   NOT NULL,
    "Lmax"            INTEGER NOT NULL,
    nb_tries          INTEGER NOT NULL,
    optimizer_config  JSONB,
    status            TEXT    NOT NULL DEFAULT 'pending',
    combos_total      INTEGER,
    combos_done       INTEGER NOT NULL DEFAULT 0,
    best_se           INTEGER,
    elapsed_seconds   DOUBLE PRECISION,
    error_message     TEXT,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at        TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    started_at        TIMESTAMPTZ,
    finished_at       TIMESTAMPTZ
);

CREATE INDEX IF NOT EXISTS idx_tsr_status ON tempering_search_run(status);

CREATE TABLE IF NOT EXISTS tempering_search_component (
    id                BIGSERIAL PRIMARY KEY,
    search_run_id     BIGINT NOT NULL REFERENCES tempering_search_run(id) ON DELETE CASCADE,
    component_index   INTEGER NOT NULL,
    shared_with_component INTEGER,
    tempering_config  JSONB   NOT NULL,
    UNIQUE(search_run_id, component_index)
);

CREATE TABLE IF NOT EXISTS tempering_search_generator (
    id                BIGSERIAL PRIMARY KEY,
    component_id      BIGINT NOT NULL REFERENCES tempering_search_component(id) ON DELETE CASCADE,
    generator_id      BIGINT NOT NULL REFERENCES primitive_generator(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_tsg_component ON tempering_search_generator(component_id);

-- ================================================================
-- TESTED GENERATORS
-- ================================================================

CREATE TABLE IF NOT EXISTS tested_generator (
    id                BIGSERIAL PRIMARY KEY,
    search_run_id     BIGINT REFERENCES tempering_search_run(id) ON DELETE SET NULL,
    "Lmax"            INTEGER NOT NULL,
    k_g               INTEGER NOT NULL,
    "J"               INTEGER NOT NULL,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    library_id        TEXT
);

CREATE INDEX IF NOT EXISTS idx_tg_search_run ON tested_generator(search_run_id);
CREATE INDEX IF NOT EXISTS idx_tg_k_g ON tested_generator(k_g);
CREATE UNIQUE INDEX IF NOT EXISTS idx_tg_library_id
    ON tested_generator(library_id) WHERE library_id IS NOT NULL;

CREATE TABLE IF NOT EXISTS tested_generator_component (
    id                BIGSERIAL PRIMARY KEY,
    tested_gen_id     BIGINT NOT NULL REFERENCES tested_generator(id) ON DELETE CASCADE,
    component_index   INTEGER NOT NULL,
    generator_id      BIGINT REFERENCES primitive_generator(id) ON DELETE SET NULL,
    family            TEXT    NOT NULL,
    "L"               INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    all_params        JSONB   NOT NULL,
    tempering_params  JSONB   NOT NULL,
    UNIQUE(tested_gen_id, component_index)
);

CREATE TABLE IF NOT EXISTS test_result (
    id                BIGSERIAL PRIMARY KEY,
    tested_gen_id     BIGINT NOT NULL REFERENCES tested_generator(id) ON DELETE CASCADE,
    test_type         TEXT    NOT NULL,
    test_config       JSONB   NOT NULL,
    se                INTEGER,
    is_me             BOOLEAN,
    secf              INTEGER,
    is_cf             BOOLEAN,
    score             DOUBLE PRECISION,
    detail            JSONB   NOT NULL,
    elapsed_seconds   DOUBLE PRECISION,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_tr_tested_gen ON test_result(tested_gen_id);
CREATE INDEX IF NOT EXISTS idx_tr_test_type ON test_result(test_type);
CREATE INDEX IF NOT EXISTS idx_tr_se ON test_result(se);

-- ================================================================
-- TYPED TEST RESULTS (schema v2)
-- ================================================================

CREATE TABLE IF NOT EXISTS equidistribution_result (
    id                BIGSERIAL PRIMARY KEY,
    tested_gen_id     BIGINT NOT NULL REFERENCES tested_generator(id) ON DELETE CASCADE,
    test_config       JSONB   NOT NULL,
    kg                INTEGER NOT NULL,
    "L"               INTEGER NOT NULL,
    "Lmax"            INTEGER NOT NULL,
    ecart_json        JSONB   NOT NULL,
    se                INTEGER NOT NULL,
    verified          BOOLEAN NOT NULL,
    elapsed_seconds   DOUBLE PRECISION,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_er_tested_gen ON equidistribution_result(tested_gen_id);
CREATE INDEX IF NOT EXISTS idx_er_se ON equidistribution_result(se);

CREATE TABLE IF NOT EXISTS collision_free_result (
    id                BIGSERIAL PRIMARY KEY,
    tested_gen_id     BIGINT NOT NULL REFERENCES tested_generator(id) ON DELETE CASCADE,
    test_config       JSONB   NOT NULL,
    kg                INTEGER NOT NULL,
    "L"               INTEGER NOT NULL,
    ecart_cf_json     JSONB   NOT NULL,
    secf              INTEGER NOT NULL,
    verified          BOOLEAN NOT NULL,
    elapsed_seconds   DOUBLE PRECISION,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_cfr_tested_gen ON collision_free_result(tested_gen_id);
CREATE INDEX IF NOT EXISTS idx_cfr_secf ON collision_free_result(secf);

CREATE TABLE IF NOT EXISTS tuplets_result (
    id                BIGSERIAL PRIMARY KEY,
    tested_gen_id     BIGINT NOT NULL REFERENCES tested_generator(id) ON DELETE CASCADE,
    test_config       JSONB   NOT NULL,
    kg                INTEGER NOT NULL,
    "L"               INTEGER NOT NULL,
    tupd              INTEGER NOT NULL,
    testtype          INTEGER NOT NULL,
    threshold         DOUBLE PRECISION NOT NULL,
    tuph_json         JSONB   NOT NULL,
    gap_json          JSONB   NOT NULL,
    "DELTA_json"      JSONB   NOT NULL,
    pourcentage_json  JSONB   NOT NULL,
    firstpart_max     DOUBLE PRECISION NOT NULL,
    firstpart_sum     DOUBLE PRECISION NOT NULL,
    secondpart_max    DOUBLE PRECISION NOT NULL,
    secondpart_sum    DOUBLE PRECISION NOT NULL,
    elapsed_seconds   DOUBLE PRECISION,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_tplr_tested_gen ON tuplets_result(tested_gen_id);

-- ================================================================
-- SEARCH PROGRESS (live updates)
-- ================================================================

CREATE TABLE IF NOT EXISTS search_progress (
    id                BIGSERIAL PRIMARY KEY,
    search_type       TEXT    NOT NULL,
    search_run_id     BIGINT  NOT NULL,
    tries_done        INTEGER NOT NULL DEFAULT 0,
    found_count       INTEGER NOT NULL DEFAULT 0,
    current_info      TEXT,
    message           TEXT,
    updated_at        TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_sp_lookup ON search_progress(search_type, search_run_id);

-- ================================================================
-- YAML IMPORT TRACKING
-- ================================================================

CREATE TABLE IF NOT EXISTS yaml_import (
    id                BIGSERIAL PRIMARY KEY,
    file_path         TEXT    NOT NULL UNIQUE,
    import_type       TEXT    NOT NULL,
    row_count         INTEGER NOT NULL,
    imported_at       TIMESTAMPTZ NOT NULL DEFAULT NOW()
);
"""


def apply(conn: psycopg.Connection) -> None:
    """Apply the initial PG schema. Idempotent (every CREATE uses IF NOT EXISTS)."""
    conn.execute(SQL)
