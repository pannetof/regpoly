# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Add ``library_test_run`` table (dockerize-plan Phase 2).

Replaces the in-memory ``app.state.run_test_jobs`` dict that drives the
"Run test" button on library detail pages. Backed by a real DB row so
job state survives web container restarts and so the worker container
can pick up pending jobs from the queue.
"""

from __future__ import annotations

import psycopg

VERSION = 2

SQL = r"""
CREATE TABLE IF NOT EXISTS library_test_run (
    id            BIGSERIAL PRIMARY KEY,
    library_id    TEXT    NOT NULL,
    test_type     TEXT    NOT NULL,
    method        TEXT,
    test_config   JSONB   NOT NULL,
    status        TEXT    NOT NULL DEFAULT 'pending',
    result        JSONB,
    error_code    TEXT,
    error_message TEXT,
    created_at    TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at    TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    finished_at   TIMESTAMPTZ
);

CREATE INDEX IF NOT EXISTS idx_ltr_status_id
    ON library_test_run(status, id)
    WHERE status = 'pending';
"""


def apply(conn: psycopg.Connection) -> None:
    conn.execute(SQL)
