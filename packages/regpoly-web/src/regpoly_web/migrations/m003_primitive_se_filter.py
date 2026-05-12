# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Add SE post-filter columns to ``primitive_search_run``; widen ``pis_se``.

Lets a primitive-search run filter full-period candidates by their
sum-of-gaps (SE = Σ Δ_v) inline in the worker hot loop, keeping only
generators with ``se ≤ max_se``. The ``rejected_count`` column tracks
full-period candidates that failed the SE threshold so the count
survives pause/resume and page reload.

``pis_se`` is widened from INTEGER to BIGINT because SE = Σ gaps can
exceed 2³¹-1 for large (L, k) generators.
"""

from __future__ import annotations

import psycopg

VERSION = 3

SQL = r"""
ALTER TABLE primitive_search_run
    ADD COLUMN IF NOT EXISTS max_se         INTEGER,
    ADD COLUMN IF NOT EXISTS rejected_count BIGINT NOT NULL DEFAULT 0;

ALTER TABLE primitive_generator
    ALTER COLUMN pis_se TYPE BIGINT;
"""


def apply(conn: psycopg.Connection) -> None:
    conn.execute(SQL)
