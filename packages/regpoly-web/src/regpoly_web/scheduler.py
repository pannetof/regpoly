# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Pending-row dispatcher used by the web (dev mode) and worker (prod).

Routes only ever insert ``status='pending'`` rows. Dispatch happens in
this module: a periodic poll uses ``UPDATE ... FOR UPDATE SKIP LOCKED``
to atomically claim a batch of pending rows and bumps them to
``status='running'``. Each claimed row is then handed to a
:class:`regpoly_web.tasks.pool.TaskPool` for the actual computation.

Two schedulers live here:

- :func:`search_scheduler` — polls ``primitive_search_run`` /
  ``tempering_search_run`` / ``library_test_run``.
- :func:`analysis_scheduler` — polls ``primitive_generator`` for rows
  needing PIS computation.

Both are coroutines designed to run as ``asyncio.create_task(...)``;
they exit cleanly on ``CancelledError``.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Callable

from psycopg_pool import AsyncConnectionPool

from regpoly_web.tasks.analysis import analyze_generator
from regpoly_web.tasks.library_test import run_library_test
from regpoly_web.tasks.pool import TaskPool
from regpoly_web.tasks.primitive import run_primitive_search
from regpoly_web.tasks.tempering import run_tempering_search

logger = logging.getLogger(__name__)


_DISPATCH_TABLE: tuple[tuple[str, str, Callable], ...] = (
    ("primitive_search_run", "primitive", run_primitive_search),
    ("tempering_search_run", "tempering", run_tempering_search),
    ("library_test_run",     "library_test", run_library_test),
)


async def _claim_pending(
    dbpool: AsyncConnectionPool, table: str, batch: int,
) -> list[int]:
    """Atomically flip up to ``batch`` rows from ``pending`` →
    ``running``; return the claimed ids.

    Two callers (web's dev-mode embed + a worker container) can race on
    the same DB safely thanks to ``FOR UPDATE SKIP LOCKED``.
    """
    sql = f"""
        UPDATE {table}
           SET status='running', updated_at=NOW()
         WHERE id IN (
             SELECT id FROM {table}
              WHERE status='pending'
              ORDER BY id
              LIMIT %s
              FOR UPDATE SKIP LOCKED
         )
         RETURNING id
    """
    async with dbpool.connection() as conn:
        async with conn.cursor() as cur:
            await cur.execute(sql, (batch,))
            rows = await cur.fetchall()
    return [int(r[0]) for r in rows]


async def _present_tables(
    dbpool: AsyncConnectionPool, candidates: list[str],
) -> set[str]:
    async with dbpool.connection() as conn:
        async with conn.cursor() as cur:
            await cur.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = current_schema() "
                "AND table_name = ANY(%s)",
                (candidates,),
            )
            rows = await cur.fetchall()
    return {r[0] for r in rows}


async def search_scheduler(
    dbpool: AsyncConnectionPool,
    pool: TaskPool,
    db_url: str,
    *,
    poll_seconds: float = 0.5,
    batch: int = 5,
) -> None:
    """Dispatch pending search/library-test rows to ``pool``.

    For each claimed row id, schedules ``fn(db_url, row_id)`` on the
    :class:`TaskPool`. Loop body is wrapped in a broad except so a
    transient PG hiccup doesn't kill the scheduler.
    """
    present = await _present_tables(
        dbpool, [t for t, _, _ in _DISPATCH_TABLE]
    )
    dispatch = [(t, k, f) for t, k, f in _DISPATCH_TABLE if t in present]
    if not dispatch:
        logger.warning("search_scheduler: no dispatch tables present; idling")
        return
    in_flight: set[tuple[str, int]] = set()
    try:
        while True:
            try:
                for table, kind, fn in dispatch:
                    ids = await _claim_pending(dbpool, table, batch)
                    for rid in ids:
                        key = (kind, rid)
                        if key in in_flight:
                            continue
                        in_flight.add(key)
                        fut = pool.submit(kind, rid, fn)
                        fut.add_done_callback(
                            lambda _f, k=key: in_flight.discard(k)
                        )
            except Exception:
                logger.exception("search_scheduler tick failed")
            await asyncio.sleep(poll_seconds)
    except asyncio.CancelledError:
        pass


async def analysis_scheduler(
    dbpool: AsyncConnectionPool,
    pool: TaskPool,
    db_url: str,
    *,
    poll_seconds: float = 2.0,
    batch: int = 50,
) -> None:
    """Poll ``primitive_generator`` for un-analysed rows and dispatch
    them to ``pool``."""
    in_flight: set[int] = set()
    try:
        while True:
            try:
                async with dbpool.connection() as conn:
                    async with conn.cursor() as cur:
                        await cur.execute(
                            "SELECT id FROM primitive_generator "
                            "WHERE pis_computed_at IS NULL "
                            "ORDER BY id ASC LIMIT %s",
                            (batch,),
                        )
                        rows = await cur.fetchall()
                for r in rows:
                    gid = int(r[0])
                    if gid in in_flight:
                        continue
                    in_flight.add(gid)
                    fut = pool.submit("analysis", gid, analyze_generator)
                    fut.add_done_callback(
                        lambda _f, g=gid: in_flight.discard(g)
                    )
            except Exception:
                logger.exception("analysis_scheduler tick failed")
            await asyncio.sleep(poll_seconds)
    except asyncio.CancelledError:
        pass


