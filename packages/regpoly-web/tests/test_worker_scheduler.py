# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Tests for the search-row dispatcher (regpoly_web.scheduler).

The scheduler claims ``status='pending'`` rows via
``UPDATE … FOR UPDATE SKIP LOCKED`` and dispatches the matching worker
function on the supplied :class:`TaskPool`. The dockerize plan uses
this in two places — the worker container's own loop (worker.py) and
the dev-mode embed inside the FastAPI lifespan (app.py).

These tests exercise the scheduler against a real ephemeral PG but
patch the dispatch table to a no-op function so the real
multi-second C++ search loop doesn't enter the picture.
"""

from __future__ import annotations

import asyncio
import json
import os
from concurrent.futures import Future
from pathlib import Path

import psycopg
import pytest

from regpoly_web.database import open_pool
from regpoly_web.tasks.pool import TaskPool


def _noop_worker(db_url: str, run_id: int) -> None:
    """Stand-in for run_primitive_search.

    Marks the claimed row 'completed' synchronously so we can assert
    the full pending → running → completed transition without
    dragging the real search loop in.

    Also bumps ``found_count`` so multi-dispatch tests can detect
    double-firing — each call increments ``found_count`` by 1.
    """
    conn = psycopg.connect(db_url, autocommit=True)
    try:
        conn.execute(
            "UPDATE primitive_search_run "
            "SET status='completed', "
            "    found_count = found_count + 1, "
            "    updated_at=NOW(), finished_at=NOW() "
            "WHERE id = %s",
            (run_id,),
        )
    finally:
        conn.close()


def _insert_pending(db_url: str) -> int:
    conn = psycopg.connect(db_url, autocommit=True)
    try:
        cur = conn.execute(
            "INSERT INTO primitive_search_run "
            "(family, l, k, structural_params, fixed_params, "
            " max_tries, status) "
            "VALUES (%s, %s, %s, %s, %s, %s, 'pending') "
            "RETURNING id",
            ("MTGen", 64, 32,
             json.dumps({"w": 32, "n": 2, "m": 1, "r": 31, "u": 11}),
             json.dumps({}), 1),
        )
        return int(cur.fetchone()[0])
    finally:
        conn.close()


def _read_status(db_url: str, run_id: int) -> str:
    conn = psycopg.connect(db_url, autocommit=True)
    try:
        row = conn.execute(
            "SELECT status FROM primitive_search_run WHERE id = %s",
            (run_id,),
        ).fetchone()
    finally:
        conn.close()
    assert row is not None
    return row[0]


@pytest.mark.asyncio
async def test_scheduler_claims_and_dispatches_pending_row(
    tmp_db_url: str, monkeypatch,
) -> None:
    """A `pending` row gets claimed (UPDATE → 'running') and dispatched.

    With the no-op worker the row lands at 'completed' within a few
    poll intervals.
    """
    from regpoly_web import scheduler as sched

    # Patch the dispatch table to route 'primitive_search_run' to our
    # no-op worker. Other tables stay attached to their real workers
    # but nothing's pending in them, so the scheduler skips.
    monkeypatch.setattr(
        sched, "_DISPATCH_TABLE",
        (("primitive_search_run", "primitive", _noop_worker),),
    )

    run_id = _insert_pending(tmp_db_url)
    assert _read_status(tmp_db_url, run_id) == "pending"

    dbpool = await open_pool(tmp_db_url, min_size=1, max_size=4)
    pool = TaskPool(db_url=tmp_db_url, max_workers=1)
    task = asyncio.create_task(sched.search_scheduler(
        dbpool, pool, tmp_db_url, poll_seconds=0.05, batch=5,
    ))
    try:
        # Poll up to ~5s for the row to reach 'completed'.
        deadline = asyncio.get_event_loop().time() + 5.0
        while asyncio.get_event_loop().time() < deadline:
            status = _read_status(tmp_db_url, run_id)
            if status == "completed":
                break
            await asyncio.sleep(0.1)
        assert status == "completed", (
            f"row stuck in {status!r} — scheduler did not dispatch"
        )
    finally:
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass
        pool.shutdown()
        await dbpool.close()


@pytest.mark.asyncio
async def test_scheduler_skip_locked_avoids_double_dispatch(
    tmp_db_url: str, monkeypatch,
) -> None:
    """``FOR UPDATE SKIP LOCKED`` must prevent two scheduler instances
    from claiming the same row.

    Simulates the multi-worker future by running TWO scheduler tasks
    against the same DB and counting how many times the no-op fires
    per row.
    """
    from regpoly_web import scheduler as sched

    monkeypatch.setattr(
        sched, "_DISPATCH_TABLE",
        (("primitive_search_run", "primitive", _noop_worker),),
    )

    ids = [_insert_pending(tmp_db_url) for _ in range(3)]

    dbpool = await open_pool(tmp_db_url, min_size=2, max_size=8)
    pool_a = TaskPool(db_url=tmp_db_url, max_workers=1)
    pool_b = TaskPool(db_url=tmp_db_url, max_workers=1)
    task_a = asyncio.create_task(sched.search_scheduler(
        dbpool, pool_a, tmp_db_url, poll_seconds=0.05, batch=5,
    ))
    task_b = asyncio.create_task(sched.search_scheduler(
        dbpool, pool_b, tmp_db_url, poll_seconds=0.05, batch=5,
    ))
    try:
        deadline = asyncio.get_event_loop().time() + 5.0
        while asyncio.get_event_loop().time() < deadline:
            done = sum(
                1 for rid in ids
                if _read_status(tmp_db_url, rid) == "completed"
            )
            if done == 3:
                break
            await asyncio.sleep(0.1)
        assert done == 3, f"only {done}/3 rows completed"

        # Each row's noop_worker bumped found_count by 1 per dispatch.
        # SKIP LOCKED guarantees exactly one dispatch per row across
        # the two schedulers, so found_count must be 1 everywhere.
        # Give a small grace period in case a second scheduler tick
        # is in flight for a stale row.
        await asyncio.sleep(0.5)
        conn = psycopg.connect(tmp_db_url, autocommit=True)
        try:
            rows = conn.execute(
                "SELECT id, found_count FROM primitive_search_run "
                "WHERE id = ANY(%s) ORDER BY id",
                (ids,),
            ).fetchall()
        finally:
            conn.close()
        assert all(r[1] == 1 for r in rows), (
            f"each row must be dispatched exactly once across both "
            f"schedulers; got found_count per row: {dict(rows)}"
        )
    finally:
        for t in (task_a, task_b):
            t.cancel()
            try:
                await t
            except asyncio.CancelledError:
                pass
        pool_a.shutdown()
        pool_b.shutdown()
        await dbpool.close()


@pytest.mark.asyncio
async def test_scheduler_skips_missing_optional_table(
    tmp_db_url: str, monkeypatch,
) -> None:
    """If a dispatch-table entry references a non-existent table, the
    scheduler logs a warning and idles instead of crash-looping.

    Reproduces the dev scenario where someone runs the scheduler
    against an old DB that hasn't had m002_library_test_run applied.
    """
    from regpoly_web import scheduler as sched

    # All three entries point at a fake table that doesn't exist.
    monkeypatch.setattr(
        sched, "_DISPATCH_TABLE",
        (("does_not_exist_run", "fake", _noop_worker),),
    )

    dbpool = await open_pool(tmp_db_url, min_size=1, max_size=2)
    pool = TaskPool(db_url=tmp_db_url, max_workers=1)
    task = asyncio.create_task(sched.search_scheduler(
        dbpool, pool, tmp_db_url, poll_seconds=0.05,
    ))
    try:
        # Give the scheduler a moment to do its present-tables probe
        # and exit early. It should NOT raise.
        await asyncio.sleep(0.5)
        # Task should have exited cleanly (returned, not cancelled).
        assert task.done(), "scheduler must idle-and-exit when no dispatch tables match"
        assert task.exception() is None, (
            f"scheduler must not raise on missing table; got {task.exception()!r}"
        )
    finally:
        if not task.done():
            task.cancel()
            try:
                await task
            except asyncio.CancelledError:
                pass
        pool.shutdown()
        await dbpool.close()
