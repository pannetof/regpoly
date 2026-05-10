# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""SIGTERM cleanup: ``regpoly_web.worker`` cancels in-flight rows on
graceful shutdown.

The dockerize plan promises that ``docker compose stop --timeout 70
worker`` drains cleanly: scheduled jobs that are still ``running``
when the SIGTERM lands get their status set to ``cancelled`` with
``error_message='shutdown'`` so the front end stops showing them as
live. This test spawns the worker as a subprocess against the
ephemeral PG, plants a ``running`` row, sends SIGTERM, and reads
back the post-shutdown state.

The worker subprocess hits the real `run_primitive_search` C++
path; we keep work small (``max_tries=1``) and short — the test
asserts the FINAL state matches the shutdown contract, not any
particular intermediate timing.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time

import psycopg
import pytest


def _wait_for_status(
    db_url: str, run_id: int, statuses: set[str], deadline: float,
) -> str | None:
    """Poll until the row's status is in ``statuses`` or the deadline
    elapses. Returns the matching status or ``None`` on timeout."""
    while time.time() < deadline:
        conn = psycopg.connect(db_url, autocommit=True)
        try:
            row = conn.execute(
                "SELECT status FROM primitive_search_run WHERE id = %s",
                (run_id,),
            ).fetchone()
        finally:
            conn.close()
        if row and row[0] in statuses:
            return row[0]
        time.sleep(0.1)
    return None


@pytest.mark.slow
def test_worker_sigterm_marks_running_rows_cancelled_with_shutdown(
    tmp_db_url: str,
) -> None:
    """A SIGTERM during a `running` job leaves the row `cancelled`
    with error_message='shutdown'."""
    # Plant a row directly in 'running' state, simulating a worker
    # that has just claimed it but hasn't completed. We could go
    # through `pending` and rely on the worker's scheduler to claim,
    # but that races with the spawn time; planting a pre-claimed row
    # makes the test deterministic.
    conn = psycopg.connect(tmp_db_url, autocommit=True)
    try:
        cur = conn.execute(
            "INSERT INTO primitive_search_run "
            "(family, l, k, structural_params, fixed_params, "
            " max_tries, status, started_at) "
            "VALUES (%s, %s, %s, %s, %s, %s, 'running', NOW()) "
            "RETURNING id",
            ("MTGen", 64, 32,
             json.dumps({"w": 32, "n": 2, "m": 1, "r": 31, "u": 11}),
             json.dumps({}), 100000),
        )
        run_id = int(cur.fetchone()[0])
    finally:
        conn.close()

    env = os.environ.copy()
    env["REGPOLY_DB_URL"] = tmp_db_url
    env["REGPOLY_WORKER"] = "1"
    env["REGPOLY_POOL_SIZE"] = "1"
    env["REGPOLY_ANALYSIS_POOL_SIZE"] = "1"

    proc = subprocess.Popen(
        [sys.executable, "-m", "regpoly_web.worker"],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    try:
        # Wait for the worker's heartbeat to fire — confirms the event
        # loop is up and signal handlers are installed.
        time.sleep(2.0)
        # SIGTERM and wait for clean exit (drain budget is 60s but
        # an unscheduled 'running' row exits immediately).
        proc.send_signal(signal.SIGTERM)
        try:
            proc.wait(timeout=75)
        except subprocess.TimeoutExpired:
            proc.kill()
            raise AssertionError(
                "worker did not exit within 75 s of SIGTERM"
            )
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()

    # Read the final state.
    conn = psycopg.connect(tmp_db_url, autocommit=True)
    try:
        row = conn.execute(
            "SELECT status, error_message FROM primitive_search_run "
            "WHERE id = %s",
            (run_id,),
        ).fetchone()
    finally:
        conn.close()

    assert row is not None
    status, err = row
    assert status == "cancelled", (
        f"row must be cancelled after SIGTERM; got status={status!r}, "
        f"error_message={err!r}"
    )
    assert err == "shutdown", (
        f"shutdown reason must be exactly 'shutdown' (Phase 8 step 10); "
        f"got {err!r}"
    )
