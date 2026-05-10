# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Unauthenticated `/healthz` probe (dockerize-plan Phase 2).

External monitors (UptimeRobot) and the Docker container healthcheck
both hit this endpoint. The Caddyfile exempts `/healthz` from basic
auth via an exact-path matcher, so credentials never need to be
shared with a third-party service.
"""

from __future__ import annotations

import asyncio

import psycopg
from fastapi import APIRouter, Request, Response

router = APIRouter()


@router.get("/healthz", response_model=None)
async def healthz(request: Request) -> dict | Response:
    """Return 200 with ``{"status":"ok","db":"up"}`` if the asyncpg pool
    can do a ``SELECT 1`` within 2 seconds. 503 with
    ``{"status":"degraded","db":"down"}`` otherwise. The pool-acquire
    timeout prevents a saturated pool from hanging the healthcheck.
    """
    pool = request.app.state.dbpool
    if pool is None:
        return {"status": "ok", "db": "skipped"}
    try:
        async with asyncio.timeout(2):
            async with pool.connection() as conn:
                async with conn.cursor() as cur:
                    await cur.execute("SELECT 1")
                    await cur.fetchone()
        return {"status": "ok", "db": "up"}
    except (asyncio.TimeoutError, psycopg.OperationalError):
        return Response(
            content='{"status":"degraded","db":"down"}',
            status_code=503,
            media_type="application/json",
        )
