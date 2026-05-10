# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Worker container entrypoint.

Single asyncio process that owns:

- One :class:`~regpoly_web.tasks.pool.TaskPool` for search and
  library-test dispatch and one for analysis (sized via
  ``REGPOLY_POOL_SIZE`` and ``REGPOLY_ANALYSIS_POOL_SIZE``).
- A :func:`~regpoly_web.scheduler.search_scheduler` task that polls
  ``primitive_search_run`` / ``tempering_search_run`` /
  ``library_test_run`` for ``status='pending'`` rows and dispatches
  them to the search pool.
- A :func:`~regpoly_web.scheduler.analysis_scheduler` task for
  ``primitive_generator`` PIS computation.
- A heartbeat that touches ``/tmp/worker.alive`` every 10 s so the
  container's healthcheck can detect a wedged scheduler.
- SIGTERM/SIGINT handlers that stop scheduling, mark every locally-
  ``running`` row as ``cancelled`` (single-worker invariant), drain
  the pools, and exit cleanly.
"""

from __future__ import annotations

import argparse
import asyncio
import logging
import signal
import sys
from pathlib import Path

from regpoly_web.config import Settings
from regpoly_web.database import open_pool
from regpoly_web.scheduler import analysis_scheduler, search_scheduler
from regpoly_web.tasks.pool import TaskPool

logger = logging.getLogger(__name__)

_HEARTBEAT_PATH = Path("/tmp/worker.alive")
_HEARTBEAT_INTERVAL = 10.0


async def _heartbeat() -> None:
    try:
        while True:
            try:
                _HEARTBEAT_PATH.touch()
            except OSError:
                logger.exception("failed to touch heartbeat file")
            await asyncio.sleep(_HEARTBEAT_INTERVAL)
    except asyncio.CancelledError:
        pass


async def _cancel_running(dbpool, db_url: str) -> int:
    """Mark every still-``running`` row as ``cancelled``.

    Single-worker invariant: this worker is the only process that
    could own them. When scaling to multi-worker this needs to narrow
    to "rows owned by THIS worker_id" via a heartbeat table.
    """
    total = 0
    for table in ("primitive_search_run", "tempering_search_run",
                  "library_test_run"):
        async with dbpool.connection() as conn:
            async with conn.cursor() as cur:
                await cur.execute(
                    f"UPDATE {table} "
                    f"SET status='cancelled', "
                    f"    error_message=COALESCE(error_message,'shutdown'), "
                    f"    updated_at=NOW() "
                    f"WHERE status='running'"
                )
                total += cur.rowcount
    return total


async def _amain() -> int:
    settings = Settings.from_env()
    logger.info("worker starting; db=%s pool=%d analysis=%d",
                settings.db_url.split("@")[-1],
                settings.pool_size, settings.analysis_pool_size)

    dbpool = await open_pool(
        settings.db_url,
        min_size=settings.db_pool_min_size,
        max_size=settings.db_pool_max_size,
    )
    search_pool = TaskPool(
        db_url=settings.db_url, max_workers=settings.pool_size,
    )
    analysis_pool = TaskPool(
        db_url=settings.db_url, max_workers=settings.analysis_pool_size,
    )

    stop_event = asyncio.Event()

    def _request_stop(signame: str) -> None:
        if not stop_event.is_set():
            logger.info("received %s; draining", signame)
            stop_event.set()

    loop = asyncio.get_running_loop()
    for sig in (signal.SIGTERM, signal.SIGINT):
        try:
            loop.add_signal_handler(sig, _request_stop, sig.name)
        except NotImplementedError:  # pragma: no cover -- Windows
            pass

    tasks = [
        asyncio.create_task(search_scheduler(
            dbpool, search_pool, settings.db_url,
            poll_seconds=settings.worker_poll_seconds,
        )),
        asyncio.create_task(analysis_scheduler(
            dbpool, analysis_pool, settings.db_url,
        )),
        asyncio.create_task(_heartbeat()),
    ]

    await stop_event.wait()
    for t in tasks:
        t.cancel()
    await asyncio.gather(*tasks, return_exceptions=True)

    n = await _cancel_running(dbpool, settings.db_url)
    if n:
        logger.info("cancelled %d in-flight rows", n)

    # Drain the pools so the C++ children get a chance to exit cleanly.
    search_pool.shutdown()
    analysis_pool.shutdown()
    await dbpool.close()
    logger.info("worker exit")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="regpoly_web.worker")
    parser.parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    try:
        return asyncio.run(_amain())
    except Exception as e:
        logger.exception("worker crashed: %s", e)
        return 1


if __name__ == "__main__":
    sys.exit(main())
