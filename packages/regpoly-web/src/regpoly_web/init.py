# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Init container entrypoint.

Waits for PostgreSQL to be reachable, applies any pending schema
migrations, and (with ``--reap-orphans``) cancels any rows left in
``status='running'`` from a prior process. The reap is opt-in so an
operator can run ``python -m regpoly_web.init`` against a live
database for an ad-hoc schema check without disturbing in-flight
work; the compose stack passes ``--reap-orphans`` because it has
already drained the worker container via
``docker compose stop --timeout 70 worker`` before bringing init up.
"""

from __future__ import annotations

import argparse
import asyncio
import logging
import sys

from regpoly_web.config import Settings
from regpoly_web.database import init_db, reap_orphans

logger = logging.getLogger(__name__)


async def _wait_for_pg(db_url: str, budget_seconds: int = 60) -> None:
    """Poll ``db_url`` until ``SELECT 1`` succeeds or the budget elapses."""
    import psycopg

    delay = 0.5
    deadline = asyncio.get_event_loop().time() + budget_seconds
    while True:
        try:
            async with await psycopg.AsyncConnection.connect(
                db_url, connect_timeout=2
            ) as c:
                await c.execute("SELECT 1")
            return
        except psycopg.OperationalError as e:
            if asyncio.get_event_loop().time() > deadline:
                raise RuntimeError(
                    f"PG not reachable after {budget_seconds}s: {e}"
                ) from e
            await asyncio.sleep(delay)
            delay = min(delay * 1.5, 5.0)


async def _amain(reap: bool) -> None:
    settings = Settings.from_env()
    logger.info("waiting for PG at %s", _redact(settings.db_url))
    await _wait_for_pg(settings.db_url, budget_seconds=60)
    logger.info("applying migrations")
    await init_db(settings.db_url)
    if reap:
        # Cold-start invariant: deploy.sh has drained workers via
        # `docker compose stop --timeout 70 worker`, so any remaining
        # 'running' row is definitionally orphaned.
        n = await reap_orphans(settings.db_url, stale_seconds=0)
        logger.info("reaped %d orphan rows", n)
    logger.info("done")


def _redact(dsn: str) -> str:
    """Hide the password component when logging the DSN."""
    import re
    return re.sub(r"(:)([^:@/]+)(@)", r"\1***\3", dsn)


def main() -> int:
    parser = argparse.ArgumentParser(prog="regpoly_web.init")
    parser.add_argument(
        "--reap-orphans", action="store_true",
        help="Cancel any 'running' rows (use only after worker drain).",
    )
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    try:
        asyncio.run(_amain(reap=args.reap_orphans))
        return 0
    except Exception as e:
        logger.exception("init failed: %s", e)
        return 1


if __name__ == "__main__":
    sys.exit(main())
