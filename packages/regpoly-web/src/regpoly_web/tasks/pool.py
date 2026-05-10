# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""ProcessPoolExecutor wrapper for running long searches in worker processes.

Each search task receives the PG DSN and its search_run_id. It opens
its own psycopg connection, runs the work, and writes results +
progress back to the shared database.

Cancellation is cooperative: the API sets `status='cancelled'` on the
run row; the worker checks periodically and exits cleanly.
"""

from __future__ import annotations

import logging
from concurrent.futures import Future, ProcessPoolExecutor
from typing import Callable

logger = logging.getLogger(__name__)


class TaskPool:
    """Thin wrapper over ProcessPoolExecutor tracking live futures."""

    def __init__(self, db_url: str, max_workers: int = 4) -> None:
        # Parameter renamed from db_path → db_url with the SQLite → PG
        # cutover (dockerize-plan Phase 1.2). The wrapped function
        # signature (db_url, run_id, *args, **kwargs) is unchanged
        # because workers always receive the connection target as the
        # first positional argument.
        self.db_url = db_url
        self.executor = ProcessPoolExecutor(max_workers=max_workers)
        self._futures: dict[tuple[str, int], Future] = {}

    def submit(
        self,
        kind: str,
        run_id: int,
        func: Callable[..., None],
        *args,
        **kwargs,
    ) -> Future:
        """Submit a background task.

        `kind` is a free-form label ("primitive", "tempering") used only
        for bookkeeping.
        """
        fut = self.executor.submit(func, self.db_url, run_id, *args, **kwargs)
        key = (kind, run_id)
        self._futures[key] = fut

        def _cleanup(_f: Future, key=key) -> None:
            self._futures.pop(key, None)
            if _f.exception() is not None:
                logger.exception(
                    "Task %s/%d failed", key[0], key[1], exc_info=_f.exception()
                )

        fut.add_done_callback(_cleanup)
        return fut

    def is_running(self, kind: str, run_id: int) -> bool:
        fut = self._futures.get((kind, run_id))
        return fut is not None and not fut.done()

    def shutdown(self, *, wait: bool = False, timeout: float = 0.0) -> None:
        """Shut the pool down.

        - ``wait=False`` (default): immediately SIGKILL every child so
          the port releases instantly. Suitable for dev-mode app
          shutdown where the caller has already marked running rows
          ``status='cancelled'``.
        - ``wait=True, timeout=N``: wait up to ``N`` seconds for
          in-flight jobs to finish naturally, then SIGKILL whatever is
          still alive. Used by the worker container's signal handler
          so a clean ``docker compose stop --timeout 70`` lets running
          C++ jobs commit their results before the kill.
        """
        procs = list(self.executor._processes.values())
        if wait and timeout > 0:
            import time as _time
            self.executor.shutdown(wait=False, cancel_futures=True)
            deadline = _time.monotonic() + timeout
            for proc in procs:
                remaining = max(0.0, deadline - _time.monotonic())
                proc.join(timeout=remaining)
        else:
            self.executor.shutdown(wait=False, cancel_futures=True)
        for proc in procs:
            if proc.is_alive():
                proc.kill()
        for proc in procs:
            try:
                proc.join(timeout=0.5)
            except Exception:
                pass
