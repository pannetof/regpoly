"""ProcessPoolExecutor wrapper for running long searches in worker processes.

Each search task receives the SQLite path and its search_run_id.  It opens
its own connection, runs the work, and writes results + progress back to
the shared database.

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

    def __init__(self, db_path: str, max_workers: int = 4) -> None:
        self.db_path = db_path
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
        fut = self.executor.submit(func, self.db_path, run_id, *args, **kwargs)
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

    def shutdown(self) -> None:
        # The caller has already flagged running searches as 'cancelled'
        # in the DB, so worker state is preserved.  Workers running heavy
        # C++ code may not respond to SIGTERM promptly — go straight to
        # SIGKILL so the port releases instantly.
        procs = list(self.executor._processes.values())
        self.executor.shutdown(wait=False, cancel_futures=True)
        for proc in procs:
            if proc.is_alive():
                proc.kill()
        for proc in procs:
            try:
                proc.join(timeout=0.5)
            except Exception:
                pass
