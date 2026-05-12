# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 5.5 — Playwright e2e fixtures.

Spawns the FastAPI app via uvicorn in a background thread bound to a
free port, pre-seeds an ephemeral SQLite DB with the synthetic data
the e2e tests need, and yields a `live_server_url` fixture the tests
use as their base URL.

Real searches/optimizers are NOT exercised here — they take seconds
to minutes and are covered by the C++ ctest + Python pytest layers.
The e2e suite is about UI/route plumbing on shapes that already exist.

To run:
    uv sync --group e2e
    uv run playwright install chromium
    uv run pytest -m e2e packages/regpoly-web/tests/e2e/
"""

from __future__ import annotations

import os
import socket
import sys
import threading
import time
from collections.abc import Iterator

import httpx
import pytest
import uvicorn


# ── Teardown watchdog ────────────────────────────────────────────────
#
# Known issue: pytest-playwright's session-scoped sync greenlet
# deadlocks in selector.poll during teardown when sharing a process
# with our uvicorn-thread + pgserver. Cleanly avoiding this would
# require either rewriting the fixture stack or upstream fixes; we've
# tried function-scoped browsers and workflow-level test splits and
# neither moved the deadlock — it triggers cumulatively on any
# sufficiently large e2e session.
#
# Pragmatic mitigation: a watchdog thread started at conftest import
# tracks last-test-progress and force-exits with the recorded
# testsfailed count if no progress arrives for 180 s. CI gets the
# real verdict instead of pytest-timeout's os._exit(1) trashing it.
# pytest_sessionstart in a SUBDIRECTORY conftest fires after pytest's
# own session has started, so the watchdog must start at module
# import to be armed before the first test.

_LAST_PROGRESS = [time.monotonic()]
_SESSION_REF: list[pytest.Session] = []
_NO_PROGRESS_BUDGET_S = 180.0


def _watchdog_loop() -> None:
    while True:
        time.sleep(15)
        idle = time.monotonic() - _LAST_PROGRESS[0]
        if idle > _NO_PROGRESS_BUDGET_S:
            session = _SESSION_REF[0] if _SESSION_REF else None
            failed = getattr(session, "testsfailed", 1) if session else 1
            sys.stderr.write(
                f"\n[watchdog] no test progress for >{_NO_PROGRESS_BUDGET_S}s "
                f"(idle={idle:.0f}s); force-exiting with testsfailed={failed}\n"
            )
            sys.stderr.flush()
            os._exit(1 if failed else 0)


threading.Thread(target=_watchdog_loop, daemon=True).start()


def pytest_collection_finish(session: pytest.Session) -> None:
    _SESSION_REF.append(session)
    _LAST_PROGRESS[0] = time.monotonic()


def pytest_runtest_logreport(report: pytest.TestReport) -> None:
    _LAST_PROGRESS[0] = time.monotonic()


def pytest_sessionfinish(session, exitstatus):
    """Fast-path: skip whatever else is still in the interpreter."""
    os._exit(int(exitstatus))


def _find_free_port() -> int:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


# NOTE: P2 prep — the `seeded_db` fixture moved to the root
# tests/conftest.py so non-e2e contract tests can consume it. We
# re-export it here so existing e2e tests that pin the live server
# to it keep working without code change.


@pytest.fixture(scope="session")
def live_server_url(seeded_db: str) -> Iterator[str]:
    """Run uvicorn in a daemon thread on a free port; yield the base URL."""
    from regpoly_web.app import create_app
    from regpoly_web.config import Settings

    settings = Settings(db_url=seeded_db, reload=False, pool_size=1)
    app = create_app(settings)
    port = _find_free_port()
    config = uvicorn.Config(
        app, host="127.0.0.1", port=port, log_level="warning",
        loop="asyncio",
    )
    server = uvicorn.Server(config)
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()

    base_url = f"http://127.0.0.1:{port}"
    deadline = time.time() + 15
    while time.time() < deadline:
        try:
            r = httpx.get(f"{base_url}/api/families", timeout=1.0)
            if r.status_code == 200:
                break
        except httpx.HTTPError:
            time.sleep(0.1)
    else:
        raise RuntimeError("uvicorn did not become ready in 15s")

    try:
        yield base_url
    finally:
        # `force_exit` skips connection drain so idle SSE keepalives
        # don't stall teardown indefinitely.
        server.should_exit = True
        server.force_exit = True
        thread.join(timeout=2)


# pytest-playwright already provides `page`, `browser`, `context`
# fixtures. The only customisation we want is to point new pages at the
# live server's base URL by default.
@pytest.fixture(scope="session")
def base_url(live_server_url: str) -> str:
    return live_server_url
