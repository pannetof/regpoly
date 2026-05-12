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

import socket
import threading
import time
from collections.abc import Iterator

import httpx
import pytest
import uvicorn


# ── Upstream bug patch: PipeTransport.request_stop deadlock ──────────
#
# playwright._impl._transport.PipeTransport.request_stop() only closes
# stdin to the node driver subprocess and *hopes* the driver notices
# and exits. The transport's `run()` loop is parked on
#   `await self._proc.stdout.readexactly(4)`
# which only unblocks when stdout closes. If the driver doesn't exit
# (e.g., because a browser context is still draining a long-lived
# connection like our SSE keepalives), stdout never closes, the read
# never unblocks, and pytest-playwright's session teardown wedges in
# selector.poll forever — exactly the stack trace we kept hitting.
#
# Real fix: after closing stdin, schedule a SIGKILL on the driver if
# it doesn't exit within a short grace period. This unblocks the read
# (the OS closes the stdout fd) so the transport loop exits cleanly.
#
# This needs to happen in two places:
#   - on Transport stop (per playwright instance)
#   - one shot at process exit (atexit) as a final safety net
#
# We monkey-patch instead of forking pytest-playwright/playwright-python.

def _install_pipe_transport_kill_patch() -> None:
    import os
    import signal

    try:
        from playwright._impl._transport import PipeTransport
    except ImportError:
        # The non-e2e CI lane syncs only the `test` group and still
        # discovers this conftest during collection. Skip the patch —
        # if playwright isn't installed, no e2e test will actually run.
        return

    if getattr(PipeTransport, "_regpoly_kill_patched", False):
        return

    _original_request_stop = PipeTransport.request_stop

    def _request_stop(self) -> None:
        _original_request_stop(self)
        proc = getattr(self, "_proc", None)
        if proc is None:
            return
        pid = getattr(proc, "pid", None) or getattr(getattr(proc, "_proc", None), "pid", None)
        if pid is None:
            return

        def _kill_after_grace() -> None:
            # 3 s graceful window for the driver to react to stdin EOF.
            time.sleep(3)
            try:
                os.kill(pid, 0)  # still alive?
            except ProcessLookupError:
                return  # exited cleanly
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                return
            time.sleep(2)
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

        threading.Thread(target=_kill_after_grace, daemon=True).start()

    PipeTransport.request_stop = _request_stop  # type: ignore[method-assign]
    PipeTransport._regpoly_kill_patched = True  # type: ignore[attr-defined]


_install_pipe_transport_kill_patch()


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
