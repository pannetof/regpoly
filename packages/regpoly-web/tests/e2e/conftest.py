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

    settings = Settings(db_path=seeded_db, reload=False, pool_size=1)
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
        server.should_exit = True
        thread.join(timeout=5)


# pytest-playwright already provides `page`, `browser`, `context`
# fixtures. The only customisation we want is to point new pages at the
# live server's base URL by default.
@pytest.fixture(scope="session")
def base_url(live_server_url: str) -> str:
    return live_server_url
