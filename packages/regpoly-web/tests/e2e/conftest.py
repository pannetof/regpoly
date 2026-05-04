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

import json
import socket
import sqlite3
import threading
import time
from collections.abc import Iterator
from pathlib import Path

import httpx
import pytest
import uvicorn


def _find_free_port() -> int:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


@pytest.fixture(scope="session")
def e2e_db_path(tmp_path_factory: pytest.TempPathFactory) -> str:
    """One isolated SQLite file shared across the e2e session."""
    p = tmp_path_factory.mktemp("e2e_db") / "regpoly_e2e.db"
    return str(p)


@pytest.fixture(scope="session")
def seeded_db(e2e_db_path: str) -> str:
    """Initialize the schema and seed the rows the 8 paths need."""
    from regpoly_web.database import init_sync

    init_sync(e2e_db_path)

    conn = sqlite3.connect(e2e_db_path)
    conn.row_factory = sqlite3.Row
    try:
        conn.execute(
            "INSERT INTO primitive_search_run"
            "(family, L, k, structural_params, fixed_params, "
            " status, tries_done, found_count, elapsed_seconds) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ("MTGen", 19937, 32, json.dumps({}), json.dumps({}),
             "completed", 100, 1, 0.5),
        )
        psr_id = conn.execute(
            "SELECT id FROM primitive_search_run ORDER BY id DESC LIMIT 1"
        ).fetchone()[0]

        conn.execute(
            "INSERT INTO primitive_generator"
            "(search_run_id, family, L, k, structural_params, "
            " search_params, all_params, found_at_try, char_poly, "
            " hamming_weight, pis_se, pis_computed_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))",
            (psr_id, "MTGen", 19937, 32,
             json.dumps({"L": 19937}),
             json.dumps({"a": "0x9908b0df"}),
             json.dumps({"L": 19937, "a": "0x9908b0df"}),
             1, "0xdeadbeef", 135, 0),
        )

        conn.execute(
            "INSERT INTO tested_generator"
            "(id, search_run_id, Lmax, k_g, J) VALUES (?, ?, ?, ?, ?)",
            (4242, None, 32, 19937, 1),
        )
        conn.execute(
            "INSERT INTO tested_generator_component"
            "(tested_gen_id, component_index, family, L, k, "
            " all_params, tempering_params) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (4242, 0, "MTGen", 19937, 32,
             json.dumps({"L": 19937}), json.dumps([])),
        )
        # One typed equidistribution result so the tested-generator detail
        # page has something to render.
        ecart = [0] * 33
        ecart[5] = 1
        ecart[12] = 2
        conn.execute(
            "INSERT INTO equidistribution_result"
            "(tested_gen_id, test_config, kg, L, Lmax, ecart_json, se, "
            " verified, elapsed_seconds) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (4242, json.dumps({"type": "equidistribution", "L": 32}),
             19937, 32, 32, json.dumps(ecart), 3, 1, 0.1),
        )

        # NOTE: tempering_search_run rows are NOT seeded here because
        # `create_app -> init_sync` runs an orphan-reap that flips any
        # 'running'/'pending' rows to 'cancelled' at startup. Tests
        # that need a 'running' row insert it themselves AFTER the
        # server is up.
        conn.commit()
    finally:
        conn.close()
    return e2e_db_path


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
