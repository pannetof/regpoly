"""pytest fixtures for regpoly-web.

Phase 5.3 introduces a FastAPI TestClient fixture backed by an
isolated SQLite file in a tmpdir. Subsequent phases will add a
mocked process-pool fixture and Playwright browser fixtures.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path

import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def tmp_db_path(tmp_path: Path) -> str:
    """An ephemeral SQLite file used by a single web app instance.

    SQLite + aiosqlite require a real file path (not :memory:) because
    the worker pool processes open their own connections to the same
    DB. tmp_path gives each test its own isolated DB.
    """
    return str(tmp_path / "test_regpoly.db")


@pytest.fixture
def web_settings(tmp_db_path: str):
    """Settings for an app rooted at the tmp DB.

    Disables hot-reload so the catalog is loaded once at startup
    (matches the production lifespan).
    """
    from regpoly_web.config import Settings

    return Settings(
        db_path=tmp_db_path,
        reload=False,
        pool_size=1,
    )


@pytest.fixture
def client(web_settings) -> Iterator[TestClient]:
    """A FastAPI TestClient with the app's lifespan executed
    (DB opened, catalog loaded, pool started, scheduler running)."""
    from regpoly_web.app import create_app

    app = create_app(web_settings)
    # `with TestClient(app)` runs the lifespan; without the context
    # manager startup events don't fire and app.state is empty.
    with TestClient(app) as c:
        yield c
