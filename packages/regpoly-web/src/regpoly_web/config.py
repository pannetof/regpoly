# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Runtime configuration for the regpoly web application."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Settings:
    """Web application settings, overridable via environment variables.

    The dockerize-and-deploy migration switched the backing store from
    SQLite to PostgreSQL 16; ``db_url`` is a psycopg DSN
    (``postgresql://user:pw@host:5432/dbname``).
    """

    db_url: str = "postgresql://regpoly@localhost:5432/regpoly"
    host: str = "127.0.0.1"
    port: int = 8000
    pool_size: int = 4
    analysis_pool_size: int | None = None
    progress_poll_seconds: float = 0.5
    reload: bool = False
    # P6 — confine `/api/import/generators-dir` to a known root.
    import_root: Path | None = None
    # Phase 2 worker split: when set, web's lifespan does not create
    # the in-process ProcessPoolExecutor and does not start the
    # _analysis_scheduler. Production compose enables this; dev
    # defaults to single-process behaviour.
    disable_internal_pool: bool = False
    # Worker module flag (informational; toggles a few startup
    # behaviours like the orphan-reap path).
    is_worker: bool = False
    # How often the worker polls the DB for pending searches.
    worker_poll_seconds: float = 0.5
    # asyncpg pool sizing for the web container.
    db_pool_min_size: int = 2
    db_pool_max_size: int = 10

    def __post_init__(self) -> None:
        # Default analysis pool to half of pool_size (rounded up).
        if self.analysis_pool_size is None:
            self.analysis_pool_size = max(1, self.pool_size // 2)

    @classmethod
    def from_env(cls) -> "Settings":
        root_env = os.environ.get("REGPOLY_IMPORT_ROOT")
        pool_size = int(os.environ.get("REGPOLY_POOL_SIZE", "4"))
        analysis_default = max(1, pool_size // 2)
        return cls(
            db_url=os.environ.get(
                "REGPOLY_DB_URL",
                "postgresql://regpoly@localhost:5432/regpoly",
            ),
            host=os.environ.get("REGPOLY_HOST", "127.0.0.1"),
            port=int(os.environ.get("REGPOLY_PORT", "8000")),
            pool_size=pool_size,
            analysis_pool_size=int(os.environ.get(
                "REGPOLY_ANALYSIS_POOL_SIZE", str(analysis_default),
            )),
            progress_poll_seconds=float(
                os.environ.get("REGPOLY_POLL_SECONDS", "0.5")
            ),
            reload=os.environ.get("REGPOLY_RELOAD", "0") == "1",
            import_root=Path(root_env).resolve() if root_env else None,
            disable_internal_pool=(
                os.environ.get("REGPOLY_WEB_DISABLE_INTERNAL_POOL", "0") == "1"
            ),
            is_worker=os.environ.get("REGPOLY_WORKER", "0") == "1",
            worker_poll_seconds=float(
                os.environ.get("REGPOLY_WORKER_POLL_SECONDS", "0.5")
            ),
            db_pool_min_size=int(os.environ.get("REGPOLY_DB_POOL_MIN", "2")),
            db_pool_max_size=int(os.environ.get("REGPOLY_DB_POOL_MAX", "10")),
        )


WEB_PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = WEB_PACKAGE_DIR / "templates"
STATIC_DIR = WEB_PACKAGE_DIR / "static"
SCHEMA_PATH = WEB_PACKAGE_DIR / "schema.sql"  # legacy; not used on PG path


def find_docs_dir() -> Path | None:
    """Locate the docs/generators directory.

    Looks first at an env-var override, then relative to the package.
    """
    return _find_repo_subdir("REGPOLY_DOCS_DIR", "docs/generators")


def find_library_dir() -> Path | None:
    """Locate the docs/library directory holding published-generator YAMLs."""
    return _find_repo_subdir("REGPOLY_LIBRARY_DIR", "docs/library")


def find_papers_dir() -> Path | None:
    """Locate the top-level papers/ directory holding committed PDFs."""
    return _find_repo_subdir("REGPOLY_PAPERS_DIR", "papers")


def _find_repo_subdir(env_var: str, rel_path: str) -> Path | None:
    env = os.environ.get(env_var)
    if env:
        p = Path(env)
        if p.is_dir():
            return p
    here = WEB_PACKAGE_DIR
    for _ in range(6):
        candidate = here / rel_path
        if candidate.is_dir():
            return candidate
        here = here.parent
    return None
