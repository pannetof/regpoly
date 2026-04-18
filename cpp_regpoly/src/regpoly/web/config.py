"""Runtime configuration for the regpoly web application."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Settings:
    """Web application settings, overridable via environment variables."""

    db_path: str = "regpoly.db"
    host: str = "127.0.0.1"
    port: int = 8000
    pool_size: int = 4
    progress_poll_seconds: float = 0.5
    reload: bool = False

    @classmethod
    def from_env(cls) -> "Settings":
        return cls(
            db_path=os.environ.get("REGPOLY_DB", "regpoly.db"),
            host=os.environ.get("REGPOLY_HOST", "127.0.0.1"),
            port=int(os.environ.get("REGPOLY_PORT", "8000")),
            pool_size=int(os.environ.get("REGPOLY_POOL_SIZE", "4")),
            progress_poll_seconds=float(
                os.environ.get("REGPOLY_POLL_SECONDS", "0.5")
            ),
            reload=os.environ.get("REGPOLY_RELOAD", "0") == "1",
        )


WEB_PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = WEB_PACKAGE_DIR / "templates"
STATIC_DIR = WEB_PACKAGE_DIR / "static"
SCHEMA_PATH = WEB_PACKAGE_DIR / "schema.sql"


def find_docs_dir() -> Path | None:
    """Locate the docs/generators directory.

    Looks first at an env-var override, then relative to the package.
    """
    env = os.environ.get("REGPOLY_DOCS_DIR")
    if env:
        p = Path(env)
        if p.is_dir():
            return p
    # The docs live at <repo>/docs/generators in source checkouts.
    # Walk up from the package looking for a sibling 'docs' dir.
    here = WEB_PACKAGE_DIR
    for _ in range(6):
        candidate = here / "docs" / "generators"
        if candidate.is_dir():
            return candidate
        here = here.parent
    return None
