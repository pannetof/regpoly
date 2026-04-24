"""SQLite database layer for the regpoly web application.

Provides both synchronous (for background processes) and asynchronous
(for FastAPI request handlers) access to the same SQLite file.  WAL mode
is enabled so concurrent readers and one writer can coexist.
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

import aiosqlite

from regpoly.web.config import SCHEMA_PATH


SCHEMA_VERSION = 1


def _apply_pragmas(conn: sqlite3.Connection | aiosqlite.Connection) -> None:
    """Enable WAL and a sane busy_timeout for concurrent access."""
    # Works for both sqlite3 and aiosqlite (via execute coroutine wrapping).
    pass  # handled per-connection in init_sync / get_async


def _schema_sql() -> str:
    return Path(SCHEMA_PATH).read_text()


def _migrate(conn: sqlite3.Connection) -> None:
    """Add columns introduced after v1 to existing databases, then
    create indexes that reference them."""
    expected = {
        "primitive_generator": [
            ("char_poly", "TEXT"),
            ("hamming_weight", "INTEGER"),
            ("pis_gaps", "TEXT"),
            ("pis_se", "INTEGER"),
            ("pis_elapsed", "REAL"),
            ("pis_computed_at", "TEXT"),
            ("pis_error", "TEXT"),
            ("library_id", "TEXT"),
        ],
        "primitive_search_run": [
            ("search_mode", "TEXT NOT NULL DEFAULT 'random'"),
            ("enum_index",  "INTEGER NOT NULL DEFAULT 0"),
            ("enum_total",  "TEXT"),
            ("enum_axes",   "TEXT"),
        ],
        "tested_generator": [
            ("library_id", "TEXT"),
        ],
    }
    for table, columns in expected.items():
        existing = {row[1] for row in conn.execute(
            f"PRAGMA table_info({table})"
        ).fetchall()}
        for name, col_type in columns:
            if name not in existing:
                conn.execute(
                    f"ALTER TABLE {table} ADD COLUMN {name} {col_type}"
                )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pg_pending_analysis "
        "ON primitive_generator(pis_computed_at, pis_error)"
    )
    # Indexes for library_id lookup on older DBs that lacked the column.
    conn.execute(
        "CREATE UNIQUE INDEX IF NOT EXISTS idx_tg_library_id "
        "ON tested_generator(library_id) WHERE library_id IS NOT NULL"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pg_library_id "
        "ON primitive_generator(library_id) WHERE library_id IS NOT NULL"
    )
    # One-shot rename of the TGFSR family: older DBs store the C++
    # class name TGFSRGen; the canonical name is now TGFSR.  The
    # factory still accepts the legacy name, but the sidebar + API
    # filters key off the canonical one, so normalise stored rows.
    conn.execute(
        "UPDATE primitive_generator SET family = 'TGFSR' "
        "WHERE family = 'TGFSRGen'"
    )


def init_sync(db_path: str) -> None:
    """Create the database file (if needed) and apply the schema.

    Also reaps any searches left in 'running' or 'pending' state — they
    are necessarily orphans, since no worker from a prior server
    invocation is still alive when we reach this function.

    Called once at application startup (synchronous).
    """
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("PRAGMA busy_timeout=5000")
        conn.executescript(_schema_sql())
        _migrate(conn)
        row = conn.execute(
            "SELECT MAX(version) FROM schema_version"
        ).fetchone()
        if row is None or row[0] is None:
            conn.execute(
                "INSERT INTO schema_version(version) VALUES (?)",
                (SCHEMA_VERSION,),
            )
        reap_msg = "orphaned by previous server shutdown"
        conn.execute(
            "UPDATE primitive_search_run "
            "SET status='cancelled', error_message=? "
            "WHERE status IN ('running', 'pending')",
            (reap_msg,),
        )
        conn.execute(
            "UPDATE tempering_search_run "
            "SET status='cancelled', error_message=? "
            "WHERE status IN ('running', 'pending')",
            (reap_msg,),
        )
        conn.commit()
    finally:
        conn.close()


@contextmanager
def sync_connect(db_path: str) -> Iterator[sqlite3.Connection]:
    """Open a synchronous SQLite connection with PRAGMAs applied.

    Meant for use inside worker processes.
    """
    conn = sqlite3.connect(db_path, timeout=30.0)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA busy_timeout=5000")
    try:
        yield conn
    finally:
        conn.close()


async def open_async(db_path: str) -> aiosqlite.Connection:
    """Open an async connection with PRAGMAs and Row factory applied."""
    conn = await aiosqlite.connect(db_path)
    conn.row_factory = aiosqlite.Row
    await conn.execute("PRAGMA foreign_keys=ON")
    await conn.execute("PRAGMA busy_timeout=5000")
    return conn


def json_dumps(obj: Any) -> str:
    """Deterministic JSON for parameter storage (sorted keys)."""
    return json.dumps(obj, sort_keys=True, default=_json_default)


def _json_default(o: Any) -> Any:
    # Fallbacks for numpy-style ints, etc.
    if hasattr(o, "__int__"):
        return int(o)
    raise TypeError(f"not JSON-serializable: {type(o).__name__}")


def json_loads(s: str | None) -> Any:
    if s is None or s == "":
        return None
    return json.loads(s)
