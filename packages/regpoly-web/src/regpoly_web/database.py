# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""PostgreSQL database layer for the regpoly web application.

The dockerize-and-deploy migration (commit ``d22c384..``) switched the
backing store from SQLite to PostgreSQL 16. Both async (web request
handlers + worker scheduler) and sync (worker child processes) paths
are served by **psycopg3**:

- :func:`open_pool` returns an :class:`psycopg_pool.AsyncConnectionPool`
  used by the web app (mounted on ``app.state.pool``).
- :func:`sync_connect` opens a one-shot :class:`psycopg.Connection` for
  worker child processes (one connection per dispatched job).
- :func:`init_db` applies numbered migrations (``mNNN_*.py`` modules in
  this package) in version order. Idempotent. Schema only — never
  reaps.
- :func:`reap_orphans` cancels stale ``running`` rows whose
  ``updated_at`` is older than ``stale_seconds``. Never touches
  ``pending`` rows.

Engine notes:

- JSONB columns return native Python ``dict``/``list`` via psycopg3's
  default codec (no ``set_type_codec`` dance).
- Datetime columns return timezone-aware ``datetime`` (PG's
  ``TIMESTAMPTZ``).
- All queries use ``%s`` placeholders.
"""

from __future__ import annotations

import importlib
import json
import pkgutil
from contextlib import asynccontextmanager, contextmanager
from typing import Any, AsyncIterator, Iterator, Sequence

import psycopg
from psycopg.types.json import Jsonb
from psycopg_pool import AsyncConnectionPool


def _wrap_jsonb(params):
    """Auto-wrap dict/list params with :class:`Jsonb`.

    Routes that pull JSONB columns out of a SELECT and re-insert them
    (e.g. clone-on-restart) get native ``dict``/``list`` from psycopg's
    JSONB codec; psycopg's parameter side has no automatic encoder for
    those types. Wrapping here keeps the legacy callsites unchanged.
    """
    if params is None:
        return None
    if isinstance(params, dict):
        return params  # named-style dict; psycopg handles its members natively
    out = []
    for p in params:
        if isinstance(p, (dict, list)):
            out.append(Jsonb(p))
        else:
            out.append(p)
    return tuple(out)


# ─── Hybrid row (positional + named access, like aiosqlite.Row) ─────

class HybridRow(tuple):
    """A tuple that also allows ``row['column_name']`` lookup.

    aiosqlite + sqlite3.Row support both ``row[0]`` and ``row['name']``;
    psycopg's stock ``dict_row`` only exposes the dict view, and
    ``tuple_row`` only the tuple view. This factory keeps the
    pre-cutover route code working unchanged across both access
    patterns.
    """

    __slots__ = ()
    _fields: tuple[str, ...] = ()

    def __new__(cls, values):
        return super().__new__(cls, values)

    def __getitem__(self, key):  # type: ignore[override]
        if isinstance(key, str):
            # Case-insensitive lookup so legacy callers that read
            # `row["L"]` find the PG-lowercased `l` (PG folds
            # unquoted identifiers to lowercase). Keeps the cutover
            # diff small.
            try:
                idx = self._fields.index(key)
            except ValueError:
                lo = key.lower()
                for i, f in enumerate(self._fields):
                    if f.lower() == lo:
                        return super().__getitem__(i)
                raise KeyError(key)
            return super().__getitem__(idx)
        return super().__getitem__(key)

    def get(self, key, default=None):
        try:
            return self[key]
        except (KeyError, IndexError):
            return default

    def keys(self):
        return list(self._fields)


def _qmark_to_pyformat(sql: str) -> str:
    """Convert SQLite-style ``?`` placeholders to psycopg ``%s``.

    Skips ``?`` inside single- or double-quoted string literals so JSON
    fragments / regex patterns survive untouched.

    If the SQL already contains ``%s`` / ``%(name)s`` placeholders,
    treat it as already psycopg-formatted and return unchanged. This
    avoids double-escaping the literal ``%`` in pyformat tokens that
    rewritten code paths emit.
    """
    # Already pyformat? Don't touch.
    if "?" not in sql and "%s" in sql:
        return sql
    if "?" not in sql:
        return sql

    out: list[str] = []
    in_squote = False
    in_dquote = False
    i = 0
    while i < len(sql):
        c = sql[i]
        if c == "'" and not in_dquote:
            in_squote = not in_squote
            out.append(c)
        elif c == '"' and not in_squote:
            in_dquote = not in_dquote
            out.append(c)
        elif c == "?" and not in_squote and not in_dquote:
            out.append("%s")
        elif c == "%" and not in_squote and not in_dquote:
            # Escape so psycopg doesn't read it as a placeholder
            # marker.
            out.append("%%")
        else:
            out.append(c)
        i += 1
    return "".join(out)


def hybrid_row(cursor):
    """psycopg row factory producing :class:`HybridRow` instances."""
    fields = tuple(c.name for c in cursor.description) if cursor.description else ()
    cls = type("HybridRowFactory", (HybridRow,), {"_fields": fields})

    def make(values):
        return cls(values)

    return make

from regpoly_web import migrations as _migrations_pkg

# v2 (Phase 5.4): typed result tables (equidistribution_result,
# collision_free_result, tuplets_result) mirroring the C++ structs.
SCHEMA_VERSION = 2


# ─── connection helpers ─────────────────────────────────────────────

async def open_pool(
    db_url: str,
    *,
    min_size: int = 2,
    max_size: int = 10,
) -> AsyncConnectionPool:
    """Open an asyncpg-style pool of psycopg connections.

    Each connection uses ``dict_row`` so ``cursor.fetchone()`` returns
    a dict keyed by column name (matching the legacy
    :class:`aiosqlite.Row` ergonomics).
    """
    pool = AsyncConnectionPool(
        conninfo=db_url,
        min_size=min_size,
        max_size=max_size,
        kwargs={"row_factory": hybrid_row, "autocommit": True},
        open=False,
    )
    await pool.open()
    return pool


class AsyncPoolDB:
    """Thin shim presenting an aiosqlite-style API on top of an
    :class:`AsyncConnectionPool`.

    Routes pre-cutover called ``request.app.state.db.execute(...)``.
    Rewriting all ~54 sites to ``async with pool.connection() as conn:
    async with conn.cursor() as cur: await cur.execute(...)`` would
    have been a multi-day diff. This shim presents the legacy shape so
    only the SQL placeholders (``?`` → ``%s``) need updating, and
    JSONB columns coerce to ``dict``/``list`` natively.

    ``.execute(sql, params)`` returns an async context manager whose
    enter yields a cursor with ``fetchone`` / ``fetchall`` /
    ``rowcount``. The cursor's ``__aiter__`` is also supported (some
    routes use ``async for row in db.execute(...)``).

    ``.commit()`` is a no-op — the underlying pool is autocommit, so
    every statement is committed when it returns.
    """

    def __init__(self, pool: AsyncConnectionPool) -> None:
        self._pool = pool

    def execute(self, sql: str, params: Sequence | None = None):
        """Return an async context manager yielding a cursor."""
        return _ShimCursorCM(self._pool, _qmark_to_pyformat(sql),
                             _wrap_jsonb(params))

    async def executemany(self, sql: str, seq: Sequence[Sequence]) -> None:
        async with self._pool.connection() as conn:
            async with conn.cursor() as cur:
                await cur.executemany(
                    _qmark_to_pyformat(sql),
                    [_wrap_jsonb(p) for p in seq],
                )

    async def commit(self) -> None:  # noqa: D401 -- legacy compatibility
        """No-op; the pool is autocommit."""
        return

    async def close(self) -> None:
        await self._pool.close()


class _ShimCursorCM:
    """Hybrid wrapper for ``AsyncPoolDB.execute``.

    Two usage shapes are supported, matching the legacy aiosqlite API:

    1. ``async with db.execute(sql, params) as cur:`` — connection +
       cursor are released cleanly on exit.

    2. ``cur = await db.execute(sql, params)`` — fire-and-forget
       INSERT/UPDATE; the connection is held for the lifetime of the
       cursor object, then released by garbage collection. Suitable
       for write-only statements where the caller only reads
       ``.lastrowid`` / ``.rowcount`` or doesn't read at all.

    Two parallel uses in the same handler produce two distinct
    connections from the pool.
    """

    def __init__(self, pool: AsyncConnectionPool, sql: str, params):
        self._pool = pool
        self._sql = sql
        self._params = params
        self._conn_cm = None
        self._cur_cm = None
        self._cur = None

    async def _open(self):
        self._conn_cm = self._pool.connection()
        conn = await self._conn_cm.__aenter__()
        self._cur_cm = conn.cursor()
        self._cur = await self._cur_cm.__aenter__()
        sql = self._sql
        # Auto-append RETURNING id for bare INSERTs so .lastrowid works
        # the way aiosqlite consumers expect. Skip if the caller
        # already supplied RETURNING / ON CONFLICT DO NOTHING (no PK
        # is returned in that case).
        sql_upper = sql.lstrip().upper()
        self._lastrowid: int | None = None
        if (sql_upper.startswith("INSERT")
                and " RETURNING " not in sql_upper
                and "ON CONFLICT" not in sql_upper):
            sql = sql.rstrip().rstrip(";") + " RETURNING id"
            self._auto_returning = True
        else:
            self._auto_returning = False
        await self._cur.execute(sql, self._params)
        if self._auto_returning:
            try:
                row = await self._cur.fetchone()
                if row is not None:
                    self._lastrowid = int(row[0])
            except psycopg.ProgrammingError:
                # Statement returned no rows (e.g. INSERT into a
                # sequence-less table); leave lastrowid as None.
                self._lastrowid = None
        return self._cur

    # async-context-manager shape ──────────────────────────────────
    async def __aenter__(self):
        return await self._open()

    async def __aexit__(self, exc_type, exc, tb):
        try:
            if self._cur_cm is not None:
                await self._cur_cm.__aexit__(exc_type, exc, tb)
        finally:
            if self._conn_cm is not None:
                await self._conn_cm.__aexit__(exc_type, exc, tb)
            self._cur_cm = None
            self._conn_cm = None

    # direct-await shape (cur = await db.execute(...)) ─────────────
    def __await__(self):
        # Open the connection + cursor + run the query, returning a
        # `_ShimCursorProxy` that exposes the cursor's surface plus a
        # ``lastrowid`` property derived from PG's RETURNING clause
        # (when the SQL was an INSERT). The connection stays open
        # until the proxy is GC'd.
        return self._await_impl().__await__()

    async def _await_impl(self):
        await self._open()
        return _ShimCursorProxy(self)


class _ShimCursorProxy:
    """Object returned by ``await db.execute(...)``.

    Exposes the cursor's ``fetchone`` / ``fetchall`` / ``rowcount`` /
    ``lastrowid`` surface. ``lastrowid`` is computed by re-running the
    INSERT with a RETURNING clause when the original SQL doesn't
    include one — mirroring SQLite's auto-incremented PK shape.
    """

    def __init__(self, owner: "_ShimCursorCM"):
        self._owner = owner
        self._cur = owner._cur

    async def fetchone(self):
        return await self._cur.fetchone()

    async def fetchall(self):
        return await self._cur.fetchall()

    @property
    def rowcount(self) -> int:
        return self._cur.rowcount

    @property
    def lastrowid(self) -> int | None:
        return self._owner._lastrowid

    def __del__(self):
        # Best-effort release of the held pool connection. The pool's
        # connection close happens asynchronously; if a loop is still
        # running, schedule the release. Otherwise drop on the floor
        # (the pool's GC handles it).
        try:
            if self._owner._cur_cm is not None:
                import asyncio as _asyncio
                try:
                    loop = _asyncio.get_event_loop()
                except RuntimeError:
                    return
                if loop.is_running():
                    async def _close(o=self._owner):
                        try:
                            await o.__aexit__(None, None, None)
                        except Exception:
                            pass
                    loop.create_task(_close())
        except Exception:
            pass


class _SyncCursorProxy:
    """Wrapper exposing ``.lastrowid`` on a synchronous psycopg cursor.

    sqlite3 cursors expose ``.lastrowid`` after an ``INSERT``; psycopg
    does not. To keep the worker code (which calls ``cur.lastrowid``)
    unchanged, :meth:`SyncConnShim.execute` auto-appends ``RETURNING id``
    to bare ``INSERT`` statements and stashes the id here.
    """

    def __init__(self, cur, lastrowid: int | None) -> None:
        self._cur = cur
        self._lastrowid = lastrowid

    def __getattr__(self, name):
        return getattr(self._cur, name)

    def fetchone(self):
        return self._cur.fetchone()

    def fetchall(self):
        return self._cur.fetchall()

    @property
    def rowcount(self) -> int:
        return self._cur.rowcount

    @property
    def lastrowid(self) -> int | None:
        return self._lastrowid


class SyncConnShim:
    """Sync-connection wrapper used by worker child processes.

    Mirrors :class:`AsyncPoolDB`'s placeholder/Jsonb translation but
    over a single :class:`psycopg.Connection` (one per dispatched job).
    Also bridges legacy SQLite quirks the worker code still relies on:

    - ``?`` placeholders → ``%s`` (via :func:`_qmark_to_pyformat`).
    - dict/list params → :class:`Jsonb` wrapping (via :func:`_wrap_jsonb`).
    - ``INSERT OR IGNORE`` → ``INSERT ... ON CONFLICT DO NOTHING``.
    - ``conn.total_changes`` (sqlite3) → ``rowcount`` from the last
      executed cursor.
    - ``cur.lastrowid`` after a bare ``INSERT`` (auto-``RETURNING id``).
    """

    def __init__(self, conn: psycopg.Connection) -> None:
        self._conn = conn
        self._last_rowcount = 0

    def execute(self, sql: str, params: Sequence | None = None):
        sql = _qmark_to_pyformat(_translate_sqliteisms(sql))
        # Auto-append RETURNING id for bare INSERTs so cursor.lastrowid
        # works the way sqlite3 consumers expect.
        sql_upper = sql.lstrip().upper()
        wants_returning = (
            sql_upper.startswith("INSERT")
            and " RETURNING " not in sql_upper
            and "ON CONFLICT" not in sql_upper
        )
        run_sql = sql.rstrip().rstrip(";") + " RETURNING id" if wants_returning else sql
        cur = self._conn.execute(run_sql, _wrap_jsonb(params))
        lastrowid: int | None = None
        if wants_returning:
            try:
                row = cur.fetchone()
                if row is not None:
                    lastrowid = int(row[0])
            except psycopg.ProgrammingError:
                lastrowid = None
        self._last_rowcount = cur.rowcount
        return _SyncCursorProxy(cur, lastrowid)

    def cursor(self):
        return self._conn.cursor()

    @property
    def total_changes(self) -> int:
        return self._last_rowcount

    def commit(self) -> None:
        self._conn.commit()

    def rollback(self) -> None:
        self._conn.rollback()

    def close(self) -> None:
        self._conn.close()

    @property
    def autocommit(self) -> bool:
        return self._conn.autocommit

    @autocommit.setter
    def autocommit(self, value: bool) -> None:
        self._conn.autocommit = value


def _translate_sqliteisms(sql: str) -> str:
    """Rewrite a few SQLite-only fragments to PG equivalents.

    Kept narrow on purpose — every target is a verbatim string the
    worker code emits (no parsing). The static gate in CI flags any
    new SQLite syntax that creeps back in.
    """
    if "INSERT OR IGNORE" in sql:
        # Append ON CONFLICT DO NOTHING (idempotent if already present).
        if "ON CONFLICT" not in sql.upper():
            sql = sql.rstrip().rstrip(";") + " ON CONFLICT DO NOTHING"
        sql = sql.replace("INSERT OR IGNORE", "INSERT")
    return sql


@contextmanager
def sync_connect(db_url: str) -> Iterator[SyncConnShim]:
    """Open a one-shot synchronous psycopg connection.

    Used by worker child processes — one connection per dispatched job.
    The caller is responsible for ``conn.commit()`` between writes;
    autocommit is OFF so multi-statement progress updates can be
    grouped.

    Returns a :class:`SyncConnShim` rather than the raw psycopg
    connection so the worker code's legacy ``?`` placeholders and
    ``conn.total_changes`` keep working without per-call rewrites.
    """
    conn = psycopg.connect(db_url, row_factory=hybrid_row, autocommit=False)
    try:
        yield SyncConnShim(conn)
    finally:
        conn.close()


# ─── schema management ────────────────────────────────────────────────

def _discover_migrations() -> list[tuple[int, Any]]:
    """Return ``[(VERSION, module), ...]`` sorted ascending.

    Picks up every ``mNNN_*.py`` module in the migrations package.
    """
    out: list[tuple[int, Any]] = []
    for info in pkgutil.iter_modules(_migrations_pkg.__path__):
        if not info.name.startswith("m") or not info.name[1:4].isdigit():
            continue
        mod = importlib.import_module(f"regpoly_web.migrations.{info.name}")
        if hasattr(mod, "VERSION") and hasattr(mod, "apply"):
            out.append((int(mod.VERSION), mod))
    out.sort(key=lambda p: p[0])
    return out


async def init_db(db_url: str) -> None:
    """Apply any unapplied migrations to ``db_url``. Idempotent.

    Reads ``MAX(version)`` from ``schema_version``; runs each newer
    ``mNNN_*.py`` module's ``apply(conn)`` in version order; INSERTs
    each on success.
    """
    # Use a sync connection for migrations — psycopg's autocommit
    # semantics make DDL + setval simpler than wrestling with the
    # async pool's transaction boundaries.
    with sync_connect(db_url) as conn:
        conn.autocommit = True
        # Bootstrap: ensure schema_version exists before we read it
        # (m001 itself creates this table; chicken-and-egg avoided by
        # tolerating "table missing" on the first read).
        try:
            with conn.cursor() as cur:
                cur.execute("SELECT COALESCE(MAX(version), 0) FROM schema_version")
                row = cur.fetchone()
                current = (row.get("coalesce") if row else 0) or 0
        except psycopg.errors.UndefinedTable:
            conn.rollback()
            current = 0

        for version, module in _discover_migrations():
            if version <= current:
                continue
            module.apply(conn)
            with conn.cursor() as cur:
                cur.execute(
                    "INSERT INTO schema_version (version) VALUES (%s)",
                    (version,),
                )


async def reap_orphans(db_url: str, stale_seconds: int = 60) -> int:
    """Cancel ``running`` rows whose ``updated_at`` is older than
    ``stale_seconds``. Never touches ``pending`` rows.

    Returns the total number of rows updated across all *_run tables.
    """
    msg = f"orphaned (no heartbeat for >{stale_seconds}s)"
    total = 0
    with sync_connect(db_url) as conn:
        conn.autocommit = True
        for table in ("primitive_search_run", "tempering_search_run",
                      "library_test_run"):
            with conn.cursor() as cur:
                cur.execute(
                    f"""
                    UPDATE {table}
                       SET status='cancelled',
                           error_message=%s,
                           updated_at=NOW()
                     WHERE status='running'
                       AND updated_at < NOW() - INTERVAL '%s seconds'
                    """,
                    (msg, stale_seconds),
                )
                total += cur.rowcount
    return total


# ─── JSON helpers (kept for callers that need deterministic dumps) ───

def json_dumps(obj: Any) -> str:
    """Deterministic JSON for parameter storage (sorted keys)."""
    return json.dumps(obj, sort_keys=True, default=_json_default)


def _json_default(o: Any) -> Any:
    if hasattr(o, "__int__"):
        return int(o)
    raise TypeError(f"not JSON-serializable: {type(o).__name__}")


def json_loads(s: str | dict | list | None) -> Any:
    """Deserialize a JSON column value.

    Accepts both the legacy SQLite-TEXT shape (``str``) and the PG
    JSONB shape (``dict``/``list``) so callers don't need to branch.
    """
    if s is None or s == "":
        return None
    if isinstance(s, (dict, list)):
        return s
    return json.loads(s)
