# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""FastAPI application for the regpoly web UI."""

from __future__ import annotations

import argparse
import asyncio
import logging
import sys
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from regpoly.library import Catalog

from regpoly_web.config import (
    STATIC_DIR,
    TEMPLATES_DIR,
    Settings,
    find_library_dir,
    find_papers_dir,
)
from regpoly_web.database import AsyncPoolDB, open_pool, reap_orphans
from regpoly_web.routes import (
    families,
    generators,
    health,
    import_export,
    library,
    pages,
    primitive_search,
    tempering_search,
    tested_generators,
)
from regpoly_web.routes.v2 import router as v2_router
from regpoly_web.tasks.analysis import analyze_generator
from regpoly_web.tasks.pool import TaskPool

logger = logging.getLogger(__name__)

templates = Jinja2Templates(directory=str(TEMPLATES_DIR))

# Shared Jinja globals used by every template (sidebar, highlighting,
# breadcrumbs, etc.). The v2 redesign replaces the family-hierarchy
# sidebar with a flat 6-entry nav; family browsing is driven by chips
# on /generators (P2) and the family directory on / (P2).
# Imported lazily to avoid a circular import with routes/families.py.
from regpoly_web.routes.families import KNOWN_FAMILIES as _KNOWN_FAMILIES  # noqa: F401

templates.env.globals["nav_items"] = [
    {"label": "Dashboard",          "href": "/",                   "icon": "home"},
    {"label": "Generators",         "href": "/generators",         "icon": "cpu"},
    {"label": "Combined generators","href": "/tested-generators",  "icon": "test-pipe"},
    {"label": "Searches",           "href": "/searches",           "icon": "search"},
    {"label": "Library",            "href": "/library",            "icon": "book"},
    # Tools surface (Import / Export / Imports audit trail) lands in P5.
    # The entry exists from P1 so the sidebar IA is stable from day one.
    {"label": "Tools",              "href": "/tools",              "icon": "tools"},
]
templates.env.globals.setdefault("active_family", None)
# Empty list kept (legacy templates may still iterate it without crashing).
templates.env.globals.setdefault("crumbs", [])

# Published-generators catalog — loaded once at import. The redesigned
# sidebar no longer surfaces papers; the Library page (P2) consumes the
# catalog directly from app.state.library.
_library_catalog = Catalog(find_library_dir())
_library_catalog.load()


def create_app(settings: Settings | None = None) -> FastAPI:
    if settings is None:
        settings = Settings.from_env()

    # NOTE: `init_db` is async and runs in the `init` compose container
    # (or once before web serves traffic). We do not call it from
    # module-import-time; doing so would break the `app = create_app()`
    # pattern at line 206 (no event loop yet) AND would mean every
    # restart of the web container reapplied migrations races with the
    # init container. The schema authority lives in the init container.

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        app.state.settings = settings
        app.state.templates = templates
        app.state.library = _library_catalog
        # asyncpg-style pool of psycopg connections — exposed in two
        # shapes:
        #   - app.state.dbpool is the raw psycopg_pool.AsyncConnectionPool
        #     (used by /healthz and the new worker-side schedulers).
        #   - app.state.db is an aiosqlite-shaped shim on top of the
        #     pool so the ~54 pre-cutover route sites keep their
        #     `async with db.execute(...)` shape; only the SQL
        #     placeholders (`?` → `%s`) and JSON column casts need
        #     updating.
        # NOTE: `app.state.pool` continues to mean the legacy TaskPool
        # (ProcessPoolExecutor wrapper) so existing routes' .submit()
        # calls keep working in dev mode without modification.
        app.state.dbpool = await open_pool(
            settings.db_url,
            min_size=settings.db_pool_min_size,
            max_size=settings.db_pool_max_size,
        )
        app.state.db = AsyncPoolDB(app.state.dbpool)

        # Production (compose-driven) mode: web container owns NO
        # ProcessPoolExecutor and runs NO analysis scheduler. The
        # worker container handles all dispatch. Routes that legacy-
        # called `app.state.pool.submit(...)` should null-check first.
        analysis_task = None
        if settings.disable_internal_pool:
            app.state.pool = None
            app.state.analysis_pool = None
        else:
            # Single-process dev mode: spin up the in-process pools so
            # `uv run regpoly-web` against a dev PG keeps working
            # without needing a separate worker process.
            app.state.pool = TaskPool(
                db_url=settings.db_url, max_workers=settings.pool_size,
            )
            app.state.analysis_pool = TaskPool(
                db_url=settings.db_url, max_workers=settings.analysis_pool_size,
            )
            # Reap any orphaned 'running' rows from a previous dev
            # session (the init container does this in production).
            try:
                n = await reap_orphans(settings.db_url, stale_seconds=0)
                if n:
                    logger.info("dev-mode startup reaped %d orphan rows", n)
            except Exception:
                logger.warning("dev-mode reap_orphans failed", exc_info=True)
            analysis_task = asyncio.create_task(
                _analysis_scheduler(app, settings.db_url)
            )
        try:
            yield
        finally:
            if analysis_task is not None:
                analysis_task.cancel()
            if not settings.disable_internal_pool:
                # Cancel running rows so dev-mode workers exit cleanly.
                async with app.state.dbpool.connection() as conn:
                    async with conn.cursor() as cur:
                        await cur.execute(
                            "UPDATE primitive_search_run SET status='cancelled' "
                            "WHERE status IN ('pending', 'running')"
                        )
                        await cur.execute(
                            "UPDATE tempering_search_run SET status='cancelled' "
                            "WHERE status IN ('pending', 'running')"
                        )
                if app.state.pool is not None:
                    app.state.pool.shutdown()
                if app.state.analysis_pool is not None:
                    app.state.analysis_pool.shutdown()
            await app.state.dbpool.close()

    app = FastAPI(
        title="regpoly web",
        description="Search and analysis for GF(2)-linear PRNGs",
        lifespan=lifespan,
    )

    class NoCacheStaticFiles(StaticFiles):
        async def get_response(self, path, scope):
            resp = await super().get_response(path, scope)
            # Dev UI files change constantly; force the browser to
            # revalidate so edits don't get masked by a stale cache.
            resp.headers["Cache-Control"] = "no-cache, must-revalidate"
            return resp

    app.mount(
        "/static", NoCacheStaticFiles(directory=str(STATIC_DIR)), name="static"
    )

    # Serve committed paper PDFs (see docs/library/*.yaml#reference.pdf).
    _papers_dir = find_papers_dir()
    if _papers_dir is not None:
        app.mount(
            "/papers",
            StaticFiles(directory=str(_papers_dir)),
            name="papers",
        )

    app.include_router(pages.router)
    app.include_router(families.router, prefix="/api")
    app.include_router(generators.router, prefix="/api")
    app.include_router(library.router, prefix="/api")
    app.include_router(primitive_search.router, prefix="/api")
    app.include_router(tempering_search.router, prefix="/api")
    app.include_router(tested_generators.router, prefix="/api")
    app.include_router(import_export.router, prefix="/api")
    # Unauthenticated probe; intentionally not under /api so the Caddy
    # basic-auth exemption is one exact path, not a glob.
    app.include_router(health.router)

    # v2 namespace — every contract change in the redesign mounts
    # under /api/v2/. v1 (existing /api/) remains unchanged.
    app.include_router(v2_router, prefix="/api/v2", tags=["v2"])
    app.state.v2_router_registered = True

    return app


async def _analysis_scheduler(app: FastAPI, db_url: str) -> None:
    """Poll primitive_generator for rows needing analysis and dispatch
    them to the worker pool, one at a time.

    Used only in single-process dev mode; in production the worker
    container's own scheduler (regpoly_web.worker) handles dispatch.
    """
    in_flight: set[int] = set()
    try:
        while True:
            try:
                async with app.state.dbpool.connection() as conn:
                    async with conn.cursor() as cur:
                        await cur.execute(
                            "SELECT id FROM primitive_generator "
                            "WHERE pis_computed_at IS NULL "
                            "ORDER BY id ASC LIMIT 50"
                        )
                        rows = await cur.fetchall()
                for r in rows:
                    gid = r["id"]
                    if gid in in_flight:
                        continue
                    in_flight.add(gid)
                    fut = app.state.analysis_pool.submit(
                        "analysis", gid, analyze_generator
                    )
                    fut.add_done_callback(
                        lambda _f, g=gid: in_flight.discard(g)
                    )
            except Exception:
                # Never let scheduler errors kill the loop.
                pass
            await asyncio.sleep(2.0)
    except asyncio.CancelledError:
        pass


# Exported for uvicorn/gunicorn factory-style loading
app = create_app()


def main() -> None:
    parser = argparse.ArgumentParser(description="regpoly web UI")
    parser.add_argument("--host", default=None)
    parser.add_argument("--port", type=int, default=None)
    parser.add_argument(
        "--db-url", default=None,
        help="PostgreSQL DSN (default: $REGPOLY_DB_URL)",
    )
    parser.add_argument("--reload", action="store_true")
    args = parser.parse_args()

    settings = Settings.from_env()
    if args.host:
        settings.host = args.host
    if args.port:
        settings.port = args.port
    if args.db_url:
        settings.db_url = args.db_url
    if args.reload:
        settings.reload = True

    import uvicorn

    uvicorn.run(
        "regpoly_web.app:app",
        host=settings.host,
        port=settings.port,
        reload=settings.reload,
    )


if __name__ == "__main__":
    main()
