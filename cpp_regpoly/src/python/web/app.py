"""FastAPI application for the regpoly web UI."""

from __future__ import annotations

import argparse
import asyncio
import sys
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from regpoly.web.config import (
    STATIC_DIR,
    TEMPLATES_DIR,
    Settings,
    find_library_dir,
    find_papers_dir,
)
from regpoly.library import Catalog
from regpoly.web.database import init_sync, open_async
from regpoly.web.routes import (
    families,
    generators,
    import_export,
    library,
    pages,
    primitive_search,
    tempering_search,
    tested_generators,
)
from regpoly.web.tasks.analysis import analyze_generator
from regpoly.web.tasks.pool import TaskPool


templates = Jinja2Templates(directory=str(TEMPLATES_DIR))

# Shared Jinja globals used by every template (sidebar, highlighting, etc.)
# Imported lazily to avoid a circular import with routes/families.py, which
# in turn imports from regpoly._regpoly_cpp.
from regpoly.web.routes.families import KNOWN_FAMILIES as _KNOWN_FAMILIES

# Sidebar split: the five flagship families come first in a fixed
# order; anything else the C++ registry exposes lands under an
# "Other families" heading, preserving KNOWN_FAMILIES order.
_PRIMARY_FAMILIES = [
    "TauswortheGen", "TGFSRGen", "MTGen", "WELLGen", "MELGGen",
]
_primary_set = set(_PRIMARY_FAMILIES)
templates.env.globals["families_primary"] = [
    {"name": f} for f in _PRIMARY_FAMILIES if f in _KNOWN_FAMILIES
]
templates.env.globals["families_other"] = [
    {"name": f} for f in _KNOWN_FAMILIES if f not in _primary_set
]
# Kept for templates that still iterate the flat list.
templates.env.globals["families"] = (
    templates.env.globals["families_primary"]
    + templates.env.globals["families_other"]
)
templates.env.globals.setdefault("active_family", None)

# Published-generators catalog — loaded once at import so the sidebar
# can surface papers without a per-request fetch.  The dev reload path
# lives in routes/library.py::_catalog.
_library_catalog = Catalog(find_library_dir())
_library_catalog.load()
templates.env.globals["library_papers"] = [
    {"id": p.id, "display": p.display(),
     "starred": p.starred, "year": p.year}
    for p in _library_catalog.papers()
]


def create_app(settings: Settings | None = None) -> FastAPI:
    if settings is None:
        settings = Settings.from_env()

    init_sync(settings.db_path)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        app.state.settings = settings
        app.state.db = await open_async(settings.db_path)
        app.state.templates = templates
        app.state.library = _library_catalog
        app.state.pool = TaskPool(
            db_path=settings.db_path, max_workers=settings.pool_size
        )
        analysis_task = asyncio.create_task(
            _analysis_scheduler(app, settings.db_path)
        )
        try:
            yield
        finally:
            analysis_task.cancel()
            # Cancel every still-running search so workers exit cleanly
            # on the next cancellation poll (instead of being orphaned).
            await app.state.db.execute(
                "UPDATE primitive_search_run SET status='cancelled' "
                "WHERE status IN ('pending', 'running')"
            )
            await app.state.db.execute(
                "UPDATE tempering_search_run SET status='cancelled' "
                "WHERE status IN ('pending', 'running')"
            )
            await app.state.db.commit()
            app.state.pool.shutdown()
            await app.state.db.close()

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

    return app


async def _analysis_scheduler(app: FastAPI, db_path: str) -> None:
    """Poll primitive_generator for rows needing analysis and dispatch
    them to the worker pool, one at a time."""
    import aiosqlite
    in_flight: set[int] = set()
    try:
        async with aiosqlite.connect(db_path) as conn:
            conn.row_factory = aiosqlite.Row
            while True:
                try:
                    async with conn.execute(
                        "SELECT id FROM primitive_generator "
                        "WHERE pis_computed_at IS NULL "
                        "ORDER BY id ASC LIMIT 50"
                    ) as cur:
                        rows = await cur.fetchall()
                    for r in rows:
                        gid = r["id"]
                        if gid in in_flight:
                            continue
                        in_flight.add(gid)
                        fut = app.state.pool.submit(
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
    parser.add_argument("--db", default=None, help="SQLite database path")
    parser.add_argument("--reload", action="store_true")
    args = parser.parse_args()

    settings = Settings.from_env()
    if args.host:
        settings.host = args.host
    if args.port:
        settings.port = args.port
    if args.db:
        settings.db_path = args.db
    if args.reload:
        settings.reload = True

    import uvicorn

    uvicorn.run(
        "regpoly.web.app:app",
        host=settings.host,
        port=settings.port,
        reload=settings.reload,
    )


if __name__ == "__main__":
    main()
