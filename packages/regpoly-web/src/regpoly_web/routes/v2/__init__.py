"""regpoly-web /api/v2/ namespace.

The v2 namespace is where every breaking-shape change in the redesign
lives. /api/ (v1) is preserved verbatim. Each phase's new endpoints
live here.

Mounted by app.py via:
    app.include_router(v2.router, prefix="/api/v2", tags=["v2"])
"""

from __future__ import annotations

from fastapi import APIRouter

from regpoly_web.routes.v2 import (
    dashboard,
    generators,
    library,
    publish,
    related,
    searches,
    searches_history,
    tested_generators,
    transition_matrix,
)

router = APIRouter()

# Mount every sub-router under /api/v2/.
router.include_router(dashboard.router)
router.include_router(generators.router)
router.include_router(tested_generators.router)
router.include_router(library.router)
router.include_router(searches.router)
router.include_router(searches_history.router)
router.include_router(publish.router)
router.include_router(related.router)
router.include_router(transition_matrix.router)


@router.get("/healthz")
def healthz() -> dict:
    return {"status": "ok", "version": 2}
