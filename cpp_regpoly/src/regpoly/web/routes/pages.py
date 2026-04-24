"""HTML page routes (Jinja2-rendered)."""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import HTMLResponse

from regpoly.web.routes.families import KNOWN_FAMILIES


router = APIRouter()


def _render(request: Request, template: str, **context) -> HTMLResponse:
    """Render a Jinja2 template.

    When the request carries ``?family=<name>``, we thread that value
    through as ``active_family`` so the sidebar can highlight it.
    Explicit ``active_family`` in ``context`` wins over the query param.

    When the request carries ``?embed=1`` we hide the global chrome
    (sidebar, top-bar, footer) so the page can be shown inside an
    iframe without a nested navigation column.
    """
    templates = request.app.state.templates
    context.setdefault(
        "active_family", request.query_params.get("family") or None
    )
    context.setdefault(
        "hide_chrome", request.query_params.get("embed") == "1"
    )
    return templates.TemplateResponse(
        request=request, name=template, context=context
    )


@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return _render(request, "index.html", active_family=None)


@router.get("/searches", response_class=HTMLResponse)
async def searches_page(request: Request):
    return _render(request, "searches.html", active_family=None)


@router.get("/family/{family}", response_class=HTMLResponse)
async def family_landing_page(request: Request, family: str):
    if family not in KNOWN_FAMILIES:
        raise HTTPException(404, f"Unknown generator family: {family}")
    return _render(
        request, "family/landing.html",
        family=family, active_family=family,
    )


@router.get("/generators", response_class=HTMLResponse)
async def generators_page(request: Request):
    return _render(request, "generators/list.html")


@router.get("/generators/{gen_id}", response_class=HTMLResponse)
async def generator_detail_page(request: Request, gen_id: int):
    return _render(request, "generators/detail.html", gen_id=gen_id)


@router.get("/primitive-search", response_class=HTMLResponse)
async def primitive_search_page(request: Request):
    # The unscoped full-period search page no longer exists; callers
    # must arrive via a family landing page (which supplies ?family=X).
    if not request.query_params.get("family"):
        raise HTTPException(
            404, "Full-Period searches are launched from a family page."
        )
    return _render(request, "primitive_search/form.html")


@router.get("/primitive-search/{run_id}", response_class=HTMLResponse)
async def primitive_search_detail(request: Request, run_id: int):
    return _render(
        request, "primitive_search/detail.html", run_id=run_id
    )


@router.get("/tempering-search", response_class=HTMLResponse)
async def tempering_search_page(request: Request):
    return _render(request, "tempering_search/form.html")


@router.get("/tempering-search/{run_id}", response_class=HTMLResponse)
async def tempering_search_detail(request: Request, run_id: int):
    return _render(
        request, "tempering_search/detail.html", run_id=run_id
    )


@router.get("/tested-generators", response_class=HTMLResponse)
async def tested_generators_page(request: Request):
    return _render(request, "tested_generators/list.html")


@router.get("/tested-generators/{tg_id}", response_class=HTMLResponse)
async def tested_generator_detail(request: Request, tg_id: int):
    return _render(
        request, "tested_generators/detail.html", tg_id=tg_id
    )


@router.get("/library", response_class=HTMLResponse)
async def library_page(request: Request):
    return _render(request, "library/list.html")


@router.get("/library/papers", response_class=HTMLResponse)
async def library_papers_page(request: Request):
    return _render(request, "library/papers.html")


@router.get("/library/{paper_id}", response_class=HTMLResponse)
async def library_paper_page(request: Request, paper_id: str):
    cat = request.app.state.library
    if cat.paper(paper_id) is None:
        raise HTTPException(404, f"Unknown paper id: {paper_id}")
    return _render(request, "library/paper.html", paper_id=paper_id)


@router.get("/library/{paper_id}/{gen_id}", response_class=HTMLResponse)
async def library_generator_page(
    request: Request, paper_id: str, gen_id: str,
):
    cat = request.app.state.library
    loc = cat.generator(gen_id)
    if loc is None:
        raise HTTPException(404, f"Unknown generator id: {gen_id}")
    if loc[0].id != paper_id:
        # Redirect callers who used the wrong paper slug.
        raise HTTPException(
            404, f"Generator {gen_id!r} belongs to {loc[0].id!r}")
    return _render(
        request, "library/detail.html",
        paper_id=paper_id, gen_id=gen_id)
