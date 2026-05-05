"""HTML page routes (Jinja2-rendered).

P1 redesign: every render passes ``page_title`` and ``crumbs`` so the
new Tabler page-header partial renders the breadcrumb trail. The
``_render`` helper accepts those as keyword arguments alongside the
template-specific context. ``active_family`` is still threaded for
backwards-compatible templates that highlight a family.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import HTMLResponse, RedirectResponse

from regpoly_web.routes.families import KNOWN_FAMILIES

router = APIRouter()


def _render(
    request: Request,
    template: str,
    *,
    title: str = "regpoly",
    crumbs: list[dict] | None = None,
    pretitle: str | None = None,
    **context,
) -> HTMLResponse:
    """Render a Jinja2 template under the v2 Tabler shell.

    ``crumbs`` is a list of ``{"label": str, "href": str | None}`` dicts.
    The last entry is rendered as the active (non-link) crumb.

    When the request carries ``?family=<name>``, we thread that value
    through as ``active_family`` so the page-header pretitle can
    surface it. Explicit ``active_family`` in ``context`` wins over
    the query param.

    When the request carries ``?embed=1`` we hide the global chrome
    (sidebar, page header, footer) so the page can be shown inside an
    iframe.
    """
    templates = request.app.state.templates
    context.setdefault(
        "active_family", request.query_params.get("family") or None
    )
    context.setdefault(
        "hide_chrome", request.query_params.get("embed") == "1"
    )
    context["page_title"] = title
    context["page_pretitle"] = pretitle
    context["crumbs"] = crumbs if crumbs is not None else [
        {"label": "Dashboard", "href": "/"},
    ]
    response = templates.TemplateResponse(
        request=request, name=template, context=context
    )
    response.headers["Cache-Control"] = "no-cache, must-revalidate"
    return response


# ----- root + cross-section pages -------------------------------------


@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return _render(
        request, "index.html",
        title="Dashboard",
        crumbs=[{"label": "Dashboard", "href": None}],
        active_family=None,
    )


@router.get("/searches", response_class=HTMLResponse)
async def searches_page(request: Request):
    return _render(
        request, "searches.html",
        title="Searches",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Searches", "href": None},
        ],
        active_family=None,
    )


@router.get("/family/{family}", response_class=HTMLResponse)
async def family_landing_page(request: Request, family: str):
    if family not in KNOWN_FAMILIES:
        raise HTTPException(404, f"Unknown generator family: {family}")
    return _render(
        request, "family/landing.html",
        title=family,
        pretitle="Family",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Generators", "href": "/generators"},
            {"label": family, "href": None},
        ],
        family=family, active_family=family,
    )


# ----- generators -----------------------------------------------------


@router.get("/generators", response_class=HTMLResponse)
async def generators_page(request: Request):
    return _render(
        request, "generators/list.html",
        title="Generators",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Generators", "href": None},
        ],
    )


@router.get("/generators/compare", response_class=HTMLResponse)
async def generators_compare_page(request: Request):
    return _render(
        request, "compare/index.html",
        title="Compare generators",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Generators", "href": "/generators"},
            {"label": "Compare", "href": None},
        ],
    )


@router.get("/tools", response_class=HTMLResponse)
async def tools_page(request: Request):
    return _render(
        request, "tools/index.html",
        title="Tools",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Tools", "href": None},
        ],
    )


@router.get("/generators/{gen_id}", response_class=HTMLResponse)
async def generator_detail_page(request: Request, gen_id: int):
    return _render(
        request, "generators/detail.html",
        title=f"Generator #{gen_id}",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Generators", "href": "/generators"},
            {"label": f"#{gen_id}", "href": None},
        ],
        gen_id=gen_id,
    )


# ----- primitive search ----------------------------------------------


@router.get("/primitive-search", response_class=HTMLResponse)
async def primitive_search_page(request: Request):
    family = request.query_params.get("family")
    if not family:
        raise HTTPException(
            404, "Full-Period searches are launched from a family page."
        )
    return _render(
        request, "primitive_search/form.html",
        title="New full-period search",
        pretitle=family,
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": family, "href": f"/family/{family}"},
            {"label": "Full-period search", "href": None},
        ],
    )


@router.get("/primitive-search/{run_id}", response_class=HTMLResponse)
async def primitive_search_detail(request: Request, run_id: int):
    return _render(
        request, "primitive_search/detail.html",
        title=f"Full-period search #{run_id}",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Searches", "href": "/searches"},
            {"label": f"Full-period #{run_id}", "href": None},
        ],
        run_id=run_id,
    )


# ----- tempering search -----------------------------------------------


@router.get("/tempering-search", response_class=HTMLResponse)
async def tempering_search_page(request: Request):
    return _render(
        request, "tempering_search/form.html",
        title="New combined-generator search",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Searches", "href": "/searches"},
            {"label": "New combined search", "href": None},
        ],
    )


@router.get("/tempering-search/{run_id}", response_class=HTMLResponse)
async def tempering_search_detail(request: Request, run_id: int):
    return _render(
        request, "tempering_search/detail.html",
        title=f"Combined-generator search #{run_id}",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Searches", "href": "/searches"},
            {"label": f"Combined #{run_id}", "href": None},
        ],
        run_id=run_id,
    )


# ----- tested generators ----------------------------------------------


@router.get("/tested-generators", response_class=HTMLResponse)
async def tested_generators_page(request: Request):
    return _render(
        request, "tested_generators/list.html",
        title="Tested generators",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Tested generators", "href": None},
        ],
    )


@router.get("/tested-generators/{tg_id}", response_class=HTMLResponse)
async def tested_generator_detail(request: Request, tg_id: int):
    return _render(
        request, "tested_generators/detail.html",
        title=f"Tested generator #{tg_id}",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Tested generators", "href": "/tested-generators"},
            {"label": f"#{tg_id}", "href": None},
        ],
        tg_id=tg_id,
    )


# ----- library --------------------------------------------------------


@router.get("/library", response_class=HTMLResponse)
async def library_page(request: Request):
    return _render(
        request, "library/list.html",
        title="Library",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Library", "href": None},
        ],
    )


@router.get(
    "/library/papers",
    responses={
        307: {
            "description": "Redirect to /library?tab=papers — the "
                           "library merge in the v2 redesign collapses "
                           "the standalone papers page into a tab.",
        },
    },
)
async def library_papers_page(request: Request):
    """v2 redesign — /library/papers preserved as a 307 redirect to
    /library?tab=papers so S7.2 user-story URLs keep working without
    a flag-day for any client that bookmarked them."""
    return RedirectResponse(url="/library?tab=papers", status_code=307)


@router.get("/library/{paper_id}", response_class=HTMLResponse)
async def library_paper_page(request: Request, paper_id: str):
    cat = request.app.state.library
    if cat.paper(paper_id) is None:
        raise HTTPException(404, f"Unknown paper id: {paper_id}")
    return _render(
        request, "library/paper.html",
        title=paper_id,
        pretitle="Paper",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Library", "href": "/library"},
            {"label": paper_id, "href": None},
        ],
        paper_id=paper_id,
    )


@router.get("/library/{paper_id}/{gen_id}", response_class=HTMLResponse)
async def library_generator_page(
    request: Request, paper_id: str, gen_id: str,
):
    cat = request.app.state.library
    loc = cat.generator(gen_id)
    if loc is None:
        raise HTTPException(404, f"Unknown generator id: {gen_id}")
    if loc[0].id != paper_id:
        raise HTTPException(
            404, f"Generator {gen_id!r} belongs to {loc[0].id!r}")
    return _render(
        request, "library/detail.html",
        title=gen_id,
        pretitle=f"Library entry from {paper_id}",
        crumbs=[
            {"label": "Dashboard", "href": "/"},
            {"label": "Library", "href": "/library"},
            {"label": paper_id, "href": f"/library/{paper_id}"},
            {"label": gen_id, "href": None},
        ],
        paper_id=paper_id, gen_id=gen_id,
    )
