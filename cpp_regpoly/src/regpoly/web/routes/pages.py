"""HTML page routes (Jinja2-rendered)."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse


router = APIRouter()


def _render(request: Request, template: str, **context) -> HTMLResponse:
    templates = request.app.state.templates
    return templates.TemplateResponse(
        request=request, name=template, context=context
    )


@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return _render(request, "index.html")


@router.get("/searches", response_class=HTMLResponse)
async def searches_page(request: Request):
    return _render(request, "searches.html")


@router.get("/generators", response_class=HTMLResponse)
async def generators_page(request: Request):
    return _render(request, "generators/list.html")


@router.get("/generators/{gen_id}", response_class=HTMLResponse)
async def generator_detail_page(request: Request, gen_id: int):
    return _render(request, "generators/detail.html", gen_id=gen_id)


@router.get("/primitive-search", response_class=HTMLResponse)
async def primitive_search_page(request: Request):
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
