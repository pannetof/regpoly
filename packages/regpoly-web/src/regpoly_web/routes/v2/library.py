# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""GET /api/v2/library + /api/v2/library/papers + counts.

Backed by the in-memory Catalog held on app.state.library. Sliced and
diced for the chip+pagination contract on /library.
"""

from __future__ import annotations

from fastapi import APIRouter, Query, Request
from pydantic import BaseModel

router = APIRouter()


class _GenEntry(BaseModel):
    id: str
    paper_id: str
    family: str | None
    L: int | None
    k: int | None
    starred: bool


class GeneratorListEnvelope(BaseModel):
    rows: list[_GenEntry]
    total: int


class _PaperEntry(BaseModel):
    id: str
    display: str
    starred: bool
    year: int | None
    n_generators: int


class PaperListEnvelope(BaseModel):
    rows: list[_PaperEntry]
    total: int


def _gen_to_dict(paper, gen) -> _GenEntry:
    """Catalog generator entry → dict. Tolerant of attribute names
    across catalog versions; empty values stay None."""
    family = getattr(gen, "family", None)
    L = getattr(gen, "L", None)
    k = getattr(gen, "k", None)
    return _GenEntry(
        id=getattr(gen, "id", ""),
        paper_id=paper.id,
        family=family, L=L, k=k,
        starred=bool(getattr(paper, "starred", False)),
    )


@router.get("/library", response_model=GeneratorListEnvelope)
async def v2_list_library_generators(
    request: Request,
    family: str | None = None,
    paper_id: str | None = None,
    k: int | None = None,
    starred: bool | None = None,
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
) -> GeneratorListEnvelope:
    catalog = request.app.state.library
    rows: list[_GenEntry] = []
    for paper in catalog.papers():
        if paper_id and paper.id != paper_id:
            continue
        if starred is True and not paper.starred:
            continue
        if starred is False and paper.starred:
            continue
        for gen in getattr(paper, "generators", []) or []:
            row = _gen_to_dict(paper, gen)
            if family and row.family != family:
                continue
            if k is not None and row.k != k:
                continue
            rows.append(row)
    total = len(rows)
    page = rows[offset:offset + limit]
    return GeneratorListEnvelope(rows=page, total=total)


@router.get("/library/papers", response_model=PaperListEnvelope)
async def v2_list_library_papers(
    request: Request,
    starred: bool | None = None,
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
) -> PaperListEnvelope:
    catalog = request.app.state.library
    rows: list[_PaperEntry] = []
    for paper in catalog.papers():
        if starred is True and not paper.starred:
            continue
        if starred is False and paper.starred:
            continue
        rows.append(_PaperEntry(
            id=paper.id,
            display=paper.display(),
            starred=bool(getattr(paper, "starred", False)),
            year=getattr(paper, "year", None),
            n_generators=len(getattr(paper, "generators", []) or []),
        ))
    total = len(rows)
    page = rows[offset:offset + limit]
    return PaperListEnvelope(rows=page, total=total)


@router.get("/library/families/counts")
async def v2_library_families_counts(request: Request) -> dict[str, int]:
    """Per-family count across the published catalog (parity with
    /api/v2/generators/families/counts on the local-DB side)."""
    catalog = request.app.state.library
    counts: dict[str, int] = {}
    for paper in catalog.papers():
        for gen in getattr(paper, "generators", []) or []:
            family = getattr(gen, "family", None)
            if not family:
                continue
            counts[family] = counts.get(family, 0) + 1
    return counts
