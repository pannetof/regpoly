"""API for the published-generators library (paper-centric)."""

from __future__ import annotations

import asyncio

from fastapi import APIRouter, HTTPException, Request

from regpoly.library import Catalog, Generator, Paper
from regpoly.web.database import json_dumps, sync_connect
from regpoly.web.param_format import (
    format_gen_params, format_tempering_list,
)
from regpoly.web.routes.families import _markdown_to_html


router = APIRouter()


def _catalog(request: Request) -> Catalog:
    cat: Catalog = request.app.state.library
    if getattr(request.app.state, "settings", None) and \
            request.app.state.settings.reload:
        cat.reload_if_stale()
    return cat


async def _instantiated_ids(db, library_id: str) -> dict:
    out = {"primitive_ids": [], "tested_ids": []}
    async with db.execute(
        "SELECT id FROM primitive_generator WHERE library_id = ? ORDER BY id",
        (library_id,),
    ) as cur:
        out["primitive_ids"] = [r[0] for r in await cur.fetchall()]
    async with db.execute(
        "SELECT id FROM tested_generator WHERE library_id = ? ORDER BY id",
        (library_id,),
    ) as cur:
        out["tested_ids"] = [r[0] for r in await cur.fetchall()]
    return out


# ─── Paper-level views ────────────────────────────────────────────────────

def _paper_summary(p: Paper) -> dict:
    return {
        "id": p.id,
        "display": p.display(),
        "authors": [
            {"family": a.family, "given": a.given} for a in p.authors
        ],
        "authors_short": p.author_list_short(),
        "year": p.year,
        "title": p.title,
        "venue": p.venue,
        "volume": p.volume,
        "issue": p.issue,
        "pages": p.pages,
        "doi": p.doi,
        "pdf": p.pdf,
        "bibkey": p.bibkey,
        "citation": p.acmtrans_citation(),
        "starred": p.starred,
        "tags": p.tags,
        "num_generators": len(p.generators),
        "valid": p.valid,
        "errors": p.errors,
    }


def _gen_summary(g: Generator) -> dict:
    # Per-component tempering chain, with bitmask params rendered as
    # hex strings via the shared param_format helper.  Structure:
    #   [[{"type":"tempMK2","params":{"eta":7,"mu":15,"b":"0x9D2C5680",…}}], …]
    tempering_summary: list[list[dict]] = []
    for c in g.components:
        formatted_chain = format_tempering_list(c.get("tempering") or [])
        chain: list[dict] = []
        for step in formatted_chain:
            params = {k: v for k, v in step.items() if k != "type"}
            chain.append({"type": step.get("type", ""),
                          "params": params})
        tempering_summary.append(chain)

    return {
        "id": g.id,
        "display": g.display,
        "family": g.family,
        "target": g.target,
        "combined": g.combined,
        "Lmax": g.Lmax,
        "J": len(g.components),
        "tempering_summary": tempering_summary,
        "starred": g.starred,
        "valid": g.valid,
        "errors": g.errors,
    }


@router.get("/library/papers")
async def list_papers(
    request: Request,
    starred: bool | None = None,
    tag: str | None = None,
    include_invalid: bool = False,
) -> list[dict]:
    cat = _catalog(request)
    out = []
    for p in cat.papers(starred=starred, tag=tag,
                        include_invalid=include_invalid):
        row = _paper_summary(p)
        row["generators"] = [_gen_summary(g) for g in p.generators]
        for g_row, g in zip(row["generators"], p.generators):
            g_row["instantiated"] = await _instantiated_ids(
                request.app.state.db, g.id)
        out.append(row)
    return out


@router.get("/library/papers/{paper_id}")
async def get_paper(request: Request, paper_id: str) -> dict:
    cat = _catalog(request)
    p = cat.paper(paper_id)
    if p is None:
        raise HTTPException(404, f"Unknown paper id: {paper_id}")
    row = _paper_summary(p)
    row["abstract_md"] = p.abstract_md
    row["abstract_html"] = (
        _markdown_to_html(p.abstract_md) if p.abstract_md else "")
    row["notes_md"] = p.notes_md
    row["notes_html"] = (
        _markdown_to_html(p.notes_md) if p.notes_md else "")
    gens = []
    for g in p.generators:
        gr = _gen_summary(g)
        gr["instantiated"] = await _instantiated_ids(
            request.app.state.db, g.id)
        gens.append(gr)
    row["generators"] = gens
    return row


# ─── Generator-level views ────────────────────────────────────────────────

@router.get("/library/generators")
async def list_generators(
    request: Request,
    starred: bool | None = None,
    family: str | None = None,
) -> list[dict]:
    """Flat list of library generators across every paper.

    ``starred=true`` restricts to generators explicitly flagged in
    their YAML entry — used by the ``/library`` landing page.  Each
    row carries minimal paper context so the UI can link back.
    """
    cat = _catalog(request)
    out = []
    for paper, g in cat.all_generators(family=family):
        if starred is not None and g.starred != starred:
            continue
        row = _gen_summary(g)
        row["paper"] = {
            "id": paper.id,
            "display": paper.display(),
            "year": paper.year,
            "authors_short": paper.author_list_short(),
            "title": paper.title,
            "doi": paper.doi,
        }
        row["instantiated"] = await _instantiated_ids(
            request.app.state.db, g.id)
        out.append(row)
    return out


@router.get("/library/generators/{gen_id}")
async def get_generator(request: Request, gen_id: str) -> dict:
    cat = _catalog(request)
    loc = cat.generator(gen_id)
    if loc is None:
        raise HTTPException(404, f"Unknown generator id: {gen_id}")
    paper, g = loc
    row = _gen_summary(g)
    row["components"] = [
        {
            **c,
            "params": format_gen_params(c.get("family", ""),
                                        c.get("params") or {}),
            "tempering": format_tempering_list(c.get("tempering") or []),
        }
        for c in g.components
    ]
    row["notes_md"] = g.notes_md
    row["notes_html"] = (
        _markdown_to_html(g.notes_md) if g.notes_md else "")
    row["instantiated"] = await _instantiated_ids(
        request.app.state.db, g.id)
    row["paper"] = _paper_summary(paper)
    return row


# ─── Write endpoints ─────────────────────────────────────────────────────

@router.post("/library/generators/{gen_id}/instantiate")
async def instantiate(request: Request, gen_id: str) -> dict:
    cat = _catalog(request)
    loc = cat.generator(gen_id)
    if loc is None:
        raise HTTPException(404, f"Unknown generator id: {gen_id}")
    _, g = loc
    if not g.valid:
        raise HTTPException(
            400, {"code": "invalid_entry", "errors": g.errors})
    db_path = request.app.state.settings.db_path
    tested_id, created = await asyncio.to_thread(
        _instantiate_sync, db_path, g)
    return {"library_id": g.id,
            "tested_generator_id": tested_id,
            "created": created}


@router.post("/library/generators/{gen_id}/run-test")
async def run_test(request: Request, gen_id: str) -> dict:
    cat = _catalog(request)
    loc = cat.generator(gen_id)
    if loc is None:
        raise HTTPException(404, f"Unknown generator id: {gen_id}")
    _, g = loc
    if not g.valid:
        raise HTTPException(
            400, {"code": "invalid_entry", "errors": g.errors})

    body = {}
    try:
        body = await request.json()
    except Exception:
        pass
    test_config = (body or {}).get("test") or {
        "type": "equidistribution", "method": "matricial"}

    db_path = request.app.state.settings.db_path
    tested_id, test_result_id, se, is_me, created = (
        await asyncio.to_thread(
            _run_test_sync, db_path, g, test_config))
    return {
        "library_id": g.id,
        "tested_generator_id": tested_id,
        "test_result_id": test_result_id,
        "test_type": test_config.get("type"),
        "se": se,
        "is_me": is_me,
        "instantiated": created,
    }


# ─── Sync workers ─────────────────────────────────────────────────────────

def _instantiate_sync(db_path: str, g: Generator) -> tuple[int, bool]:
    from regpoly.combinaison import Combinaison
    from regpoly.combinaison_build import build_combinaison_inputs

    with sync_connect(db_path) as conn:
        row = conn.execute(
            "SELECT id FROM tested_generator WHERE library_id = ?",
            (g.id,),
        ).fetchone()
        if row is not None:
            return int(row["id"]), False

        gen_lists, temperings = build_combinaison_inputs(
            g.components, g.Lmax)
        comb = Combinaison.CreateFromFiles(gen_lists, g.Lmax, temperings)
        comb.reset()

        cur = conn.execute(
            "INSERT INTO tested_generator"
            "(search_run_id, Lmax, k_g, J, library_id) "
            "VALUES (?, ?, ?, ?, ?)",
            (None, comb.Lmax, comb.k_g, comb.J, g.id),
        )
        tested_id = cur.lastrowid

        from regpoly.web.tasks.tempering import _find_primitive_id

        for j in range(comb.J):
            gen = comb[j]
            comp = comb.components[j]
            gen_id = _find_primitive_id(conn, gen)
            tempering = [
                {"type": t._type_name, **t._params} for t in comp.trans
            ]
            conn.execute(
                "INSERT INTO tested_generator_component"
                "(tested_gen_id, component_index, generator_id, "
                " family, L, k, all_params, tempering_params) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (tested_id, j, gen_id, gen.type_name, gen.L, gen.k,
                 json_dumps(gen.params),
                 json_dumps(tempering)),
            )
        conn.commit()
        return int(tested_id), True


def _run_test_sync(
    db_path: str, g: Generator, test_config: dict,
) -> tuple[int, int, int | None, bool, bool]:
    tested_id, created = _instantiate_sync(db_path, g)

    from regpoly.combinaison import Combinaison
    from regpoly.combinaison_build import build_combinaison_inputs
    from regpoly.web.tasks._test_build import build_test
    from regpoly.web.tasks.tempering import (
        _build_result_detail, _is_me)

    gen_lists, temperings = build_combinaison_inputs(
        g.components, g.Lmax)
    comb = Combinaison.CreateFromFiles(gen_lists, g.Lmax, temperings)
    comb.reset()
    test = build_test(test_config, g.Lmax)
    result = test.run(comb)

    se = getattr(result, "se", None)
    is_me = bool(_is_me(result)) if se is not None else False
    detail = _build_result_detail(result)

    with sync_connect(db_path) as conn:
        cur = conn.execute(
            "INSERT INTO test_result"
            "(tested_gen_id, test_type, test_config, "
            " se, is_me, secf, is_cf, score, detail) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (tested_id,
             test_config.get("type", "equidistribution"),
             json_dumps(test_config),
             se,
             1 if is_me else 0,
             None, None,
             float(se) if se is not None else None,
             json_dumps(detail)),
        )
        conn.commit()
        return (int(tested_id), int(cur.lastrowid),
                int(se) if se is not None else None,
                is_me, created)
