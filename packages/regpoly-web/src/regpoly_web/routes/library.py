# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""API for the published-generators library (paper-centric)."""

from __future__ import annotations

import asyncio
import json
import logging

from fastapi import APIRouter, HTTPException, Request
from regpoly.library import Catalog, Generator, Paper

from regpoly_web.database import json_dumps, json_loads, sync_connect
from regpoly_web.param_format import (
    format_gen_params,
    format_tempering_list,
)
from regpoly_web.routes.families import _KNOWN_TESTS, _markdown_to_html

router = APIRouter()
log = logging.getLogger(__name__)


# Methods whose recommendation indicates a reducible characteristic
# polynomial (so harase, which assumes primitive χ, must be rejected).
# Used as a guard against explicit user overrides; the rule itself
# lives in C++ Generator::default_test_method.
_REDUCIBLE_CHI_METHODS = frozenset({"notprimitive", "simd_notprimitive"})


def _default_test_methods_for(g: Generator) -> dict[str, str | None]:
    """Build the per-test-type default-method map for a catalog entry.

    Instantiates the generator (read-only — no DB writes) and asks the
    underlying C++ Generator for its recommendation per test type.
    Returns ``{test_type: method | None}`` covering every entry in
    ``_KNOWN_TESTS``. Synchronous; expensive (Berlekamp-Massey on χ
    for some non-SIMD scalar families), so callers wrap in
    asyncio.to_thread and consult ``_default_test_methods_cached``
    to amortise the cost across requests.
    """
    import regpoly._regpoly_cpp as _cpp
    from regpoly.core.combination import Combination
    from regpoly.core.combination_build import build_combinaison_inputs

    gen_lists, temperings = build_combinaison_inputs(g.components, g.Lmax)
    comb = Combination.CreateFromFiles(gen_lists, g.Lmax, temperings)
    comb.reset()

    if comb.J == 1:
        cpp_gen = comb[0]._cpp_gen
    else:
        gens = [comb[j]._cpp_gen for j in range(comb.J)]
        trans = [
            [t._cpp_trans for t in comp.trans if hasattr(t, '_cpp_trans')]
            for comp in comb.components
        ]
        cpp_gen = _cpp.CombinedGenerator(gens, trans, comb.L)

    return {t: cpp_gen.default_test_method(t) for t in _KNOWN_TESTS}


async def _default_test_methods_cached(
    request: Request, paper: Paper, g: Generator,
) -> dict[str, str | None]:
    """Process-local cache around _default_test_methods_for.

    Keyed by ``(gen_id, paper.source_mtime)`` — when the YAML file
    backing the paper is edited and reload_if_stale picks up the new
    mtime, the cache key changes and the entry is recomputed
    automatically. The cache lives on app.state for the lifetime of
    the FastAPI process; entries from older mtimes become unreachable
    once any caller re-asks (memory grows by O(catalog size) at most).
    """
    cache: dict[tuple[str, float], dict[str, str | None]] = (
        getattr(request.app.state, "_default_methods_cache", None)
        or {})
    request.app.state._default_methods_cache = cache
    key = (g.id, float(getattr(paper, "source_mtime", 0.0)))
    if key in cache:
        return cache[key]
    result = await asyncio.to_thread(_default_test_methods_for, g)
    cache[key] = result
    return result


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
    if g.valid:
        try:
            row["default_test_methods"] = await _default_test_methods_cached(
                request, paper, g)
        except (FileNotFoundError, KeyError, ValueError) as e:
            # Catalog/parameter issue: log and degrade gracefully so
            # rendering still succeeds. Other exception classes are
            # let through so they surface as 500s instead of becoming
            # invisible "no recommendation" UI states.
            log.warning(
                "default_test_methods unavailable for %s: %s", g.id, e,
            )
            row["default_test_methods"] = {
                t: None for t in _KNOWN_TESTS}
    else:
        row["default_test_methods"] = {t: None for t in _KNOWN_TESTS}
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
    db_path = request.app.state.settings.db_url
    tested_id, created = await asyncio.to_thread(
        _instantiate_sync, db_path, g)
    return {"library_id": g.id,
            "tested_generator_id": tested_id,
            "created": created}


@router.post("/library/generators/{gen_id}/run-test")
async def run_test(request: Request, gen_id: str) -> dict:
    """Kick off a test on a library generator in a worker process.

    Returns a job_id immediately; the client polls
    ``GET /library/run-test-jobs/{job_id}`` until the status is
    ``completed`` or ``failed``.  This matches how primitive- and
    tempering-search runs are dispatched (cf. routes/primitive_search.py)
    and prevents the C++ equidistribution computation from blocking the
    HTTP request thread.
    """
    cat = _catalog(request)
    loc = cat.generator(gen_id)
    if loc is None:
        raise HTTPException(404, f"Unknown generator id: {gen_id}")
    paper, g = loc
    if not g.valid:
        raise HTTPException(
            400, {"code": "invalid_entry", "errors": g.errors})

    body = {}
    try:
        body = await request.json()
    except Exception:
        pass
    # When no test config is supplied, default to equidistribution and
    # let the worker resolve the method via the generator's
    # default_test_method (no client-side method assumption).
    test_config = (body or {}).get("test") or {"type": "equidistribution"}

    # Reject harase on reducible χ at the API boundary — the primal-
    # lattice search has no finite upper bound when χ is reducible and
    # would peg a worker forever. The C++ Generator::default_test_method
    # is the canonical signal: a recommendation of "notprimitive" or
    # "simd_notprimitive" implies reducible χ. Goes through the cached
    # helper so the detail-page fetch and this check share state.
    if (test_config.get("type") == "equidistribution"
            and test_config.get("method") == "harase"):
        try:
            defaults = await _default_test_methods_cached(request, paper, g)
            recommended = defaults.get("equidistribution")
        except (FileNotFoundError, KeyError, ValueError) as e:
            log.warning(
                "harase reducibility check skipped for %s: %s", g.id, e,
            )
            recommended = None
        if recommended in _REDUCIBLE_CHI_METHODS:
            raise HTTPException(
                400,
                {"code": "method_not_supported",
                 "method": "harase",
                 "recommended": recommended,
                 "message": "harase does not terminate on this generator "
                            "(reducible characteristic polynomial). "
                            "Pick the recommended method instead."},
            )

    # Insert a pending library_test_run row; the worker scheduler
    # picks it up via FOR UPDATE SKIP LOCKED. Returns the row id as
    # job_id (opaque to the front-end).
    db = request.app.state.db
    cur = await db.execute(
        """
        INSERT INTO library_test_run
            (library_id, test_type, method, test_config, status)
        VALUES (?, ?, ?, ?::jsonb, 'pending')
        """,
        (
            g.id,
            test_config.get("type"),
            test_config.get("method"),
            json.dumps(test_config),
        ),
    )
    job_id = cur.lastrowid
    return {
        "job_id": str(job_id),
        "status": "pending",
        "library_id": g.id,
        "test_type": test_config.get("type"),
        "method": test_config.get("method"),
    }


@router.get("/library/run-test-jobs/{job_id}")
async def run_test_status(request: Request, job_id: str) -> dict:
    db = request.app.state.db
    try:
        row_id = int(job_id)
    except ValueError:
        raise HTTPException(404, f"Unknown job id: {job_id}")
    async with db.execute(
        "SELECT status, result, error_code, error_message "
        "FROM library_test_run WHERE id = ?",
        (row_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Unknown job id: {job_id}")
    status = row["status"]
    if status == "failed":
        raise HTTPException(400, {
            "code": row["error_code"] or "test_failed",
            "message": row["error_message"] or "unknown failure",
        })
    payload: dict = {"job_id": job_id, "status": status}
    if status == "completed":
        result = json_loads(row["result"]) or {}
        payload.update(result)
    return payload


# ─── Sync workers ─────────────────────────────────────────────────────────

def _instantiate_sync(db_path: str, g: Generator) -> tuple[int, bool]:
    from regpoly.core.combination import Combination
    from regpoly.core.combination_build import build_combinaison_inputs

    with sync_connect(db_path) as conn:
        row = conn.execute(
            "SELECT id FROM tested_generator WHERE library_id = ?",
            (g.id,),
        ).fetchone()
        if row is not None:
            return int(row["id"]), False

        gen_lists, temperings = build_combinaison_inputs(
            g.components, g.Lmax)
        comb = Combination.CreateFromFiles(gen_lists, g.Lmax, temperings)
        comb.reset()

        cur = conn.execute(
            "INSERT INTO tested_generator"
            "(search_run_id, Lmax, k_g, J, library_id) "
            "VALUES (?, ?, ?, ?, ?)",
            (None, comb.Lmax, comb.k_g, comb.J, g.id),
        )
        tested_id = cur.lastrowid

        from regpoly_web.tasks.tempering import _find_primitive_id

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


