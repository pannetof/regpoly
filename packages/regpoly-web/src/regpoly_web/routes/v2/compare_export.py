# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""GET /api/v2/generators/compare?ids= and /generators/export?ids=&fmt=.

The compare endpoint returns one row per id (with `null` placeholders
for missing ids); the export endpoint emits CSV or JSON for the
selected ids only — used by the "Copy selection" toolbar item.
"""

from __future__ import annotations

import csv
import io
import json

from fastapi import APIRouter, Query, Request
from fastapi.responses import PlainTextResponse, JSONResponse

from regpoly_web.database import json_loads

router = APIRouter()


def _csv_safe(v) -> str:
    """Escape leading characters that Excel/LibreOffice would interpret
    as a formula (CWE-1236). Cells starting with `=`, `+`, `-`, `@`,
    tab, or carriage-return get a single-quote prefix."""
    s = "" if v is None else str(v)
    if s and s[0] in ("=", "+", "-", "@", "\t", "\r"):
        return "'" + s
    return s


def _parse_ids(raw: str) -> list[int]:
    out: list[int] = []
    for tok in (raw or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            out.append(int(tok))
        except ValueError:
            continue
    return out


async def _rows_by_id(db, ids: list[int]) -> dict[int, dict]:
    """Hydrate primitive_generator rows by id, including pis_gaps and
    pis_se for chart overlay on the compare page."""
    if not ids:
        return {}
    placeholders = ",".join("?" for _ in ids)
    async with db.execute(
        f"SELECT * FROM primitive_generator WHERE id IN ({placeholders})",
        ids,
    ) as cur:
        rows = await cur.fetchall()
    out: dict[int, dict] = {}
    for r in rows:
        keys = r.keys()
        out[int(r["id"])] = {
            "id": int(r["id"]),
            "search_run_id": r["search_run_id"],
            "family": r["family"],
            "L": r["L"],
            "k": r["k"],
            "found_at_try": r["found_at_try"]
                if "found_at_try" in keys else None,
            "char_poly": r["char_poly"]
                if "char_poly" in keys else None,
            "hamming_weight": r["hamming_weight"]
                if "hamming_weight" in keys else None,
            "library_id": r["library_id"]
                if "library_id" in keys else None,
            "pis_se": r["pis_se"] if "pis_se" in keys else None,
            "pis_gaps": (json_loads(r["pis_gaps"])
                         if "pis_gaps" in keys and r["pis_gaps"]
                         else None),
            "all_params": json_loads(r["all_params"]) or {},
        }
    return out


@router.get("/generators/compare")
async def v2_compare_generators(
    request: Request, ids: str = Query(""),
) -> dict:
    db = request.app.state.db
    id_list = _parse_ids(ids)
    by_id = await _rows_by_id(db, id_list)
    rows = [by_id.get(i) for i in id_list]
    return {"rows": rows, "ids": id_list}


async def _tested_rows_by_id(db, ids: list[int]) -> dict[int, dict]:
    if not ids:
        return {}
    placeholders = ",".join("?" for _ in ids)
    async with db.execute(
        f"SELECT tg.id, tg.search_run_id, tg.Lmax, tg.k_g, tg.J, "
        f"       tg.library_id "
        f"FROM tested_generator AS tg "
        f"WHERE tg.id IN ({placeholders})",
        ids,
    ) as cur:
        rows = await cur.fetchall()
    out: dict[int, dict] = {}
    for r in rows:
        tid = int(r["id"])
        out[tid] = {
            "id": tid,
            "search_run_id": r["search_run_id"],
            "Lmax": r["Lmax"],
            "k_g": r["k_g"],
            "J": r["J"],
            "library_id": r["library_id"],
        }
    if out:
        async with db.execute(
            f"SELECT tested_gen_id, se, ecart_json FROM equidistribution_result "
            f"WHERE tested_gen_id IN ({placeholders})",
            ids,
        ) as cur:
            er_rows = await cur.fetchall()
        for er in er_rows:
            tid = int(er["tested_gen_id"])
            if tid not in out:
                continue
            out[tid]["equidistribution"] = {
                "se": er["se"],
                "ecart": (json_loads(er["ecart_json"])
                          if er["ecart_json"] else None),
            }
    return out


@router.get("/tested-generators/compare")
async def v2_compare_tested_generators(
    request: Request, ids: str = Query(""),
) -> dict:
    """Per-id compare payload for tested generators (parity with
    /api/v2/generators/compare). Missing ids appear as null."""
    db = request.app.state.db
    id_list = _parse_ids(ids)
    by_id = await _tested_rows_by_id(db, id_list)
    rows = [by_id.get(i) for i in id_list]
    return {"rows": rows, "ids": id_list}


@router.get("/library/compare")
async def v2_compare_library_entries(
    request: Request, ids: str = Query(""),
) -> dict:
    """Per-id compare payload for library entries. Library ids are
    string-typed (not integer), and we look them up against the
    in-memory catalog. Missing ids appear as null."""
    catalog = request.app.state.library
    id_list = [tok.strip() for tok in (ids or "").split(",") if tok.strip()]
    # Build a flat lookup once.
    by_id: dict[str, dict] = {}
    for paper in catalog.papers():
        for gen in getattr(paper, "generators", []) or []:
            gid = getattr(gen, "id", None)
            if gid is None:
                continue
            by_id[gid] = {
                "id": gid,
                "paper_id": paper.id,
                "family": getattr(gen, "family", None),
                "L": getattr(gen, "L", None),
                "k": getattr(gen, "k", None),
            }
    rows = [by_id.get(i) for i in id_list]
    return {"rows": rows, "ids": id_list}


@router.get("/generators/export")
async def v2_export_generators(
    request: Request,
    ids: str = Query(""),
    fmt: str = Query("csv"),
):
    db = request.app.state.db
    id_list = _parse_ids(ids)
    by_id = await _rows_by_id(db, id_list)
    # Skip missing ids (vs compare which preserves slot order).
    rows = [by_id[i] for i in id_list if i in by_id]

    if fmt == "json":
        return JSONResponse(content=rows)

    # CSV — flatten all_params into individual columns.
    buf = io.StringIO()
    base_cols = [
        "id", "family", "L", "k", "found_at_try", "char_poly",
        "hamming_weight", "library_id", "search_run_id",
    ]
    # Sort param keys alphabetically so column order is stable across
    # requests with the same id-set in different orders.
    seen: set[str] = set()
    for r in rows:
        for k in (r.get("all_params") or {}):
            seen.add(k)
    param_keys = sorted(seen)
    writer = csv.writer(buf)
    writer.writerow(base_cols + ["p_" + k for k in param_keys])
    for r in rows:
        params = r.get("all_params") or {}
        writer.writerow(
            [_csv_safe(r.get(c, "")) for c in base_cols]
            + [_csv_safe(params.get(k, "")) for k in param_keys]
        )
    return PlainTextResponse(
        content=buf.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": "inline; filename=generators.csv"},
    )


@router.get("/tested-generators/export")
async def v2_export_tested_generators(
    request: Request,
    ids: str = Query(""),
    fmt: str = Query("csv"),
):
    """Selection-scoped export for combined (tested) generators. Mirrors
    /generators/export — accepts a comma-separated id list and emits CSV
    or JSON. Missing ids are silently dropped."""
    db = request.app.state.db
    id_list = _parse_ids(ids)
    by_id = await _tested_rows_by_id(db, id_list)
    rows = [by_id[i] for i in id_list if i in by_id]

    if fmt == "json":
        return JSONResponse(content=rows)

    buf = io.StringIO()
    base_cols = ["id", "search_run_id", "J", "k_g", "Lmax", "library_id",
                 "se"]
    writer = csv.writer(buf)
    writer.writerow(base_cols)
    for r in rows:
        eq = r.get("equidistribution") or {}
        writer.writerow([
            _csv_safe(r.get("id")),
            _csv_safe(r.get("search_run_id")),
            _csv_safe(r.get("J")),
            _csv_safe(r.get("k_g")),
            _csv_safe(r.get("Lmax")),
            _csv_safe(r.get("library_id")),
            _csv_safe(eq.get("se", "")),
        ])
    return PlainTextResponse(
        content=buf.getvalue(),
        media_type="text/csv",
        headers={
            "Content-Disposition": "inline; filename=combined-generators.csv"
        },
    )
