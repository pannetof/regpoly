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
        out[int(r["id"])] = {
            "id": int(r["id"]),
            "search_run_id": r["search_run_id"],
            "family": r["family"],
            "L": r["L"],
            "k": r["k"],
            "found_at_try": r["found_at_try"]
                if "found_at_try" in r.keys() else None,
            "char_poly": r["char_poly"]
                if "char_poly" in r.keys() else None,
            "hamming_weight": r["hamming_weight"]
                if "hamming_weight" in r.keys() else None,
            "library_id": r["library_id"]
                if "library_id" in r.keys() else None,
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
    param_keys: list[str] = []
    seen: set[str] = set()
    for r in rows:
        for k in (r.get("all_params") or {}):
            if k not in seen:
                seen.add(k)
                param_keys.append(k)
    writer = csv.writer(buf)
    writer.writerow(base_cols + ["p_" + k for k in param_keys])
    for r in rows:
        params = r.get("all_params") or {}
        writer.writerow(
            [r.get(c, "") for c in base_cols]
            + [params.get(k, "") for k in param_keys]
        )
    return PlainTextResponse(
        content=buf.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": "inline; filename=generators.csv"},
    )
