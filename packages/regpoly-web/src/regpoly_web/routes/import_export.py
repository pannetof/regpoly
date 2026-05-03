"""YAML import/export endpoints (migration and backwards-compatibility)."""

from __future__ import annotations

import os
from pathlib import Path

import yaml
from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel

from regpoly.core.generator import Generator, resolve_family
from regpoly.web.database import json_dumps, json_loads


router = APIRouter()


class ImportRequest(BaseModel):
    file_path: str


class ImportDirRequest(BaseModel):
    directory: str


@router.post("/import/generators")
async def import_generators_file(
    request: Request, body: ImportRequest
) -> dict:
    return await _import_one_generators_file(request, Path(body.file_path))


@router.post("/import/generators-dir")
async def import_generators_dir(
    request: Request, body: ImportDirRequest
) -> dict:
    """Recursively import every generators YAML file under a directory."""
    root = Path(body.directory)
    if not root.is_dir():
        raise HTTPException(404, f"Not a directory: {root}")

    imported_files = 0
    imported_rows = 0
    errors: list[dict] = []

    for yml in sorted(root.rglob("*.yaml")):
        try:
            res = await _import_one_generators_file(request, yml)
            imported_files += 1
            imported_rows += res.get("inserted", 0)
        except Exception as exc:
            errors.append({"file": str(yml), "error": str(exc)})

    return {
        "imported_files": imported_files,
        "inserted": imported_rows,
        "errors": errors,
    }


async def _import_one_generators_file(
    request: Request, path: Path
) -> dict:
    if not path.exists():
        raise HTTPException(404, f"File not found: {path}")
    with open(path) as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict) or "family" not in data \
            or "generators" not in data:
        raise HTTPException(
            400,
            f"Not a generators YAML file (missing 'family'/'generators'): {path}",
        )

    family_raw = data["family"]
    family = resolve_family(family_raw)
    common = data.get("common") or {}
    generators = data.get("generators") or []

    db = request.app.state.db

    # Already imported?
    async with db.execute(
        "SELECT id FROM yaml_import WHERE file_path = ?",
        (str(path.resolve()),),
    ) as cur:
        prev = await cur.fetchone()
    if prev is not None:
        return {"file": str(path), "inserted": 0, "skipped": len(generators),
                "note": "already imported"}

    inserted = 0
    duplicates = 0
    failures: list[dict] = []

    for entry in generators:
        merged = {**common, **entry}
        L = merged.pop("L", None)
        if L is None:
            # Try to infer L from a 'w' field (common default)
            L = merged.get("w", 32)

        # Probe to get k; if creation fails, skip row with an error note.
        try:
            gen = Generator.create(family_raw, L, **merged)
        except Exception as exc:
            failures.append({"entry": dict(entry), "error": str(exc)})
            continue

        structural = gen.structural_params()
        search_only = {
            k: v for k, v in merged.items() if k not in structural
        }

        try:
            await db.execute(
                """
                INSERT OR IGNORE INTO primitive_generator
                    (search_run_id, family, L, k,
                     structural_params, search_params, all_params,
                     found_at_try)
                VALUES (NULL, ?, ?, ?, ?, ?, ?, NULL)
                """,
                (
                    family,
                    gen.L,
                    gen.k,
                    json_dumps(structural),
                    json_dumps(search_only),
                    json_dumps({**structural, **search_only}),
                ),
            )
            # rowcount reflects whether the INSERT actually wrote a row.
            if db.total_changes:
                inserted += 1
        except Exception as exc:
            failures.append({"entry": dict(entry), "error": str(exc)})

    await db.execute(
        "INSERT INTO yaml_import(file_path, import_type, row_count) "
        "VALUES (?, 'generators', ?)",
        (str(path.resolve()), inserted),
    )
    await db.commit()

    return {
        "file": str(path),
        "family": family,
        "inserted": inserted,
        "duplicates": duplicates,
        "failures": failures,
    }


# ── Export ──────────────────────────────────────────────────────────────

@router.get("/export/generators", response_class=PlainTextResponse)
async def export_generators_yaml(
    request: Request,
    family: str,
    structural_params: str | None = None,
) -> str:
    """Export all generators of a family (optionally filtered by structural
    params JSON) as a regpoly generators YAML file."""
    db = request.app.state.db
    where = ["family = ?"]
    params: list = [family]
    if structural_params:
        where.append("structural_params = ?")
        params.append(structural_params)

    where_clause = " AND ".join(where)
    async with db.execute(
        f"SELECT structural_params, search_params FROM primitive_generator "
        f"WHERE {where_clause} ORDER BY id",
        params,
    ) as cur:
        rows = await cur.fetchall()

    if not rows:
        raise HTTPException(404, "No generators match the filter")

    structural = json_loads(rows[0][0])
    generators = [json_loads(r[1]) for r in rows]

    data = {
        "family": family,
        "common": structural,
        "generators": generators,
    }
    return yaml.dump(data, default_flow_style=False, sort_keys=False)


@router.get(
    "/export/tested-generators/{tg_id}",
    response_class=PlainTextResponse,
)
async def export_tested_generator(request: Request, tg_id: int) -> str:
    """Export a tested generator as a regpoly tested-generator YAML file."""
    from regpoly.web.routes.tested_generators import _fetch_tested
    db = request.app.state.db
    tg = await _fetch_tested(db, tg_id)
    if tg is None:
        raise HTTPException(404, f"Tested generator {tg_id} not found")

    def _component(c):
        gen_data = {"family": c["family"], "L": c["L"], **c["all_params"]}
        out = {"generator": gen_data}
        if c["tempering_params"]:
            out["tempering"] = c["tempering_params"]
        return out

    # Build results dict matching the save_tested_generator format
    results = {}
    for r in tg["results"]:
        if r["test_type"] == "equidistribution":
            eq = {"se": r["se"]}
            if r.get("is_me"):
                eq["status"] = "ME"
            ecart = (r.get("detail") or {}).get("ecart")
            if ecart:
                eq["ecart"] = {int(k): v for k, v in ecart.items()}
            results["equidistribution"] = eq

    if tg["J"] == 1:
        data = _component(tg["components"][0])
        data["results"] = results
    else:
        data = {
            "components": [_component(c) for c in tg["components"]],
            "results": results,
        }

    return yaml.dump(data, default_flow_style=False, sort_keys=False)
