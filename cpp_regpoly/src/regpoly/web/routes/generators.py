"""API endpoints for primitive generators stored in the database."""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Request

from regpoly.web.database import json_loads


router = APIRouter()


def _row_to_generator(row) -> dict:
    keys = row.keys() if hasattr(row, "keys") else []
    def opt(name):
        return row[name] if name in keys else None
    return {
        "id": row["id"],
        "search_run_id": row["search_run_id"],
        "family": row["family"],
        "L": row["L"],
        "k": row["k"],
        "structural_params": json_loads(row["structural_params"]),
        "search_params": json_loads(row["search_params"]),
        "all_params": json_loads(row["all_params"]),
        "found_at_try": row["found_at_try"],
        "created_at": row["created_at"],
        "char_poly": opt("char_poly"),
        "hamming_weight": opt("hamming_weight"),
        "pis_gaps": json_loads(opt("pis_gaps")),
        "pis_se": opt("pis_se"),
        "pis_elapsed": opt("pis_elapsed"),
        "pis_computed_at": opt("pis_computed_at"),
        "pis_error": opt("pis_error"),
    }


@router.get("/generators")
async def list_generators(
    request: Request,
    family: str | None = None,
    k: int | None = None,
    k_min: int | None = None,
    k_max: int | None = None,
    search_run_id: int | None = None,
    available_only: bool = True,
    page: int = Query(1, ge=1),
    per_page: int = Query(50, ge=1, le=500),
) -> dict:
    """List primitive generators.

    By default only returns generators whose analysis has completed
    successfully (char poly + PIS gaps available).  Pass
    `available_only=false` to include queued/pending/errored rows too.

    Any query param prefixed with `p_` filters by an exact match on the
    corresponding parameter inside `all_params` (via json_extract).
    Example: `?family=TGFSRGen&p_w=32&p_r=3` matches generators whose
    all_params has w=32 and r=3.
    """
    db = request.app.state.db

    where = []
    params: list = []
    if family:
        where.append("family = ?")
        params.append(family)
    if k is not None:
        where.append("k = ?")
        params.append(k)
    if k_min is not None:
        where.append("k >= ?")
        params.append(k_min)
    if k_max is not None:
        where.append("k <= ?")
        params.append(k_max)
    if search_run_id is not None:
        where.append("search_run_id = ?")
        params.append(search_run_id)
    if available_only:
        where.append(
            "pis_computed_at IS NOT NULL "
            "AND pis_error IS NULL "
            "AND char_poly IS NOT NULL"
        )

    for key, val in request.query_params.multi_items():
        if not key.startswith("p_") or val == "":
            continue
        pname = key[2:]
        if not pname.replace("_", "").isalnum():
            continue  # reject anything weird as a safety measure
        where.append(
            f"json_extract(all_params, '$.{pname}') = ?"
        )
        params.append(_coerce_param_value(val))

    where_clause = ("WHERE " + " AND ".join(where)) if where else ""

    count_sql = f"SELECT COUNT(*) FROM primitive_generator {where_clause}"
    async with db.execute(count_sql, params) as cur:
        row = await cur.fetchone()
        total = row[0]

    offset = (page - 1) * per_page
    list_sql = (
        f"SELECT * FROM primitive_generator {where_clause} "
        f"ORDER BY id DESC LIMIT ? OFFSET ?"
    )
    async with db.execute(list_sql, [*params, per_page, offset]) as cur:
        rows = await cur.fetchall()

    return {
        "items": [_row_to_generator(r) for r in rows],
        "total": total,
        "page": page,
        "per_page": per_page,
    }


@router.get("/generators/families")
async def generator_families(request: Request) -> list[dict]:
    """Families that have at least one available generator in the DB.

    Used by the /generators UI to populate the Family dropdown with
    only the families the user can actually pick from.
    """
    db = request.app.state.db
    sql = (
        "SELECT family, COUNT(*) AS n "
        "FROM primitive_generator "
        "WHERE pis_computed_at IS NOT NULL "
        "AND pis_error IS NULL "
        "AND char_poly IS NOT NULL "
        "GROUP BY family ORDER BY family"
    )
    async with db.execute(sql) as cur:
        rows = await cur.fetchall()
    return [{"name": r[0], "count": r[1]} for r in rows]


@router.get("/generators/param-values")
async def generator_param_values(
    request: Request, name: str, family: str | None = None,
) -> dict:
    """Return distinct values `name` takes for rows matching the active
    filters — *excluding* any filter on `name` itself so the user can
    still switch between its options.

    `family` is optional.  If omitted, results span all families (useful
    for columns like `k` and `L` that exist on every row).

    Active filters come from the standard query params (`k`, `p_*`).
    Special names `k` and `L` read the dedicated columns instead of
    json_extract(all_params).
    """
    db = request.app.state.db

    if name not in ("k", "L") and not name.replace("_", "").isalnum():
        from fastapi import HTTPException
        raise HTTPException(400, "invalid parameter name")

    where = []
    params: list = []
    if family:
        where.append("family = ?")
        params.append(family)

    # Match the same availability default as /generators
    avail_raw = request.query_params.get("available_only", "true").lower()
    if avail_raw not in ("0", "false", "no"):
        where.append(
            "pis_computed_at IS NOT NULL "
            "AND pis_error IS NULL "
            "AND char_poly IS NOT NULL"
        )

    # Apply k filter except when computing values for k itself
    k_raw = request.query_params.get("k")
    if name != "k" and k_raw not in (None, ""):
        try:
            where.append("k = ?")
            params.append(int(k_raw))
        except ValueError:
            pass

    # Apply each p_<x> filter except when x == name
    for key, val in request.query_params.multi_items():
        if not key.startswith("p_") or val == "":
            continue
        pname = key[2:]
        if pname == name:
            continue
        if not pname.replace("_", "").isalnum():
            continue
        where.append(
            f"json_extract(all_params, '$.{pname}') = ?"
        )
        params.append(_coerce_param_value(val))

    where_clause = (" WHERE " + " AND ".join(where)) if where else ""
    if name in ("k", "L"):
        sql = (
            f"SELECT DISTINCT {name} AS v FROM primitive_generator"
            f"{where_clause} ORDER BY v"
        )
    else:
        v_not_null = (
            " AND v IS NOT NULL" if where_clause
            else " WHERE v IS NOT NULL"
        )
        sql = (
            f"SELECT DISTINCT json_extract(all_params, '$.{name}') AS v "
            f"FROM primitive_generator{where_clause}{v_not_null} "
            f"ORDER BY v"
        )
    async with db.execute(sql, params) as cur:
        rows = await cur.fetchall()
    return {"family": family, "name": name,
            "values": [r[0] for r in rows]}


def _coerce_param_value(val: str):
    """Convert a string query-param into int/hex/raw string for comparison."""
    s = val.strip()
    if s.lstrip("-").isdigit():
        return int(s)
    if s.startswith("0x") or s.startswith("0X"):
        try:
            return int(s, 16)
        except ValueError:
            pass
    return s


@router.get("/generators/{gen_id}")
async def get_generator(request: Request, gen_id: int) -> dict:
    db = request.app.state.db
    async with db.execute(
        "SELECT * FROM primitive_generator WHERE id = ?", (gen_id,)
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Generator {gen_id} not found")
    return _row_to_generator(row)


@router.delete("/generators/{gen_id}")
async def delete_generator(request: Request, gen_id: int) -> dict:
    db = request.app.state.db
    cur = await db.execute(
        "DELETE FROM primitive_generator WHERE id = ?", (gen_id,)
    )
    await db.commit()
    if cur.rowcount == 0:
        raise HTTPException(404, f"Generator {gen_id} not found")
    return {"deleted": gen_id}


@router.get("/generators/{gen_id}/transition-matrix-coords")
async def generator_transition_matrix_coords(
    request: Request, gen_id: int,
):
    """Return set-bit coordinates of the transition matrix as a binary
    blob of little-endian uint32 pairs: [row0, col0, row1, col1, ...].

    The client draws these on a canvas (one fillRect per pair).  This
    handles arbitrarily large k because only the set bits are sent.
    """
    db = request.app.state.db
    async with db.execute(
        "SELECT family, L, k, all_params FROM primitive_generator "
        "WHERE id = ?", (gen_id,)
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Generator {gen_id} not found")

    k = int(row["k"])

    from regpoly.generateur import Generateur
    from regpoly.web.database import json_loads as _jl
    params = _jl(row["all_params"]) or {}
    gen = Generateur.create(row["family"], row["L"], **params)
    mat = gen.transition_matrix()

    blob = _matrix_coords_blob(mat, k)

    from fastapi.responses import Response
    return Response(
        content=blob,
        media_type="application/octet-stream",
        headers={
            "Cache-Control": "private, max-age=60",
            "X-Matrix-Size": str(k),
            "X-Coord-Count": str(len(blob) // 8),  # bytes per pair
        },
    )


def _matrix_coords_blob(mat, k: int) -> bytes:
    """Emit all set bits as little-endian uint32 (row, col) pairs.

    BitVect convention: column j of row i is the bit at Python bit
    position (k-1-j) inside mat._rows[i].  We iterate only over the set
    bits of each row, so the cost is O(total set bits).
    """
    import array
    import sys

    coords = array.array("I")
    for i in range(k):
        row = mat._rows[i]
        while row:
            # lowest set bit
            lsb = row & -row
            bit_pos = lsb.bit_length() - 1
            col = (k - 1) - bit_pos
            coords.append(i)
            coords.append(col)
            row ^= lsb
    # array('I') is native endian — force little-endian for the browser
    if sys.byteorder == "big":
        coords.byteswap()
    return coords.tobytes()


