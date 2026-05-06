# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""GET /api/v2/generators/{id}/transition-matrix-coords — re-export of v1.

The shape of the response (binary uint32 pairs blob) and headers
(X-Matrix-Size, X-Coord-Count) are identical to v1. The v2 mount is
declared with a Pydantic-style empty response_model so the OpenAPI
schema lists the endpoint; the actual response is `application/octet-stream`.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import Response

from regpoly_web.database import json_loads

router = APIRouter()


@router.get(
    "/generators/{gen_id}/transition-matrix-coords",
    responses={
        200: {"content": {"application/octet-stream": {}}},
        404: {"description": "Generator not found"},
    },
)
async def v2_generator_transition_matrix_coords(
    request: Request, gen_id: int,
) -> Response:
    db = request.app.state.db
    async with db.execute(
        "SELECT family, L, k, all_params FROM primitive_generator "
        "WHERE id = ?", (gen_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Generator {gen_id} not found")

    k = int(row["k"])
    from regpoly.core.generator import Generator

    params = json_loads(row["all_params"]) or {}
    try:
        gen = Generator.create(row["family"], row["L"], **params)
        mat = gen.transition_matrix()
    except Exception as exc:
        raise HTTPException(422, f"Could not build matrix: {exc}")

    # Borrow the v1 helper for the blob layout.
    from regpoly_web.routes.generators import _matrix_coords_blob

    blob = _matrix_coords_blob(mat, k)
    return Response(
        content=blob,
        media_type="application/octet-stream",
        headers={
            "Cache-Control": "private, max-age=60",
            "X-Matrix-Size": str(k),
            "X-Coord-Count": str(len(blob) // 8),
        },
    )
