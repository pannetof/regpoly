"""POST /api/v2/import/generators/preview — dry-run import.

Reads the YAML payload (uploaded file or path on disk), reports what
would be added / skipped / conflicts, without touching the DB.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from fastapi import APIRouter, HTTPException, Request, UploadFile, File
from pydantic import BaseModel

router = APIRouter()


class _Conflict(BaseModel):
    file_path: str
    last_imported_at: str | None


class PreviewResponse(BaseModel):
    would_add: int
    would_skip: int
    conflicts: list[_Conflict]


async def _previously_imported(db, file_path: str | None) -> str | None:
    if not file_path:
        return None
    async with db.execute(
        "SELECT created_at FROM yaml_import WHERE file_path = ?",
        (file_path,),
    ) as cur:
        row = await cur.fetchone()
    return row["created_at"] if row else None


def _count_in_yaml(data: dict) -> int:
    if not isinstance(data, dict):
        return 0
    gens = data.get("generators") or []
    return len(gens)


@router.post(
    "/import/generators/preview", response_model=PreviewResponse,
)
async def v2_import_preview(
    request: Request, file: UploadFile | None = File(default=None),
) -> PreviewResponse:
    if file is None:
        raise HTTPException(400, "Upload a YAML file as 'file'")
    blob = await file.read()
    try:
        data = yaml.safe_load(blob)
    except yaml.YAMLError as exc:
        raise HTTPException(400, f"YAML parse error: {exc}")

    if not isinstance(data, dict):
        raise HTTPException(400, "Top-level YAML must be a mapping")
    n = _count_in_yaml(data)
    if not data.get("family") or not data.get("generators"):
        # The v1 import requires both, so a preview tells the user
        # it would skip everything.
        return PreviewResponse(would_add=0, would_skip=n, conflicts=[])

    # The dry-run check uses file_path = the uploaded filename;
    # there's no on-disk path for an upload. Check by name only.
    db = request.app.state.db
    fname = file.filename or "<upload>"
    prev = await _previously_imported(db, fname)
    conflicts: list[_Conflict] = []
    would_add = n
    would_skip = 0
    if prev is not None:
        conflicts.append(_Conflict(file_path=fname, last_imported_at=prev))
        would_add = 0
        would_skip = n

    return PreviewResponse(
        would_add=would_add, would_skip=would_skip, conflicts=conflicts,
    )
