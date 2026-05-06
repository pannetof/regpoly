# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Generator script for tests/fixtures/openapi.v2.json.

Run via:
    uv run python packages/regpoly-web/tests/_gen_openapi_snapshot.py

Boots a TestClient against an in-memory app and writes the v2 path
slice of /openapi.json to the committed fixture. Updates require a
deliberate refresh — `test_v2_openapi_matches_committed_snapshot`
diffs against this file in CI.
"""

from __future__ import annotations

import json
from pathlib import Path

from fastapi.testclient import TestClient

from regpoly_web.app import create_app
from regpoly_web.config import Settings


def main() -> None:
    here = Path(__file__).resolve().parent
    fixture = here / "fixtures" / "openapi.v2.json"
    fixture.parent.mkdir(parents=True, exist_ok=True)

    settings = Settings(db_path=":memory:", reload=False, pool_size=1)
    app = create_app(settings)
    with TestClient(app) as c:
        spec = c.get("/openapi.json").json()
    paths = {p: ops for p, ops in spec.get("paths", {}).items()
             if p.startswith("/api/v2/")}
    snapshot = {"openapi": spec.get("openapi"), "paths": paths}
    fixture.write_text(json.dumps(snapshot, indent=2, sort_keys=True))
    print(f"wrote {fixture} ({len(paths)} v2 paths)")


if __name__ == "__main__":
    main()
