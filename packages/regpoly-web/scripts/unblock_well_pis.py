#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""One-shot DB fix: mark currently-pending PIS rows that match the
known WELL hang signature as errored, so the analysis scheduler stops
re-dispatching them into the same trap.

Scope (intentionally narrow):

  * family = 'WELLGen'
  * pis_computed_at IS NULL
  * structural_params match the {p:1, r:19, w:32}  k=607 fingerprint

To retry after the underlying C++ hang is fixed, clear
``pis_computed_at`` and ``pis_error`` on the affected rows.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _load_db_url() -> str:
    """Mirror deploy/run-web-local.sh: prefer ``REGPOLY_DB_URL`` if
    already in env; otherwise construct the host-side DSN from
    ``DB_PASSWORD`` in ``deploy/.env``."""
    if (url := os.environ.get("REGPOLY_DB_URL")):
        return url
    env_path = Path(__file__).resolve().parents[3] / "deploy" / ".env"
    pw: str | None = None
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            if k.strip() == "REGPOLY_DB_URL":
                return v.strip().strip('"').strip("'")
            if k.strip() == "DB_PASSWORD":
                pw = v.strip().strip('"').strip("'")
    if pw is None:
        sys.exit("Neither REGPOLY_DB_URL nor DB_PASSWORD found")
    return f"postgresql://regpoly:{pw}@127.0.0.1:5432/regpoly"


def main() -> None:
    import psycopg
    db_url = _load_db_url()
    msg = ("WELL PIS analysis hang under investigation (k=607, "
           "matrix-dependent infinite loop in PISCache). "
           "Clear pis_computed_at to retry once fixed.")
    with psycopg.connect(db_url) as conn, conn.cursor() as cur:
        cur.execute(
            """
            SELECT COUNT(*) FROM primitive_generator
            WHERE family = 'WELLGen'
              AND pis_computed_at IS NULL
              AND k = 607
            """
        )
        n = cur.fetchone()[0]
        print(f"matched {n} WELL k=607 row(s) with pis_computed_at IS NULL")
        if n == 0:
            return
        cur.execute(
            """
            UPDATE primitive_generator
               SET pis_error = %s, pis_computed_at = NOW()
             WHERE family = 'WELLGen'
               AND pis_computed_at IS NULL
               AND k = 607
            """,
            (msg,),
        )
        conn.commit()
        print(f"marked {cur.rowcount} row(s) errored")


if __name__ == "__main__":
    main()
