#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Standalone repro for the WELL k=607 PIS analysis hang.

Loads one (or a few) ``primitive_generator`` rows by id, reconstructs
the Generator from its persisted ``all_params``, and runs
``analyze_single_generator`` foreground — single-process, no pool, no
scheduler. So you can:

  * py-spy dump --pid $(pgrep -f repro_pis_hang.py) to see where it
    loops without touching the running web app.
  * Ctrl-C cleanly: foreground means the SIGINT lands on the worker.
  * Run under gdb if you want a native stack: ``gdb --args uv run …``.

Usage:

    cd /path/to/regpoly_monorepo
    # Loads DSN from deploy/.env (REGPOLY_DB_URL); same DB as the web.
    uv run python packages/regpoly-web/scripts/repro_pis_hang.py 12

    # Multiple ids:
    uv run python packages/regpoly-web/scripts/repro_pis_hang.py 12 13 14

The script prints elapsed time and either the PIS summary or the
exception. If it doesn't print anything within a few seconds for a row
the web pool also got stuck on, you've reproduced the hang.
"""

from __future__ import annotations

import json
import os
import sys
import time
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


def _fetch_row(db_url: str, gen_id: int) -> dict:
    import psycopg
    with psycopg.connect(db_url) as conn, conn.cursor() as cur:
        cur.execute(
            "SELECT family, L, k, all_params FROM primitive_generator "
            "WHERE id = %s",
            (gen_id,),
        )
        row = cur.fetchone()
    if row is None:
        sys.exit(f"primitive_generator #{gen_id} not found")
    family, L, k, all_params = row
    if isinstance(all_params, str):
        all_params = json.loads(all_params)
    return {"family": family, "L": L, "k": k, "all_params": all_params}


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(f"usage: {sys.argv[0]} <gen_id> [gen_id ...]")
    ids = [int(x) for x in sys.argv[1:]]

    from regpoly.analyses.pis import analyze_single_generator
    from regpoly.core.generator import Generator

    db_url = _load_db_url()
    print(f"# pid={os.getpid()}  db={db_url.split('@')[-1]}")
    print(f"# attach a profiler with:  py-spy dump --pid {os.getpid()}")
    print()

    keys_of_interest = ("w", "r", "p", "nb_terms", "t", "q", "step")
    # PYTHONUNBUFFERED isn't auto-set when launched via nohup, so flush
    # at every checkpoint — otherwise an early hang shows up as an empty
    # log file and the user can't tell which call we're stuck in.
    def say(msg: str) -> None:
        print(msg, flush=True)

    # Bisect by stage so a hang under nohup can be pinned down by
    # tailing the log: Generator.create / gen.char_poly / PISCache.
    from regpoly_cpp import _regpoly_cpp as _cpp  # noqa: E402
    from regpoly.core.combination import Combination  # noqa: E402

    for gen_id in ids:
        meta = _fetch_row(db_url, gen_id)
        struct = {k: meta["all_params"].get(k) for k in keys_of_interest}
        say(f"=== #{gen_id}  family={meta['family']}  "
            f"L={meta['L']}  k={meta['k']}  struct={struct} ===")
        t0 = time.monotonic()
        try:
            say("  [1/5] Generator.create …")
            gen = Generator.create(
                meta["family"], meta["L"], **meta["all_params"])
            t1 = time.monotonic()
            say(f"  [1/5] Generator.create ok  ({t1 - t0:.3f}s)")

            say("  [2/5] gen.char_poly …")
            cpoly_bv = gen.char_poly()
            cpoly_int = int(cpoly_bv._val)
            t2 = time.monotonic()
            say(f"  [2/5] gen.char_poly ok  ({t2 - t1:.3f}s)  "
                f"cpoly=0x{cpoly_int:x}  hw={bin(cpoly_int).count('1') + 1}")

            say("  [3/5] Combination(J=1) …")
            comb = Combination(J=1, Lmax=gen.L)
            comb.components[0].add_gen(gen)
            comb.reset()
            t3 = time.monotonic()
            say(f"  [3/5] Combination ok  ({t3 - t2:.3f}s)  k_g={comb.k_g}")

            say("  [4/5] PISCache(…) ctor …")
            cache = _cpp.PISCache(
                [comb[0]._cpp_gen], [[]], comb.k_g, comb.L,
            )
            t4 = time.monotonic()
            say(f"  [4/5] PISCache ctor ok  ({t4 - t3:.3f}s)")

            say("  [5/5] cache.compute_all() …")
            ecart = cache.compute_all()
            t5 = time.monotonic()
            gaps = [int(ecart[v]) for v in range(1, comb.L + 1)]
            say(f"  [5/5] compute_all ok  ({t5 - t4:.3f}s)  "
                f"se={sum(gaps)}  gaps={gaps}")
            say(f"  total: {t5 - t0:.3f}s")

        except KeyboardInterrupt:
            print(f"  ^C after {time.monotonic() - t0:.3f}s "
                  "(this is the repro — attach py-spy next time)")
            raise
        except Exception as exc:
            print(f"  FAILED in {time.monotonic() - t0:.3f}s: "
                  f"{type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
