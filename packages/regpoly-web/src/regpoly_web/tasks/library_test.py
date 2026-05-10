# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Background worker for "Run test" actions on library generators.

The web request handler INSERTs a ``library_test_run`` row with
``status='pending'`` and returns its id immediately; the worker
scheduler picks it up via ``FOR UPDATE SKIP LOCKED`` and dispatches
this function in a child process. The actual computation never
blocks the HTTP request thread.
"""

from __future__ import annotations

import json
import time

from regpoly_web.database import json_dumps, json_loads, sync_connect
from regpoly_web.results import find_existing_typed_result, save_typed_result


def run_library_test(db_url: str, run_id: int) -> None:
    """Top-level worker function (must be picklable for the pool).

    Reads ``library_id`` + ``test_config`` from the
    ``library_test_run`` row; on success writes
    ``status='completed'`` + ``result`` (JSONB); on failure writes
    ``status='failed'`` + ``error_code``/``error_message``.
    """
    job = _fetch_job(db_url, run_id)
    if job is None:
        return  # row deleted between dispatch and start; nothing to do

    library_id = job["library_id"]
    test_config = job["test_config"]

    try:
        result = _execute(db_url, library_id, test_config)
        _mark_completed(db_url, run_id, result)
    except Exception as exc:
        _mark_failed(db_url, run_id, "test_failed", str(exc))


def _execute(db_url: str, library_id: str, test_config: dict) -> dict:
    from regpoly.library import Catalog

    from regpoly_web.config import find_library_dir

    library_dir = find_library_dir()
    if library_dir is None:
        raise RuntimeError("Could not locate docs/library directory")

    cat = Catalog(library_dir)
    cat.load()
    loc = cat.generator(library_id)
    if loc is None:
        raise RuntimeError(f"Unknown generator id: {library_id}")
    _, g = loc
    if not g.valid:
        raise RuntimeError(f"Library entry has validation errors: {g.errors}")

    tested_id, created = _instantiate_sync(db_url, g)

    test_type = test_config.get("type", "equidistribution")
    method = test_config.get("method")
    config_json = json_dumps(test_config)

    with sync_connect(db_url) as conn:
        prior = find_existing_typed_result(conn, tested_id, test_type, config_json)
        if prior is not None:
            return {
                "library_id": g.id,
                "tested_generator_id": int(tested_id),
                "test_result_id": int(prior["id"]),
                "se": int(prior["se"]) if prior["se"] is not None else None,
                "is_me": bool(prior["is_me"]),
                "instantiated": created,
                "test_type": test_type,
                "method": method,
            }

    from regpoly.core.combination import Combination
    from regpoly.core.combination_build import build_combinaison_inputs

    from regpoly_web.tasks._test_build import build_test
    from regpoly_web.tasks.tempering import _is_me

    gen_lists, temperings = build_combinaison_inputs(g.components, g.Lmax)
    comb = Combination.CreateFromFiles(gen_lists, g.Lmax, temperings)
    comb.reset()
    test = build_test(test_config, g.Lmax)
    t0 = time.time()
    result = test.run(comb)
    elapsed = time.time() - t0

    se = getattr(result, "se", None)
    is_me = bool(_is_me(result)) if se is not None else False

    with sync_connect(db_url) as conn:
        rid = save_typed_result(
            conn, tested_id, test_config, result,
            kg=comb.k_g, L=comb.Lmax, elapsed_seconds=elapsed,
        )
        conn.commit()
        return {
            "library_id": g.id,
            "tested_generator_id": int(tested_id),
            "test_result_id": int(rid) if rid is not None else None,
            "se": int(se) if se is not None else None,
            "is_me": is_me,
            "instantiated": created,
            "test_type": test_type,
            "method": method,
        }


# ── library_test_run row helpers ───────────────────────────────────

def _fetch_job(db_url: str, run_id: int) -> dict | None:
    with sync_connect(db_url) as conn:
        conn.autocommit = True
        row = conn.execute(
            "SELECT library_id, test_type, method, test_config "
            "FROM library_test_run WHERE id = ?",
            (run_id,),
        ).fetchone()
    if row is None:
        return None
    return {
        "library_id": row["library_id"],
        "test_type": row["test_type"],
        "method": row["method"],
        "test_config": json_loads(row["test_config"]) or {},
    }


def _mark_completed(db_url: str, run_id: int, result: dict) -> None:
    with sync_connect(db_url) as conn:
        conn.execute(
            """
            UPDATE library_test_run
               SET status='completed',
                   result=?::jsonb,
                   updated_at=NOW(),
                   finished_at=NOW()
             WHERE id = ?
            """,
            (json.dumps(result), run_id),
        )
        conn.commit()


def _mark_failed(db_url: str, run_id: int, code: str, message: str) -> None:
    with sync_connect(db_url) as conn:
        conn.execute(
            """
            UPDATE library_test_run
               SET status='failed',
                   error_code=?, error_message=?,
                   updated_at=NOW(),
                   finished_at=NOW()
             WHERE id = ?
            """,
            (code, message, run_id),
        )
        conn.commit()


def _instantiate_sync(db_path, g) -> tuple[int, bool]:
    """Mirror of routes/library._instantiate_sync, duplicated here so the
    worker doesn't import the FastAPI route module (which pulls in
    APIRouter/Request and is unnecessary in a child process)."""
    from regpoly.core.combination import Combination
    from regpoly.core.combination_build import build_combinaison_inputs

    from regpoly_web.tasks.tempering import _find_primitive_id

    with sync_connect(db_path) as conn:
        row = conn.execute(
            "SELECT id FROM tested_generator WHERE library_id = ?",
            (g.id,),
        ).fetchone()
        if row is not None:
            return int(row["id"]), False

        gen_lists, temperings = build_combinaison_inputs(g.components, g.Lmax)
        comb = Combination.CreateFromFiles(gen_lists, g.Lmax, temperings)
        comb.reset()

        cur = conn.execute(
            "INSERT INTO tested_generator"
            "(search_run_id, Lmax, k_g, J, library_id) "
            "VALUES (?, ?, ?, ?, ?)",
            (None, comb.Lmax, comb.k_g, comb.J, g.id),
        )
        tested_id = cur.lastrowid

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


__all__ = ["run_library_test"]
