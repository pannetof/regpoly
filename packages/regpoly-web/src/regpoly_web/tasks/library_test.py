# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Background worker for "Run test" actions on library generators.

The web request handler returns a job_id immediately; the actual test
is computed inside a ProcessPoolExecutor worker so the HTTP response
never blocks on the C++ equidistribution computation.

The worker is fully self-contained: it re-loads the YAML library
catalog from disk (so the parent's in-memory catalog does not need to
be pickled across processes), looks up the generator by id, then runs
the same instantiate-then-test logic that previously lived inline in
``routes/library.py``.
"""

from __future__ import annotations

import time

from regpoly_web.database import json_dumps, sync_connect
from regpoly_web.results import find_existing_typed_result, save_typed_result


def run_library_test(
    db_path: str,
    library_id: str,
    test_config: dict,
) -> dict:
    """Top-level worker function (must be picklable for the pool).

    Returns a result dict with keys ``library_id``,
    ``tested_generator_id``, ``test_result_id``, ``se``, ``is_me``,
    ``instantiated``.  Raises on hard failures so the caller sees the
    original exception via ``Future.exception()``.
    """
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

    tested_id, created = _instantiate_sync(db_path, g)

    test_type = test_config.get("type", "equidistribution")
    config_json = json_dumps(test_config)

    with sync_connect(db_path) as conn:
        prior = find_existing_typed_result(conn, tested_id, test_type, config_json)
        if prior is not None:
            return {
                "library_id": g.id,
                "tested_generator_id": int(tested_id),
                "test_result_id": int(prior["id"]),
                "se": int(prior["se"]) if prior["se"] is not None else None,
                "is_me": bool(prior["is_me"]),
                "instantiated": created,
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

    with sync_connect(db_path) as conn:
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
        }


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
