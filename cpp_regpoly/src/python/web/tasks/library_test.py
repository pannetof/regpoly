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

from regpoly.web.database import json_dumps, sync_connect


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
    from regpoly.web.config import find_library_dir

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
        row = conn.execute(
            "SELECT id, se, is_me FROM test_result "
            "WHERE tested_gen_id = ? AND test_type = ? AND test_config = ? "
            "ORDER BY id DESC LIMIT 1",
            (tested_id, test_type, config_json),
        ).fetchone()
        if row is not None:
            return {
                "library_id": g.id,
                "tested_generator_id": int(tested_id),
                "test_result_id": int(row["id"]),
                "se": int(row["se"]) if row["se"] is not None else None,
                "is_me": bool(row["is_me"]),
                "instantiated": created,
            }

    from regpoly.core.combination import Combination
    from regpoly.core.combination_build import build_combinaison_inputs
    from regpoly.web.tasks._test_build import build_test
    from regpoly.web.tasks.tempering import _build_result_detail, _is_me

    gen_lists, temperings = build_combinaison_inputs(g.components, g.Lmax)
    comb = Combination.CreateFromFiles(gen_lists, g.Lmax, temperings)
    comb.reset()
    test = build_test(test_config, g.Lmax)
    result = test.run(comb)

    se = getattr(result, "se", None)
    is_me = bool(_is_me(result)) if se is not None else False
    detail = _build_result_detail(result)

    with sync_connect(db_path) as conn:
        cur = conn.execute(
            "INSERT INTO test_result"
            "(tested_gen_id, test_type, test_config, "
            " se, is_me, secf, is_cf, score, detail) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (tested_id, test_type, config_json,
             se, 1 if is_me else 0,
             None, None,
             float(se) if se is not None else None,
             json_dumps(detail)),
        )
        conn.commit()
        return {
            "library_id": g.id,
            "tested_generator_id": int(tested_id),
            "test_result_id": int(cur.lastrowid),
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
    from regpoly.web.tasks.tempering import _find_primitive_id

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
