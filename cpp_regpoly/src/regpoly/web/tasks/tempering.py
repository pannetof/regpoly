"""Tempering-search worker task.

Runs inside a child process.  Loads the search configuration from the DB
(components, tempering configs, generator pools, test config, optimizer
config), builds the Combinaison, and runs the search loop — mirroring
TemperingSearch.run() but writing results and progress to SQLite
instead of YAML.
"""

from __future__ import annotations

import json
import time
from typing import Any

from regpoly.combinaison import Combinaison
from regpoly.generateur import Generateur
from regpoly.transformation import Transformation
from regpoly.web.database import json_dumps, json_loads, sync_connect


_CANCEL_POLL_EVERY = 1   # check cancellation once per combo


def run_tempering_search(db_path: str, run_id: int) -> None:
    with sync_connect(db_path) as conn:
        job = _fetch_job(conn, run_id)
        if job is None:
            return

        try:
            _mark_running(conn, run_id)
            components = _load_components(conn, run_id)
            test = _build_test(job["test_config"], job["Lmax"])
            gen_lists, temperings = _build_gen_and_trans(components, job["Lmax"])
        except Exception as exc:
            _mark_failed(conn, run_id, f"setup failed: {exc}")
            return

        try:
            comb = Combinaison.CreateFromFiles(
                gen_lists, job["Lmax"], temperings
            )
            if not comb.reset():
                pool_sizes = [len(pool) for pool in gen_lists]
                msg = (
                    "No valid generator combination — check that each "
                    "component has at least one generator and that pools "
                    "share at most one generator at a time.  Pool sizes: "
                    f"{pool_sizes}"
                )
                _mark_failed(conn, run_id, msg)
                return
        except Exception as exc:
            _mark_failed(conn, run_id, f"no valid combination: {exc}")
            return

        # Fast-forward the iterator to skip combinations already processed
        # in a previous run segment (for resume after pause/crash).
        already_done = int(job["combos_done"] or 0)
        for _ in range(already_done):
            try:
                next(comb)
            except StopIteration:
                _mark_completed(
                    conn, run_id, already_done, float(job["elapsed_seconds"] or 0.0),
                    best_se=job["best_se"],
                )
                return

        optimizer = _build_optimizer(job["optimizer_config"])
        use_optimizer = optimizer is not None and _has_optimizable(comb)

        t_start = time.time()
        prior_elapsed = float(job["elapsed_seconds"] or 0.0)
        combo_idx = already_done
        best_overall_se: int | None = job["best_se"]

        try:
            while True:
                combo_idx += 1

                status = _read_status(conn, run_id)
                if status == "cancelled":
                    _mark_cancelled(
                        conn, run_id, combo_idx - 1,
                        prior_elapsed + (time.time() - t_start),
                        best_overall_se,
                    )
                    return
                if status == "paused":
                    # Preserve combos_done and best_se for resume
                    conn.execute(
                        "UPDATE tempering_search_run SET combos_done=?, "
                        "best_se=?, elapsed_seconds=? WHERE id = ?",
                        (combo_idx - 1, best_overall_se,
                         prior_elapsed + (time.time() - t_start), run_id),
                    )
                    conn.commit()
                    return

                result = _search_one_combo(
                    conn, run_id, comb, test, optimizer, use_optimizer,
                    job["nb_tries"], combo_idx, t_start, best_overall_se,
                )
                if result is not None:
                    se = result["se"]
                    if best_overall_se is None or se < best_overall_se:
                        best_overall_se = se
                    _save_tested_result(
                        conn, run_id, comb, result, job["test_config"]
                    )

                conn.execute(
                    "UPDATE tempering_search_run SET combos_done=?, "
                    "best_se=?, elapsed_seconds=? WHERE id = ?",
                    (combo_idx, best_overall_se,
                     prior_elapsed + (time.time() - t_start), run_id),
                )
                conn.commit()

                try:
                    next(comb)
                except StopIteration:
                    break
        except Exception as exc:
            _mark_failed(conn, run_id, str(exc))
            return

        _mark_completed(
            conn, run_id, combo_idx,
            prior_elapsed + (time.time() - t_start),
            best_se=best_overall_se,
        )


# ── Combination loop ─────────────────────────────────────────────────────

def _search_one_combo(conn, run_id: int, comb: Combinaison, test,
                      optimizer, use_optimizer: bool,
                      nb_tries: int, combo_idx: int,
                      t_total_start: float,
                      best_overall_se: int | None) -> dict | None:
    best_se: int | None = None
    best_result = None
    best_params: list[list[dict]] | None = None

    for t in range(nb_tries):
        if t % 10 == 0 and _read_status(conn, run_id) in ("cancelled", "paused"):
            return None

        for comp in comb.components:
            for trans in comp.trans:
                trans.randomize_params()

        if use_optimizer:
            try:
                optimizer.run(comb)
            except ValueError:
                pass

        test_result = test.run(comb)
        se = _score(test_result)

        if best_se is None or se < best_se:
            best_se = se
            best_result = test_result
            best_params = _save_params(comb)

            _write_progress(
                conn, run_id,
                current_info={
                    "combo_idx": combo_idx,
                    "try": t + 1,
                    "nb_tries": nb_tries,
                    "se": se,
                    "best_overall_se": best_overall_se
                        if best_overall_se is None or se >= best_overall_se
                        else se,
                },
                message=f"combo {combo_idx} try {t + 1}: se={se}",
            )

            if se == 0:
                break

    if best_params is None:
        return None

    # Honour the test's acceptance criterion — for equidistribution
    # that's is_presque_me() (all per-v gaps within delta[v] AND
    # se <= mse).  If the best try across nb_tries still fails the
    # criterion, do not save the combination.  This is what makes
    # mse=0 act as a hard filter.
    if not _result_passes(best_result):
        return None

    _restore_params(comb, best_params)

    return {
        "se": best_se,
        "result": best_result,
        "params": best_params,
    }


def _result_passes(result) -> bool:
    """True iff the result meets the test's canonical acceptance
    predicate.  For equidistribution this is ``is_presque_me``
    (gaps + mse); for collision-free, ``is_cf``; for tuplets,
    ``is_ok``.  A result with none of these predicates is accepted
    unconditionally."""
    for name in ("is_presque_me", "is_ok", "is_cf"):
        fn = getattr(result, name, None)
        if callable(fn):
            try:
                return bool(fn())
            except Exception:
                continue
    return True


# ── Saving a single tested-generator record ──────────────────────────────

def _save_tested_result(conn, run_id: int, comb: Combinaison,
                         outcome: dict, test_config: dict) -> None:
    best_result = outcome["result"]

    cur = conn.execute(
        """
        INSERT INTO tested_generator(search_run_id, Lmax, k_g, J)
        VALUES (?, ?, ?, ?)
        """,
        (run_id, comb.Lmax, comb.k_g, comb.J),
    )
    tested_id = cur.lastrowid

    # Look up which primitive_generator each active gen came from
    for j in range(comb.J):
        gen = comb[j]
        comp = comb.components[j]
        gen_id = _find_primitive_id(conn, gen)
        tempering = [
            {"type": t._type_name, **t._params}
            for t in comp.trans
        ]
        conn.execute(
            """
            INSERT INTO tested_generator_component
                (tested_gen_id, component_index, generator_id,
                 family, L, k, all_params, tempering_params)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                tested_id, j, gen_id, gen.type_name, gen.L, gen.k,
                json_dumps(gen.params),
                json_dumps(tempering),
            ),
        )

    detail = _build_result_detail(best_result)
    is_me = 1 if _is_me(best_result) else 0
    conn.execute(
        """
        INSERT INTO test_result
            (tested_gen_id, test_type, test_config,
             se, is_me, secf, is_cf, score, detail)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            tested_id,
            test_config.get("type", "unknown"),
            json_dumps(test_config),
            outcome["se"],
            is_me,
            None,
            None,
            float(outcome["se"]),
            json_dumps(detail),
        ),
    )
    conn.commit()


def _find_primitive_id(conn, gen: Generateur) -> int | None:
    """Try to locate the primitive_generator row that corresponds to gen."""
    structural = gen.structural_params()
    search = {
        k: v for k, v in gen.params.items() if k not in structural
    }
    row = conn.execute(
        """
        SELECT id FROM primitive_generator
        WHERE family = ? AND structural_params = ? AND search_params = ?
        """,
        (gen.type_name, json_dumps(structural), json_dumps(search)),
    ).fetchone()
    return row["id"] if row else None


def _build_result_detail(result) -> dict:
    detail: dict[str, Any] = {}
    if hasattr(result, "se"):
        detail["se"] = result.se
        if hasattr(result, "ecart") and hasattr(result, "L"):
            # `ecart[l]` is only meaningful when psi12[l] is True — other
            # entries carry a sentinel (≈2^63-1 / INT_MAX) from uncomputed
            # resolutions.  Filter them out so they never reach the UI.
            psi12 = getattr(result, "psi12", None)
            SENTINEL_MAX = 1 << 40   # any gap this large is not real
            ecart = {}
            for l in range(1, result.L + 1):
                if psi12 is not None and not psi12[l]:
                    continue
                v = int(result.ecart[l])
                if v == 0 or v >= SENTINEL_MAX or v < 0:
                    continue
                ecart[str(l)] = v
            if ecart:
                detail["ecart"] = ecart
    return detail


def _is_me(result) -> bool:
    return bool(getattr(result, "is_me", lambda: False)())


# ── Setup helpers ────────────────────────────────────────────────────────

def _fetch_job(conn, run_id: int) -> dict | None:
    row = conn.execute(
        "SELECT * FROM tempering_search_run WHERE id = ?", (run_id,)
    ).fetchone()
    if row is None:
        return None
    return {
        "test_type": row["test_type"],
        "test_config": json_loads(row["test_config"]),
        "Lmax": row["Lmax"],
        "nb_tries": row["nb_tries"],
        "optimizer_config": json_loads(row["optimizer_config"]),
        "combos_done": row["combos_done"],
        "best_se": row["best_se"],
        "elapsed_seconds": row["elapsed_seconds"],
    }


def _load_components(conn, run_id: int) -> list[dict]:
    rows = conn.execute(
        "SELECT * FROM tempering_search_component WHERE search_run_id = ? "
        "ORDER BY component_index",
        (run_id,),
    ).fetchall()
    components = []
    for row in rows:
        gen_rows = conn.execute(
            "SELECT pg.* FROM tempering_search_generator tsg "
            "JOIN primitive_generator pg ON pg.id = tsg.generator_id "
            "WHERE tsg.component_id = ? ORDER BY pg.id",
            (row["id"],),
        ).fetchall()
        components.append({
            "component_index": row["component_index"],
            "shared_with_component": row["shared_with_component"],
            "tempering_config": json_loads(row["tempering_config"]),
            "generators": [
                {
                    "id": g["id"],
                    "family": g["family"],
                    "L": g["L"],
                    "k": g["k"],
                    "all_params": json_loads(g["all_params"]),
                }
                for g in gen_rows
            ],
        })
    return components


def _build_gen_and_trans(components: list[dict], Lmax: int):
    """Build per-component lists of Generateur and Transformation instances."""
    gen_pools: list[list[Generateur]] = []
    pool_by_index: dict[int, list[Generateur]] = {}

    for c in components:
        if c["shared_with_component"] is not None \
                and c["shared_with_component"] in pool_by_index:
            pool = pool_by_index[c["shared_with_component"]]
        else:
            pool = []
            for g in c["generators"]:
                pool.append(
                    # Always rebuild each component at Lmax so every
                    # component of the combined generator emits the
                    # same output width as the one the test analyses.
                    # Using the stored per-generator L (e.g. 64 for
                    # Tausworthe) while the combined Lmax is 32 would
                    # produce a different sequence than the one the
                    # paper / test assumes.
                    Generateur.create(g["family"], Lmax, **g["all_params"])
                )
        gen_pools.append(pool)
        pool_by_index[c["component_index"]] = pool

    temperings: list[list[Transformation]] = []
    for c in components:
        chain = []
        for tconf in c["tempering_config"]:
            ttype = tconf["type"]
            tparams = {k: v for k, v in tconf.items() if k != "type"}
            chain.append(Transformation.create(ttype, **tparams))
        temperings.append(chain)

    return gen_pools, temperings


def _build_test(test_config: dict, Lmax: int):
    from regpoly.analyses.equidistribution_test import EquidistributionTest
    from regpoly.analyses.collision_free_test import CollisionFreeTest
    from regpoly.analyses.tuplets_test import TupletsTest

    dispatch = {
        "equidistribution": EquidistributionTest,
        "collision_free": CollisionFreeTest,
        "tuplets": TupletsTest,
    }
    ttype = test_config.get("type", "equidistribution")
    cls = dispatch.get(ttype)
    if cls is None:
        raise ValueError(f"Unknown test type: {ttype}")
    return cls._from_params(test_config, Lmax)


def _build_optimizer(cfg: dict | None):
    if cfg is None:
        return None
    from regpoly.tempering_optimizer import TemperingOptimizer
    return TemperingOptimizer(
        max_essais=cfg.get("max_essais", 400),
        delta=cfg.get("delta"),
        mse=cfg.get("mse"),
        n_restarts=cfg.get("n_restarts", 1),
        verbose=False,
    )


def _has_optimizable(comb) -> bool:
    for comp in comb.components:
        for trans in comp.trans:
            if trans.optimizable_params():
                return True
    return False


def _score(result) -> int:
    if hasattr(result, "se"):
        return int(result.se)
    if hasattr(result, "verified") and result.verified:
        return 0
    return 1


def _save_params(comb):
    return [[dict(trans.params) for trans in comp.trans]
            for comp in comb.components]


def _restore_params(comb, saved):
    for comp, comp_params in zip(comb.components, saved):
        for trans, params in zip(comp.trans, comp_params):
            for name, value in params.items():
                if value != trans.get_param(name):
                    trans.set_param(name, value)


# ── Status updates ───────────────────────────────────────────────────────

def _read_status(conn, run_id: int) -> str | None:
    row = conn.execute(
        "SELECT status FROM tempering_search_run WHERE id = ?", (run_id,)
    ).fetchone()
    return row["status"] if row else None


def _mark_running(conn, run_id: int) -> None:
    conn.execute(
        "UPDATE tempering_search_run SET status='running', "
        "started_at=datetime('now') WHERE id = ?",
        (run_id,),
    )
    conn.commit()


def _mark_cancelled(conn, run_id: int, combos_done: int, elapsed: float,
                    best_se: int | None) -> None:
    conn.execute(
        """
        UPDATE tempering_search_run
        SET status='cancelled', combos_done=?, elapsed_seconds=?,
            best_se=?, finished_at=datetime('now')
        WHERE id = ?
        """,
        (combos_done, elapsed, best_se, run_id),
    )
    conn.commit()


def _mark_completed(conn, run_id: int, combos_done: int, elapsed: float,
                    best_se: int | None) -> None:
    conn.execute(
        """
        UPDATE tempering_search_run
        SET status='completed', combos_done=?, elapsed_seconds=?,
            best_se=?, finished_at=datetime('now')
        WHERE id = ?
        """,
        (combos_done, elapsed, best_se, run_id),
    )
    conn.commit()


def _mark_failed(conn, run_id: int, message: str) -> None:
    conn.execute(
        "UPDATE tempering_search_run SET status='failed', error_message=?, "
        "finished_at=datetime('now') WHERE id = ?",
        (message, run_id),
    )
    conn.commit()


def _write_progress(conn, run_id: int, current_info: dict | None = None,
                    message: str | None = None) -> None:
    conn.execute(
        """
        INSERT INTO search_progress
            (search_type, search_run_id, tries_done, found_count,
             current_info, message)
        VALUES ('tempering', ?, 0, 0, ?, ?)
        """,
        (
            run_id,
            json_dumps(current_info) if current_info else None,
            message,
        ),
    )
    conn.execute(
        """
        DELETE FROM search_progress
        WHERE search_type='tempering' AND search_run_id=?
          AND id NOT IN (
            SELECT id FROM search_progress
            WHERE search_type='tempering' AND search_run_id=?
            ORDER BY id DESC LIMIT 10
          )
        """,
        (run_id, run_id),
    )
    conn.commit()
