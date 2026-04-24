"""Primitive-search worker task.

Runs inside a child process (ProcessPoolExecutor).  Opens its own SQLite
connection, fetches the job description from primitive_search_run, runs
the search loop, and writes results + progress back to the database.
"""

from __future__ import annotations

import json
import time

import regpoly._regpoly_cpp as _cpp
from regpoly.generateur import Generateur, resolve_family
from regpoly.parametric import build_gen_enumerator
from regpoly.web.database import json_dumps, json_loads, sync_connect


_CANCEL_POLL_EVERY = 200   # check status for cancellation every N tries
_PROGRESS_EVERY = 100      # write a progress row at least every N tries
_PROGRESS_MIN_SECONDS = 2.0  # ...or every N seconds, whichever comes first


def run_primitive_search(db_path: str, run_id: int) -> None:
    """Entry point for a primitive search worker process."""
    with sync_connect(db_path) as conn:
        job = _fetch_job(conn, run_id)
        if job is None:
            return

        _mark_running(conn, run_id)

        try:
            if (job.get("search_mode") or "random") == "exhaustive":
                _run_exhaustive(conn, run_id, job)
            else:
                _run_random(conn, run_id, job)
        except Exception as exc:
            _mark_failed(conn, run_id, str(exc))


# ── Random mode (existing behaviour) ─────────────────────────────────────

def _run_random(conn, run_id: int, job: dict) -> None:
    tries = int(job["tries_done"] or 0)
    found = int(job["found_count"] or 0)
    resume = tries > 0

    family_raw = job["family"]
    family = resolve_family(family_raw, job["structural_params"])
    L = job["L"]
    structural = job["structural_params"]
    fixed_user = job["fixed_params"]
    max_tries = job["max_tries"]
    max_seconds = job["max_seconds"]

    try:
        specs = _cpp.get_gen_param_specs(family)
    except Exception as exc:
        _mark_failed(conn, run_id, f"Unknown family: {family} ({exc})")
        return

    fixed = dict(structural)
    for key, val in fixed_user.items():
        if val is not None:
            fixed[key] = val

    t_start = time.time()
    prior_elapsed = float(job["elapsed_seconds"] or 0.0)
    t_last_progress = t_start
    last_reported_found = found
    _write_progress(
        conn, run_id, tries, found,
        message="resuming" if resume else "starting",
        current_info={"phase": "resuming" if resume else "starting"},
    )

    while True:
        if max_tries is not None and tries >= max_tries:
            break
        if max_seconds is not None and (time.time() - t_start) >= max_seconds:
            break

        if tries and tries % _CANCEL_POLL_EVERY == 0:
            status = _read_status(conn, run_id)
            elapsed = prior_elapsed + (time.time() - t_start)
            if status == "cancelled":
                _mark_cancelled(conn, run_id, tries, found, elapsed)
                return
            if status == "paused":
                _update_run_stats(conn, run_id, tries, found, elapsed)
                return

        tries += 1

        try:
            gen = Generateur.create(family_raw, L, **fixed)
        except Exception:
            continue
        params = gen.params

        if gen.is_full_period():
            found += _insert_found(conn, run_id, family, L, gen,
                                    params, structural, tries)

        now = time.time()
        due_by_tries = tries % _PROGRESS_EVERY == 0
        due_by_time = (now - t_last_progress) >= _PROGRESS_MIN_SECONDS
        due_by_found = found != last_reported_found
        if due_by_tries or due_by_time or due_by_found:
            elapsed = prior_elapsed + (now - t_start)
            _update_run_stats(conn, run_id, tries, found, elapsed)
            _write_progress(
                conn, run_id, tries, found,
                current_info={
                    "k": _safe_k(gen),
                    "rate": tries / max(elapsed, 1e-9),
                },
            )
            t_last_progress = now
            last_reported_found = found

    elapsed = prior_elapsed + (time.time() - t_start)
    _mark_completed(conn, run_id, tries, found, elapsed)
    _write_progress(conn, run_id, tries, found,
                    message="completed",
                    current_info={"rate": tries / max(elapsed, 1e-9)})


# ── Exhaustive mode ──────────────────────────────────────────────────────

def _run_exhaustive(conn, run_id: int, job: dict) -> None:
    family_raw = job["family"]
    family = resolve_family(family_raw, job["structural_params"])
    L = job["L"]
    structural = job["structural_params"]
    fixed_user = job["fixed_params"]
    max_tries = job["max_tries"]
    max_seconds = job["max_seconds"]

    # `poly` is the enumerated output; exclude any lingering value.
    resolved = dict(structural)
    for key, val in fixed_user.items():
        if val is not None and key != "poly":
            resolved[key] = val

    enumerator = build_gen_enumerator(family_raw, L, resolved)
    if enumerator is None:
        _mark_failed(conn, run_id,
                     f"Family {family_raw!r} has no exhaustive enumerator.")
        return
    total = enumerator.total

    tries = int(job["tries_done"] or 0)
    found = int(job["found_count"] or 0)
    idx = int(job["enum_index"] or 0)
    resume = idx > 0

    t_start = time.time()
    prior_elapsed = float(job["elapsed_seconds"] or 0.0)
    t_last_progress = t_start
    last_reported_found = found
    _write_progress(
        conn, run_id, tries, found,
        message="resuming" if resume else "starting",
        current_info={
            "phase": "resuming" if resume else "starting",
            "enum_index": idx,
            "enum_total": str(total),
        },
    )

    while idx < total:
        if max_tries is not None and tries >= max_tries:
            break
        if max_seconds is not None and (time.time() - t_start) >= max_seconds:
            break

        if tries and tries % _CANCEL_POLL_EVERY == 0:
            status = _read_status(conn, run_id)
            elapsed = prior_elapsed + (time.time() - t_start)
            if status == "cancelled":
                _mark_cancelled(conn, run_id, tries, found, elapsed,
                                enum_index=idx)
                return
            if status == "paused":
                _update_run_stats(conn, run_id, tries, found, elapsed,
                                  enum_index=idx)
                return

        try:
            gen = Generateur.create_at_index(
                family_raw, L, enumerator, idx, **resolved)
            params = gen.params
        except Exception:
            idx += 1
            tries += 1
            continue

        if gen.is_full_period():
            found += _insert_found(conn, run_id, family, L, gen,
                                    params, structural, tries)

        idx += 1
        tries += 1

        now = time.time()
        due_by_tries = tries % _PROGRESS_EVERY == 0
        due_by_time = (now - t_last_progress) >= _PROGRESS_MIN_SECONDS
        due_by_found = found != last_reported_found
        if due_by_tries or due_by_time or due_by_found:
            elapsed = prior_elapsed + (now - t_start)
            _update_run_stats(conn, run_id, tries, found, elapsed,
                              enum_index=idx)
            _write_progress(
                conn, run_id, tries, found,
                current_info={
                    "enum_index": idx,
                    "enum_total": str(total),
                    "rate": tries / max(elapsed, 1e-9),
                },
            )
            t_last_progress = now
            last_reported_found = found

    elapsed = prior_elapsed + (time.time() - t_start)
    _update_run_stats(conn, run_id, tries, found, elapsed, enum_index=idx)
    _mark_completed(conn, run_id, tries, found, elapsed)
    _write_progress(conn, run_id, tries, found,
                    message="completed",
                    current_info={"enum_index": idx,
                                  "enum_total": str(total),
                                  "rate": tries / max(elapsed, 1e-9)})


def _insert_found(conn, run_id: int, family: str, L: int, gen,
                   params: dict, structural: dict, tries: int) -> int:
    """Insert a found full-period generator.  Returns 1 on fresh insert,
    0 on duplicate."""
    search_only = {
        n: _py_int(v) for n, v in params.items()
        if n not in structural
    }
    all_params = {**structural, **search_only}
    try:
        conn.execute(
            """
            INSERT OR IGNORE INTO primitive_generator
                (search_run_id, family, L, k,
                 structural_params, search_params,
                 all_params, found_at_try)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                run_id, family, L, gen.k,
                json_dumps(structural),
                json_dumps(search_only),
                json_dumps(all_params),
                tries,
            ),
        )
        inserted = 1 if conn.total_changes else 0
        conn.commit()
        return inserted
    except Exception:
        conn.rollback()
        return 0


# ── Helpers ──────────────────────────────────────────────────────────────

def _fetch_job(conn, run_id: int) -> dict | None:
    row = conn.execute(
        "SELECT * FROM primitive_search_run WHERE id = ?",
        (run_id,),
    ).fetchone()
    if row is None:
        return None
    return {
        "family": row["family"],
        "L": row["L"],
        "k": row["k"],
        "structural_params": json_loads(row["structural_params"]) or {},
        "fixed_params": json_loads(row["fixed_params"]) or {},
        "max_tries": row["max_tries"],
        "max_seconds": row["max_seconds"],
        "tries_done": row["tries_done"],
        "found_count": row["found_count"],
        "elapsed_seconds": row["elapsed_seconds"],
        "search_mode": _row_get(row, "search_mode", "random"),
        "enum_index":  _row_get(row, "enum_index", 0) or 0,
        "enum_total":  _row_get(row, "enum_total"),
    }


def _row_get(row, name, default=None):
    try:
        val = row[name]
    except (KeyError, IndexError):
        return default
    return val if val is not None else default


def _read_status(conn, run_id: int) -> str | None:
    row = conn.execute(
        "SELECT status FROM primitive_search_run WHERE id = ?",
        (run_id,),
    ).fetchone()
    return row["status"] if row else None


def _mark_running(conn, run_id: int) -> None:
    conn.execute(
        """
        UPDATE primitive_search_run
        SET status='running', started_at=datetime('now')
        WHERE id = ?
        """,
        (run_id,),
    )
    conn.commit()


def _mark_cancelled(conn, run_id: int, tries: int, found: int,
                    elapsed: float, enum_index: int | None = None) -> None:
    if enum_index is None:
        conn.execute(
            """
            UPDATE primitive_search_run
            SET status='cancelled', tries_done=?, found_count=?,
                elapsed_seconds=?, finished_at=datetime('now')
            WHERE id = ?
            """,
            (tries, found, elapsed, run_id),
        )
    else:
        conn.execute(
            """
            UPDATE primitive_search_run
            SET status='cancelled', tries_done=?, found_count=?,
                elapsed_seconds=?, enum_index=?,
                finished_at=datetime('now')
            WHERE id = ?
            """,
            (tries, found, elapsed, enum_index, run_id),
        )
    conn.commit()


def _mark_completed(conn, run_id: int, tries: int, found: int,
                    elapsed: float) -> None:
    conn.execute(
        """
        UPDATE primitive_search_run
        SET status='completed', tries_done=?, found_count=?,
            elapsed_seconds=?, finished_at=datetime('now')
        WHERE id = ?
        """,
        (tries, found, elapsed, run_id),
    )
    conn.commit()


def _mark_failed(conn, run_id: int, message: str) -> None:
    conn.execute(
        """
        UPDATE primitive_search_run
        SET status='failed', error_message=?, finished_at=datetime('now')
        WHERE id = ?
        """,
        (message, run_id),
    )
    conn.commit()


def _update_run_stats(conn, run_id: int, tries: int, found: int,
                      elapsed: float,
                      enum_index: int | None = None) -> None:
    if enum_index is None:
        conn.execute(
            """
            UPDATE primitive_search_run
            SET tries_done=?, found_count=?, elapsed_seconds=?
            WHERE id = ?
            """,
            (tries, found, elapsed, run_id),
        )
    else:
        conn.execute(
            """
            UPDATE primitive_search_run
            SET tries_done=?, found_count=?, elapsed_seconds=?,
                enum_index=?
            WHERE id = ?
            """,
            (tries, found, elapsed, enum_index, run_id),
        )
    conn.commit()


def _write_progress(conn, run_id: int, tries: int, found: int,
                    current_info: dict | None = None,
                    message: str | None = None) -> None:
    conn.execute(
        """
        INSERT INTO search_progress
            (search_type, search_run_id, tries_done, found_count,
             current_info, message)
        VALUES ('primitive', ?, ?, ?, ?, ?)
        """,
        (
            run_id, tries, found,
            json_dumps(current_info) if current_info else None,
            message,
        ),
    )
    # Keep only the ~10 most recent rows for this search
    conn.execute(
        """
        DELETE FROM search_progress
        WHERE search_type='primitive' AND search_run_id=?
          AND id NOT IN (
            SELECT id FROM search_progress
            WHERE search_type='primitive' AND search_run_id=?
            ORDER BY id DESC LIMIT 10
          )
        """,
        (run_id, run_id),
    )
    conn.commit()


def _py_int(v):
    if isinstance(v, int) and not isinstance(v, bool):
        return int(v)
    if isinstance(v, list):
        return [_py_int(x) for x in v]
    return v


def _safe_k(gen) -> int | None:
    try:
        return int(gen.k)
    except Exception:
        return None
