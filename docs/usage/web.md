# Using the REGPOLY web UI

`regpoly-web` is a FastAPI server backed by a SQLite result store. It exposes a browser UI for kicking off full-period and tempering searches, watching their progress over Server-Sent Events, browsing tested generators, and inspecting equidistribution results.

## Start the server

```bash
uv run regpoly-web --db packages/regpoly-web/var/regpoly.db
```

| Flag | Default | Notes |
|---|---|---|
| `--host HOST` | `127.0.0.1` | Bind address. |
| `--port PORT` | `8000` | TCP port. |
| `--db DB` | `var/regpoly.db` | SQLite file. Auto-created on first start; v1 databases are migrated to v2 in place at startup (see [Database schema](#database-schema)). |
| `--reload` | off | Hot-reload Python sources (development only). |

Then point a browser at <http://localhost:8000/>.

## Pages

- `/` — family directory + published-configurations catalog.
- `/searches` — list of every full-period and combined-generator search, current and past. Pause / resume / cancel / restart from here.
- `/primitive-search` — start a new full-period search.
- `/tempering-search` — start a new combined-generator (tempering) search.
- `/tested-generators` — every tested generator on file, filterable by family, k_g, test type, and equidistribution class.
- `/tested-generators/{id}` — one tested generator with its components, tempering chain, and per-test results.
- `/library/{paper_id}` — papers in the catalog and the published parameter sets they contain.

## Database schema

The web app keeps state in a SQLite file. Schema v2 (Phase 5.4) introduces three typed result tables that mirror the C++ analysis structs 1:1 — `equidistribution_result`, `collision_free_result`, `tuplets_result` — so the UI receives the full typed result without going through opaque JSON.

Old v1 databases are auto-upgraded the first time the server opens them: `init_sync` runs `regpoly_web.migrations.v2.migrate_v2_inplace`, which backfills the typed tables from the legacy `test_result.detail` blobs. Migration is idempotent.

To run the migration offline (e.g. before deploying a new server build):

```bash
uv run python packages/regpoly-web/scripts/migrate_v2.py path/to/regpoly.db
```

The legacy `test_result` table is preserved for backward-compat reads through Phase 5; Phase 8 drops it.

## Headless smoke and end-to-end tests

The default test lane already exercises the FastAPI app via `TestClient`:

```bash
uv run pytest packages/regpoly-web/tests/
```

The 8 Playwright golden-path tests are opt-in (gated by `@pytest.mark.e2e`) so they don't run on every CI build. To run them locally:

```bash
uv sync --group test --group e2e
uv run playwright install chromium
uv run pytest -m e2e packages/regpoly-web/tests/e2e/
```

## See also

- [Python usage](python.md) — same algorithms without the web layer.
- [C++ usage](cpp.md) — the standalone `regpoly-cli`.
- [Architecture](../dev/architecture.md) — how the web app talks to the C++ core through `regpoly`.
