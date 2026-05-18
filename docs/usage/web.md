# Using the REGPOLY web UI

`regpoly-web` is a FastAPI server backed by a PostgreSQL result store (via `psycopg3`). It exposes a browser UI for kicking off full-period and tempering searches, watching their progress over Server-Sent Events, browsing tested generators, and inspecting equidistribution results.

## Start the server

```bash
uv run regpoly-web --db-url postgresql://user:pw@host/regpoly
# or via the environment:
REGPOLY_DB_URL=postgresql://user:pw@host/regpoly uv run regpoly-web
```

| Flag | Default | Notes |
|---|---|---|
| `--host HOST` | `127.0.0.1` | Bind address (env: `REGPOLY_HOST`). |
| `--port PORT` | `8000` | TCP port (env: `REGPOLY_PORT`). |
| `--db-url URL` | *(env: `REGPOLY_DB_URL`)* | PostgreSQL DSN. Schema is applied at lifespan startup from `schema.sql`. |
| `--reload` | off | Hot-reload Python sources (development only). |

For local development without a system Postgres, the `pgserver` package (installed via the `test` dependency group) can spin up an ephemeral instance. See [`docs/dev/postgres.md`](../dev/postgres.md).

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

The web app keeps state in PostgreSQL. The schema (`packages/regpoly-web/src/regpoly_web/schema.sql`) is applied at FastAPI lifespan startup via `psycopg3`. Three typed result tables mirror the C++ analysis structs 1:1 — `equidistribution_result`, `collision_free_result`, `tuplets_result` — so the UI receives the full typed result without opaque JSON.

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
