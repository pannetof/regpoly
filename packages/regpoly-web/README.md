# regpoly-web

FastAPI-based web UI for running the regpoly primitive-generator,
tempering, and equidistribution searches. Persists results to
**PostgreSQL** via `psycopg3`.

## Install

The package is a workspace member of the `regpoly_monorepo`. The
canonical install path is `uv sync` from the workspace root, which
builds every package and installs them editable:

```bash
cd regpoly_monorepo
uv sync
```

This pulls in FastAPI, uvicorn, `psycopg[binary,pool]`, Jinja2,
`python-multipart`, and `markdown` on top of `regpoly` and
`regpoly-cpp`.

## Run

```bash
regpoly-web                                            # defaults; reads $REGPOLY_DB_URL
regpoly-web --port 9000
regpoly-web --db-url postgresql://user:pw@host/regpoly
regpoly-web --reload                                   # dev mode (auto-reload on edits)
```

For local development without a system Postgres, the `pgserver` package
(installed via the `test` dependency group) can spin up an ephemeral
Postgres instance. See [docs/dev/postgres.md](../../docs/dev/postgres.md).

## Environment variables (alternative to CLI flags)

| Var | Purpose | Default |
|---|---|---|
| `REGPOLY_DB_URL`               | PostgreSQL DSN — required for production runs | *(unset → app errors)* |
| `REGPOLY_HOST`                 | bind address                                  | `127.0.0.1` |
| `REGPOLY_PORT`                 | port                                          | `8000` |
| `REGPOLY_POOL_SIZE`            | max concurrent search processes               | `4` |
| `REGPOLY_ANALYSIS_POOL_SIZE`   | analysis worker pool                          | derived from cpu count |
| `REGPOLY_POLL_SECONDS`         | SSE progress polling interval                 | `0.5` |
| `REGPOLY_DB_POOL_MIN`          | psycopg pool min size                         | `2` |
| `REGPOLY_DB_POOL_MAX`          | psycopg pool max size                         | `10` |
| `REGPOLY_DOCS_DIR`             | override path to generator docs markdown      | repo-relative |
| `REGPOLY_LIBRARY_DIR`          | override path to catalog YAMLs                | repo-relative |
| `REGPOLY_PAPERS_DIR`           | override path to reference PDFs               | repo-relative |
| `REGPOLY_RELOAD`               | set to `1` to enable hot reload               | `0` |
| `REGPOLY_WORKER`               | set to `1` when the process is a worker       | `0` |

## Architecture

- **Backend**: FastAPI + `psycopg3` (sync + async from one library).
- **Database**: PostgreSQL 16; schema lives in `schema.sql` and is
  applied at app startup via the FastAPI lifespan.
- **Frontend**: server-rendered Jinja2 templates with HTMX and Alpine.js
  — no Node.js build step.
- **Background searches**: `ProcessPoolExecutor` — each search runs in
  its own OS process and writes results + progress to the shared
  PostgreSQL database. A scheduler component dispatches queued jobs;
  see `tasks/`.
- **Live updates**: Server-Sent Events (SSE) streamed from
  `/api/*/progress` endpoints.

## Files

```
src/regpoly_web/
├── app.py              FastAPI factory + entry point (also defines main())
├── config.py           Runtime settings; reads REGPOLY_* env vars
├── database.py         psycopg3 connection helpers + pool management
├── schema.sql          Database schema (applied at lifespan startup)
├── models.py           Pydantic request/response models
├── routes/             HTTP endpoints (one module per logical area)
├── tasks/              Background search workers + scheduler
├── templates/          Jinja2 HTML
└── static/             CSS and JS
```

The existing `regpoly` library classes (`PrimitiveSearch`,
`TemperingSearch`, `TemperingOptimizer`, etc.) are not modified by the
web app — the worker tasks invoke the core search loops with database
writes in place of the YAML file output.
