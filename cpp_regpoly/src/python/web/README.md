# regpoly web

A FastAPI-based web UI for running the regpoly primitive-generator and
tempering searches against a SQLite database.

## Install

From the `cpp_regpoly` directory:

```bash
pip install -e '.[web]'
```

This installs FastAPI, uvicorn, aiosqlite, Jinja2, python-multipart, and
markdown in addition to the base regpoly package.

## Run

```bash
regpoly-web                      # default: http://127.0.0.1:8000
regpoly-web --port 9000
regpoly-web --db /tmp/my.db      # choose the SQLite file
regpoly-web --reload             # dev mode (auto-reload on code changes)
```

Environment variables (alternative to CLI flags):

- `REGPOLY_DB`           — SQLite database path (default `regpoly.db`)
- `REGPOLY_HOST`         — bind address (default `127.0.0.1`)
- `REGPOLY_PORT`         — port (default `8000`)
- `REGPOLY_POOL_SIZE`    — max concurrent search processes (default `4`)
- `REGPOLY_POLL_SECONDS` — SSE progress polling interval (default `0.5`)
- `REGPOLY_DOCS_DIR`     — override path to the generator docs markdown
- `REGPOLY_RELOAD`       — set to `1` to enable hot reload

## First-time data

To populate the database from an existing regpoly YAML tree:

```bash
curl -X POST http://127.0.0.1:8000/api/import/generators-dir \
     -H 'Content-Type: application/json' \
     -d '{"directory": "/path/to/yaml/generators"}'
```

## Architecture

- **Backend**: FastAPI + aiosqlite.
- **Frontend**: server-rendered Jinja2 templates with HTMX and Alpine.js;
  no Node.js build step.
- **Background searches**: `ProcessPoolExecutor` — each search runs in
  its own OS process and writes results + progress to the shared SQLite
  database.
- **Live updates**: Server-Sent Events (SSE) streamed from
  `/api/*/progress` endpoints.

Files:

```
src/regpoly/web/
    app.py             FastAPI factory, entry point
    config.py          Runtime settings
    database.py        SQLite connection helpers
    schema.sql         Database schema
    models.py          Pydantic request/response models
    routes/            HTTP endpoints
    tasks/             Background search tasks
    templates/         Jinja2 HTML
    static/            CSS and JS
```

The existing regpoly library classes (`PrimitiveSearch`, `TemperingSearch`,
`TemperingOptimizer`, etc.) are not modified by the web app — the worker
tasks replicate the core search loops with SQLite writes in place of the
YAML file output.
