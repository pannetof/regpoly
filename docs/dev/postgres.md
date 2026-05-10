# Local PostgreSQL for regpoly-web development

Since the dockerize-and-deploy work (commit `d200a94…`) the web app uses
PostgreSQL 16 instead of SQLite. Local development now needs a Postgres
instance.

## Quickstart (Docker)

```sh
docker run --rm -d \
    --name regpoly-pg \
    -p 5432:5432 \
    -e POSTGRES_USER=regpoly \
    -e POSTGRES_PASSWORD=devpw \
    -e POSTGRES_DB=regpoly \
    postgres:16.6-alpine

export REGPOLY_DB_URL=postgresql://regpoly:devpw@localhost:5432/regpoly
uv run regpoly-web
```

To stop and discard the dev DB:

```sh
docker stop regpoly-pg
```

`--rm` discards the container on stop; data does not persist between
restarts. For a persistent dev DB, drop `--rm` and add
`-v regpoly-pg-data:/var/lib/postgresql/data`.

## Apt-installed PG (alternative)

```sh
sudo apt-get install -y postgresql-16
sudo -u postgres psql -c "CREATE USER regpoly WITH PASSWORD 'devpw';"
sudo -u postgres psql -c "CREATE DATABASE regpoly OWNER regpoly;"
export REGPOLY_DB_URL=postgresql://regpoly:devpw@localhost:5432/regpoly
```

## Tests

The test fixtures use `pytest-postgresql` to spin up an ephemeral
Postgres per test session. It needs `pg_ctl` on `PATH`, or a configured
`postgresql_proc` fixture pointing at a running instance.

If you have the apt PG installed, tests work out-of-the-box. With the
Docker PG above, configure `pytest-postgresql` to reuse it via
`--postgresql-host=localhost --postgresql-port=5432
--postgresql-user=regpoly --postgresql-password=devpw` (see the
[plugin docs](https://github.com/ClearcodeHQ/pytest-postgresql)).

## Migrating an existing SQLite DB

If you have rows in `packages/regpoly-web/var/regpoly.db` you want to
preserve, run the one-shot migration script after the PG instance is up:

```sh
uv run python packages/regpoly-web/scripts/migrate_sqlite_to_pg.py \
    packages/regpoly-web/var/regpoly.db
```

The script is idempotent and reports per-table counts. Use `--dry-run`
to preview without writing.
