# deploy/

Build and runtime artifacts for the regpoly Docker stack.

**The deployment guide lives at the repo root: [`../DEPLOY.md`](../DEPLOY.md).** Setup, configuration, access patterns, ops, backups, and rationale are all there. This file is just a directory index.

## Files

- `Dockerfile` — multi-stage build for `ghcr.io/pannetof/regpoly`. Builds portably (`-O3`, no host-CPU tuning) so the same Dockerfile produces both `linux/amd64` and `linux/arm64` images.
- `compose.yml` — three services: `db` (Postgres), `web` (FastAPI), `worker` (schedulers). Single bridge network. Named `pgdata` volume.
- `.env.example` — env template. Copy to `.env`, fill in `DB_PASSWORD`, `chmod 600 .env`.
- `scripts/initial-data-import.sh` — operator-run, one-time migration of a legacy SQLite database into the PG stack. Not part of the routine deploy path.
