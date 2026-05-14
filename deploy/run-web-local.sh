#!/usr/bin/env bash
# Run regpoly-web on the host against the stack's docker Postgres.
# Loads deploy/.env, builds the host-side DSN, then exec's uvicorn.
set -euo pipefail
cd "$(dirname "$0")/.."

set -a
# shellcheck disable=SC1091
source deploy/.env
set +a

export REGPOLY_DB_URL="postgresql://regpoly:${DB_PASSWORD}@127.0.0.1:5432/regpoly"
exec uv run regpoly-web "$@"
