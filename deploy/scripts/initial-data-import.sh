#!/usr/bin/env bash
# initial-data-import.sh — one-shot SQLite → PG import for migrating
# existing regpoly.db rows on first deploy (dockerize-plan Phase 7).
#
# Operator-run: AFTER bootstrap-cax21.sh + deploy.sh have brought up
# the empty PG-backed stack. SCPs your local SQLite DB to the CAX21
# and runs the migrate_sqlite_to_pg.py tool inside the web container.
#
# Usage (from your dev machine, NOT the CAX21):
#   bash deploy/scripts/initial-data-import.sh \
#       /path/to/local/regpoly.db                     \
#       regpoly@regpoly.frpanneton.ca                  \
#       [--dry-run]
#
# Exit 0 on success.

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "usage: $0 <local-sqlite-path> <ssh-target> [--dry-run]" >&2
    exit 1
fi

LOCAL_DB="$1"
HOST="$2"
EXTRA="${3:-}"

if [ ! -f "$LOCAL_DB" ]; then
    echo "error: source SQLite file not found: $LOCAL_DB" >&2
    exit 2
fi

REMOTE_PATH="/tmp/regpoly-import-$(date +%s).db"

echo "== scp $LOCAL_DB → $HOST:$REMOTE_PATH"
scp "$LOCAL_DB" "${HOST}:${REMOTE_PATH}"

echo "== copy DB into the web container"
ssh "$HOST" "docker compose -f /home/regpoly/regpoly/deploy/compose.yml \
    cp ${REMOTE_PATH} web:/tmp/regpoly-import.db"

echo "== run migrate_sqlite_to_pg.py inside web ${EXTRA}"
# shellcheck disable=SC2029
ssh "$HOST" "docker compose -f /home/regpoly/regpoly/deploy/compose.yml \
    exec -T web python /app/packages/regpoly-web/scripts/migrate_sqlite_to_pg.py \
    /tmp/regpoly-import.db ${EXTRA}"

echo "== cleanup remote tmp"
# shellcheck disable=SC2029
ssh "$HOST" "rm -f ${REMOTE_PATH}"

echo
echo "✓ import complete. Verify row counts via:"
echo "  ssh $HOST 'docker compose -f /home/regpoly/regpoly/deploy/compose.yml exec db \\"
echo "      psql -U regpoly -d regpoly -c \"SELECT count(*) FROM primitive_search_run;\"'"
