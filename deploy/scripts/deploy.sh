#!/usr/bin/env bash
# deploy.sh — pull the latest GHCR image and (re)launch the stack
# (dockerize-plan Phase 7).
#
# Run as the `regpoly` user on the CAX21, with deploy/.env populated
# and the working directory at the repo root.
#
# Drains the worker first so the signal handler cancels any in-flight
# rows cleanly (status=cancelled, error_message=shutdown). Without
# this, init's --reap-orphans would clobber the user's just-submitted
# searches on every redeploy.

set -euo pipefail

cd "$(dirname "$0")/../.."  # repo root

if [ ! -f deploy/.env ]; then
    echo "error: deploy/.env missing — copy from .env.example and fill in values" >&2
    exit 1
fi

# shellcheck disable=SC1091
source deploy/.env

COMPOSE_FILE=deploy/compose.yml
DC="docker compose -f $COMPOSE_FILE"

echo "== drain worker (timeout 70s for in-flight searches)"
if $DC ps -q worker 2>/dev/null | grep -q .; then
    $DC stop --timeout 70 worker
else
    echo "  (worker not running; nothing to drain)"
fi

echo "== pull latest images"
# Public GHCR image — no docker login needed.
$DC pull

echo "== bring up the stack"
$DC up -d --remove-orphans

echo "== wait for healthchecks (~120s)"
for i in $(seq 1 24); do
    pending=$($DC ps --format json 2>/dev/null \
        | grep -E '"Service":"(web|worker|db)"' \
        | grep -cv '"Health":"healthy"' || true)
    if [ "${pending:-0}" -eq 0 ]; then
        echo "  ✓ all services healthy"
        break
    fi
    echo "  ($pending services not yet healthy; sleeping 5s)"
    sleep 5
done

echo "== wait for cert acquisition (~180s budget)"
for i in $(seq 1 36); do
    if $DC logs caddy 2>&1 \
        | grep -qE "certificate obtained successfully|certificate is up to date"; then
        echo "  ✓ Let's Encrypt cert ready"
        break
    fi
    sleep 5
done

echo "== HTTPS smoke test"
if curl -fsS -o /dev/null "https://${DOMAIN}/healthz"; then
    echo "  ✓ /healthz reachable over HTTPS"
else
    echo "  ✗ /healthz failed — check 'docker compose logs caddy web'"
fi

echo "== prune unused images"
docker image prune -f >/dev/null

echo
echo "✓ deploy complete. Stack:"
$DC ps

echo
echo "Hint: if you updated WELL data and need the DB-payload migration:"
echo "  $DC exec worker python /app/packages/regpoly-web/scripts/migrate_well_matrices.py /var/regpoly/regpoly.db"
echo "  (opt-in per CLAUDE.md auto-memory; not run automatically)"
