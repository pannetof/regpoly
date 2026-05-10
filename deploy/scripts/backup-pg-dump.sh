#!/bin/sh
# backup-pg-dump.sh — emit a custom-format pg_dump on stdout
# (dockerize-plan Phase 7).
#
# Lives at /home/regpoly/bin/backup-pg-dump.sh on the CAX21 and is the
# only command the NAS-side SSH key is allowed to invoke (via
# `command="..."` in authorized_keys). Wrapping `docker compose exec`
# in a stable script path means renaming the deploy directory doesn't
# break the backup pipeline — only this script's COMPOSE_FILE env
# needs to track the location.
#
# The NAS-side script captures stdout, gzip-9s it, and atomically
# renames into the daily backup file. Restore is `pg_restore` of the
# uncompressed dump.

set -euo pipefail
COMPOSE_FILE="${COMPOSE_FILE:-/home/regpoly/regpoly/deploy/compose.yml}"
exec docker compose -f "$COMPOSE_FILE" exec -T db \
    pg_dump -U regpoly --no-password --format=custom regpoly
