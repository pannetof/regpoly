#!/usr/bin/env bash
# Activate the venv and start regpoly-web against the repo's SQLite DB.
# Usage: ./regpoly-web.sh [extra args forwarded to regpoly-web]

set -eu

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# shellcheck disable=SC1091
source "$HERE/.venv/bin/activate"

exec regpoly-web --db "$HERE/regpoly.db" "$@"
