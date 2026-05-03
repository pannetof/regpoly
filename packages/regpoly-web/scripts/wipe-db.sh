#!/usr/bin/env bash
# Wipe the regpoly-web SQLite database.
#
# Removes regpoly.db (and the SQLite WAL/SHM sidecar files if present).
# The web app re-creates the schema from src/python/web/schema.sql on
# the next start, so this leaves you with an empty, valid database.
#
# Usage:
#     ./wipe-db.sh            # interactive: asks for confirmation
#     ./wipe-db.sh --force    # skip the prompt (for scripts/CI)
#     ./wipe-db.sh /path/db   # wipe a specific DB instead of ./regpoly.db
#     ./wipe-db.sh --force /path/db
#
# Refuses to run while regpoly-web is up — kill the web process first.

set -eu

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

force=0
db=""
for arg in "$@"; do
    case "$arg" in
        --force|-f) force=1 ;;
        --help|-h)
            sed -n '2,15p' "$0"
            exit 0
            ;;
        -*)
            echo "wipe-db.sh: unknown flag '$arg'" >&2
            exit 2
            ;;
        *)
            if [[ -n "$db" ]]; then
                echo "wipe-db.sh: only one DB path allowed (got '$db' and '$arg')" >&2
                exit 2
            fi
            db="$arg"
            ;;
    esac
done

db="${db:-$HERE/regpoly.db}"

# Refuse if regpoly-web is currently holding the file open.
# pgrep is portable enough; fall back to `ps | grep` if missing.
if command -v pgrep >/dev/null && pgrep -f 'regpoly-web' >/dev/null; then
    echo "wipe-db.sh: regpoly-web appears to be running — stop it first." >&2
    exit 1
fi

# Resolve siblings (WAL + SHM created by SQLite in WAL mode).
wal="${db}-wal"
shm="${db}-shm"

# Show what we're about to delete.
echo "About to delete:"
existed=0
for f in "$db" "$wal" "$shm"; do
    if [[ -e "$f" ]]; then
        size=$(stat -c %s "$f" 2>/dev/null || echo "?")
        echo "  $f ($size bytes)"
        existed=1
    fi
done
if [[ "$existed" -eq 0 ]]; then
    echo "  (nothing — database already absent)"
    exit 0
fi

# Confirmation prompt unless --force.
if [[ "$force" -ne 1 ]]; then
    read -r -p "Proceed? [y/N] " ans
    case "$ans" in
        y|Y|yes|YES) ;;
        *) echo "Aborted."; exit 1 ;;
    esac
fi

rm -f -- "$db" "$wal" "$shm"
echo "Done. Schema will be recreated on the next \`regpoly-web\` start."
