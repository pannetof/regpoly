#!/bin/bash
# Compare transition matrix output between C and Python programs.
# Usage: ./compare_transition.sh <generator_file>

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <generator_file>" >&2
    exit 1
fi

FILE="$1"
DIR="$(cd "$(dirname "$0")" && pwd)"

C_PROG="$DIR/test_transition_matrix_generator"

if [ ! -x "$C_PROG" ]; then
    echo "Error: C program not found at $C_PROG (run 'make test_transition_matrix_generator')" >&2
    exit 1
fi

C_OUT=$(mktemp)
PY_OUT=$(mktemp)
trap 'rm -f "$C_OUT" "$PY_OUT"' EXIT

"$C_PROG" "$FILE" > "$C_OUT" 2>&1
python3 -m regpoly.tools.transition_matrix "$FILE" > "$PY_OUT" 2>&1

if diff -q "$C_OUT" "$PY_OUT" > /dev/null 2>&1; then
    echo "OK: C and Python outputs are identical for $FILE"
else
    echo "MISMATCH: C and Python outputs differ for $FILE"
    diff "$C_OUT" "$PY_OUT"
    exit 1
fi
