# Contributing

> **Status**: Phase 0 placeholder.

Code lives in three packages under `packages/`. Tests live alongside each package under `tests/`. The active redesign roadmap is at `~/.claude/plans/i-want-to-have-rippling-starlight.md` (local, not committed).

When adding a new generator family:

1. Add the C++ implementation under `packages/regpoly-cpp/src/generators/` and register it in `factory.cpp`.
2. Add a GoogleTest case under `packages/regpoly-cpp/tests/`.
3. Add a doc page under `docs/generators/<family>.md` from the template (Phase 6+).
4. Add a notebook under `docs/notebooks/families/<family>.ipynb` from the template (Phase 7+).
5. If the family has known parameter sets from the literature, add a catalog entry under `docs/library/`.

CI is gated by:

- `ctest` (C++ unit tests).
- `pytest` default lane (no `slow` / `e2e`).
- `ruff check`, `ruff format --check`.
- `import-linter` (enforces `regpoly-web → regpoly → regpoly_cpp`).

The slow lane runs nightly and on `slow`-labeled PRs.
