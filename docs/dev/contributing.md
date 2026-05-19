# Contributing

Code lives in four packages under `packages/`. Tests live alongside each package under `tests/`.

When adding a new generator family:

1. Add the C++ implementation under `packages/regpoly-cpp/src/generators/` and register it in `factory.cpp`.
2. Add a GoogleTest case under `packages/regpoly-cpp/tests/`.
3. Add a doc page under `docs/generators/<family>.md` (copy [`docs/generators/_template.md`](../generators/_template.md) and fill it in).
4. Add a per-family notebook by extending `docs/notebooks/families/_stamp.py` (the parametric driver that stamps `<family>.ipynb`).
5. If the family has known parameter sets from the literature, add a catalog entry under `docs/library/`.

CI is gated by:

- `ctest` (C++ unit tests).
- `pytest` default lane (no `slow` / `e2e`).
- `ruff check`, `ruff format --check`.
- `import-linter` — two layered contracts:
    - `regpoly_web → regpoly → regpoly_cpp`
    - `regpoly_legacy → regpoly → regpoly_cpp`

The slow lane (`-m slow`) runs MTToolBox cross-checks and `nbmake` notebook execution; the e2e lane (`-m e2e`) runs Playwright golden-path browser tests for the web app.

## API documentation

Every new public symbol (Python `__all__` member, or header declaration
under `packages/regpoly-cpp/src/include/`) ships with a docstring or
Doxygen block. See the [API documentation guide](api-docs.md) for the
required tags, examples, and how to render the docs locally.
