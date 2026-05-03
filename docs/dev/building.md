# Building

```bash
# Workspace-wide install (builds the C++ extension):
uv sync

# C++ tests (after Phase 0):
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build

# Python tests (per package):
uv run pytest packages/regpoly/tests
uv run pytest packages/regpoly-cpp/tests/python
uv run pytest packages/regpoly-web/tests

# Slow lane (MTToolBox cross-checks, nbmake notebooks, Playwright e2e):
uv run pytest -m slow
uv run pytest -m e2e

# Docs site:
uv sync --group docs
cd docs && uv run mkdocs serve
```
