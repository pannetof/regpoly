# Building

```bash
# Workspace-wide install (builds the C++ extension + every workspace member):
uv sync

# C++ tests via the default (Python wheel) build:
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build

# Pure-C++ build, no Python / pybind11 needed at all
# (produces libregpoly_core.a + regpoly-cli + headers + find_package(regpoly) config):
cmake -S packages/regpoly-cpp -B build-cpp \
      -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF \
      -DREGPOLY_BUILD_TESTS=ON
cmake --build build-cpp
ctest --test-dir build-cpp
sudo cmake --install build-cpp        # optional: makes find_package(regpoly) available system-wide

# Python tests (per package):
uv run pytest packages/regpoly/tests
uv run pytest packages/regpoly-cpp/tests/python
uv run pytest packages/regpoly-web/tests
uv run pytest packages/regpoly-legacy/tests/python

# Slow lane (MTToolBox cross-checks, nbmake notebooks, Playwright e2e):
uv run pytest -m slow
uv run pytest -m e2e

# Docs site (Sphinx + Doxygen + Breathe + Exhale + PyData theme):
uv sync --group docs
sudo apt install -y doxygen graphviz                      # one-time
uv run --group docs sphinx-build -b html docs site        # one-shot build
uv run --group docs sphinx-autobuild docs site            # live-reload preview
```
