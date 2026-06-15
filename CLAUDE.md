# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**regpoly** is the canonical REGPOLY engine — a research/academic toolkit for analysing and searching for high-quality combined pseudo-random number generators (PRNGs) based on modulo-2 linear recurrences (LFSRs over GF(2)). It is a **single distribution** built by scikit-build-core, runnable two ways:

- **From Python / Jupyter** — `import regpoly`; the native core is the `_regpoly_cpp` pybind11 extension (importable as both `regpoly._regpoly_cpp` and `regpoly_cpp._regpoly_cpp`).
- **Standalone from C++** — `cmake -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF` builds `libregpoly_core.a` + the `regpoly-cli` binary + `find_package(regpoly)` config, with **no pybind11/Python** needed.

This repo was carved out of `regpoly_monorepo` (kept as a frozen archive). The FastAPI web app + Postgres/docker deployment live in a **separate** repo (`pannetof/regpoly-web`), which depends on this one as a published `regpoly` distribution. Nothing here imports the web layer.

## Layout

```
regpoly/
├── pyproject.toml          scikit-build-core build + ruff/pytest config (single package)
├── CMakeLists.txt          C++ build (pybind11 wheel mode OR pure-C++ via -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF)
├── src/
│   ├── regpoly/            pure-Python wrapper: core/ analyses/ search/ io/ data/ library/ tools/
│   │   └── _data/          PACKAGE DATA — runtime source of truth, shipped in wheel + sdist:
│   │       ├── library/    catalog (paper YAML) read by the C++ + Python Catalog
│   │       └── papers/     reference PDFs (the catalog `pdf:` field is `papers/<file>`)
│   ├── regpoly_cpp/        thin Python package; the _regpoly_cpp extension lands here
│   ├── {algebra,analyses,core,generators,transforms,lattice,library,search,yaml_config}/  C++ sources
│   ├── bindings/           pybind11 registration
│   ├── cli/                regpoly-cli (main.cpp)
│   └── include/            C++ headers
├── include/regpoly/        public C++ API headers (find_package consumers)
├── cmake/                  regpolyConfig.cmake.in
├── cpp-tests/              GoogleTest C++ tests (ctest) + python/ (pybind11 ABI tests)
├── tests/                  pytest for the Python wrapper layer
├── shared/yaml/            example search configs + generator pools (not shipped in the wheel)
└── docs/                   Sphinx site (theory, generators, notebooks, C++/Python API)
```

## Build & test

```bash
uv sync                                  # editable install + builds the extension (yaml-cpp via FetchContent)
uv run regpoly shared/yaml/equidist/mt19937.yaml

# C++ tests (Python-wheel mode):
cmake -S . -B build -DREGPOLY_BUILD_TESTS=ON && cmake --build build -j
ctest --test-dir build --output-on-failure

# Pure-C++ mode (no pybind11):
cmake -S . -B build-cpp -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF -DREGPOLY_BUILD_TESTS=ON
cmake --build build-cpp -j && ./build-cpp/regpoly-cli catalog list

uv run pytest                            # default lane (skips slow + e2e)
uv run pytest -m slow                    # MTToolBox cross-checks, nbmake notebooks
```

## Catalog resolution

The catalog lives in package data (`src/regpoly/_data/library`). Resolution is two-tier — explicit override then bundled data:

- **Python**: `Catalog()` (no arg) → `regpoly.library.bundled_library_dir()` (`importlib.resources`). Also `bundled_papers_dir()`.
- **C++ CLI**: `--library DIR` / `-l` → `$REGPOLY_CATALOG_DIR` → the build-time `REGPOLY_SOURCE_CATALOG` define (points at the in-tree `_data/library`). `cmake --install` copies the catalog/papers to `share/regpoly/{catalog,papers}`.
- Paper `pdf:` paths are resolved as a sibling of the catalog dir (`<_data>/papers/<file>`), so they work in the source tree, the wheel, and an installed prefix alike.

## Conventions

- **One distribution.** Do not re-split `regpoly_cpp` into a separate PyPI package. It is a sub-package of `regpoly`; the layering (web → regpoly → regpoly_cpp) is enforced in the *App* repo's import-linter, not here.
- **No git submodules.** `yaml-cpp` is pinned `FetchContent` (set `FETCHCONTENT_SOURCE_DIR_YAML-CPP` for offline builds). The MTToolBox cross-check clones on demand in the slow CI lane.
- **Publishing.** `v*` tag → `.github/workflows/publish.yml` builds an **sdist** and publishes to PyPI (needs `PYPI_API_TOKEN`). sdist-only by design (no binary wheels yet); consumers compile (needs NTL/GMP/gf2x + cmake). **The bundled PDFs are copyrighted — run a license check before the first public publish.**
- **Math in docs** renders as proper math (`$...$`), never backticked code. Mirror `docs/theory/notation.md`.

## Known issue

`tests/test_mttoolbox_crosscheck.py::test_regpoly_matches_mttoolbox_d5` was silently dormant for the monorepo's whole life (a path bug); the split activated it and revealed `mt-mt19937` is slow (~150s) and disagrees with the frozen reference. It is `@pytest.mark.slow` (default lane unaffected). Re-baseline the `*_mttoolbox_d5.py` references before relying on it.

## Specs (read when relevant)

- `docs/theory/equidistribution-spec.md` — design of record for matricial equidistribution on non-full-period F₂-linear generators.
- `docs/theory/antithetic-check.md` — local-antitheticity test for a linear RNG point set.
