# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**regpoly** is the canonical REGPOLY engine — a research/academic toolkit for analysing and searching for high-quality combined pseudo-random number generators (PRNGs) based on modulo-2 linear recurrences (LFSRs over GF(2)). It is a **single distribution** built by scikit-build-core, runnable two ways:

- **From Python / Jupyter** — `import regpoly`; the native core is the `_regpoly_cpp` pybind11 extension (importable as both `regpoly._regpoly_cpp` and `regpoly_cpp._regpoly_cpp`).
- **Standalone from C++** — `cmake -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF` builds `libregpoly_core.a` + the `regpoly-cli` binary + `find_package(regpoly)` config, with **no pybind11/Python** needed.

This repo was carved out of `regpoly_monorepo` (kept as a frozen archive). The FastAPI web app + Postgres/docker deployment live in a **separate** repo (`pannetof/regpoly-web`), which depends on this one as a published `regpoly` distribution. Nothing here imports the web layer.

## Layout

**Clean C++ / Python demarcation:** all C++ is under `cpp/`, all Python under
`python/`, and the shared catalog data is the neutral top-level `data/`.

```
regpoly/
├── pyproject.toml          scikit-build-core (cmake.source-dir="cpp", wheel.packages=python/*)
├── cpp/                    ← ALL C++, a self-contained CMake project (see cpp/README.md)
│   ├── CMakeLists.txt      pybind11 wheel mode OR pure-C++ via -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF
│   ├── CMakePresets.json   `default` (+ extension) / `cpp-only` (pure C++) presets
│   ├── include/regpoly/    public umbrella header (find_package consumers)
│   ├── cmake/              regpolyConfig.cmake.in + FindNTL.cmake
│   ├── src/                {algebra,analyses,core,generators,transforms,lattice,library,
│   │                         search,yaml_config}/ + bindings/ (pybind11) + cli/ + include/
│   ├── tests/              GoogleTest C++ tests (ctest) + fixtures/
│   └── examples/           standalone find_package(regpoly) consumer (CI acceptance gate)
├── python/                 ← ALL Python
│   ├── regpoly/            wrapper: core/ analyses/ search/ io/ data/ library/ tools/
│   └── regpoly_cpp/        thin package; the built _regpoly_cpp.so lands here (gitignored)
├── data/                   ← neutral catalog data (read by C++ CLI + Python)
│   ├── library/            catalog (paper YAML)
│   └── papers/             reference PDFs (catalog `pdf:` field is `papers/<file>`)
├── tests/                  pytest (wrapper layer) + tests/bindings/ (pybind11 ABI)
├── shared/yaml/            example search configs + generator pools (not in the wheel)
└── docs/                   Sphinx site (theory, generators, notebooks, C++/Python API)
```

The catalog/papers ship in the wheel as `regpoly/_data/` (CMake installs them there
when `REGPOLY_CATALOG_WHEEL_LAYOUT=ON`); a pure-C++ install puts them at
`share/regpoly/{catalog,papers}`. The Python resolver
(`regpoly.library.bundled_library_dir`) uses the packaged `_data` when installed,
falling back to top-level `data/` in an editable/dev tree.

## Build & test

```bash
uv sync                                  # editable install + builds the extension (yaml-cpp via FetchContent)
uv run regpoly shared/yaml/equidist/mt19937.yaml

# C++ tests (Python-wheel mode):
cmake -S cpp -B build -DREGPOLY_BUILD_TESTS=ON && cmake --build build -j
ctest --test-dir build --output-on-failure

# Pure-C++ mode (no pybind11) — or use: cmake --preset cpp-only
cmake -S cpp -B build-cpp -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF -DREGPOLY_BUILD_TESTS=ON
cmake --build build-cpp -j && ./build-cpp/regpoly-cli catalog list
# install + consume from another project: see cpp/README.md + cpp/examples/

uv run pytest                            # default lane (skips slow + e2e)
uv run pytest -m slow                    # MTToolBox cross-checks, nbmake notebooks
```

## Catalog resolution

The catalog is the neutral top-level `data/` tree. Resolution is two-tier — explicit override then bundled data:

- **Python**: `Catalog()` (no arg) → `regpoly.library.bundled_library_dir()`: the packaged `regpoly/_data` when installed from a wheel, else the top-level `data/` (editable/dev). Also `bundled_papers_dir()`.
- **C++ CLI**: `--library DIR` / `-l` → `$REGPOLY_CATALOG_DIR` → the installed `share/regpoly/catalog` (resolved exe-relative) → the build-time `REGPOLY_SOURCE_CATALOG` define (points at top-level `data/library`, when present).
- Paper `pdf:` paths are resolved as a sibling of the catalog dir (`.../papers/<file>`), so they work in the source tree, the wheel, and an installed prefix alike.

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
