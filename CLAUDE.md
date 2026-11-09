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

**Resolved (2026-11-08):** `tests/test_mttoolbox_crosscheck.py::test_regpoly_matches_mttoolbox_d5[mt-mt19937]`
disagreed with the frozen MTToolBox reference. The frozen references were
*correct* — no re-baseline was needed. Two independent bugs:

1. **Every leaf `copy()` dropped the inherited tempering chain.** 17 of 19
   generator families implemented `copy()` by constructing a fresh object and
   hand-copying their own fields, silently omitting `F2LinearSource::tempering_`
   (only `sobol`/`niederreiter_f2` used the prescribed
   `make_unique<X>(*this)` form). Since `GaussMatrix::prepare` and `find_polys`
   call `clone_source()` per matrix row, **every** equidistribution analysis ran
   on *untempered* output regardless of the configured tempering. MT19937
   reported `SE = 41650` instead of 6750. Fixed by switching all leaf `copy()`
   overrides to the copy constructor, which also immunises the whole class of
   bug against future field additions. Verified: all 32 `DIM`/`ECART` values now
   match Table II of Matsumoto–Nishimura 1998 exactly, `SE = 6750`.
2. **The harness drove MT through the wrong method.** `MersenneTwister` was
   absent from `_METHOD_FOR_FAMILY`, so it defaulted to `METHOD_NOTPRIMITIVE`;
   it is full-period, so `METHOD_DUALLATTICE` is correct (and is what
   `shared/yaml/equidist/mt19937.yaml` uses). The harness also passed
   `EquidistributionTest(L=MAX_V)`, conflating comparison granularity with the
   generator's output word width (that `L` drives the column stride `j*L + cl`);
   the matricial path rejects such a mismatch via `self.L < gen_L`, but the
   lattice/notprimitive paths return before that guard.

**Still open — `recover_char_poly` mis-recovers χ:** for MT19937 it returns a
**degree-1** polynomial instead of degree 19937, so `select_phi` yields `p = 1`,
the `if (p == kg) return test_me_lat(...)` fast path in
`cpp/src/lattice/me_notprimitive.cpp` never fires, and the generator falls into
the non-primitive branch reporting spurious **all-zero** gaps. `char_poly()`
itself is correct (degree 19937, Hamming weight 135 per the paper) and
`is_full_period()` returns `True`, so the defect is isolated to the BM-based
state-sequence recovery (`build_state_seq` / `step_once` round-trip). Routing MT
to `lattice` sidesteps it; any genuinely non-full-period generator still depends
on this path, so it needs a real fix.

Since then (2026-09-19), running it in CI surfaced two more issues, both now worked around rather than fixed:
- Cases above a catalog's `slow_threshold_mexp` (e.g. `dsfmt44497` and bigger, `mt19937`) run the `O(p^3/W)`-per-phase notprimitive/SIMD-notprimitive reduction (`docs/theory/equidistribution-spec.md` §6.7) for minutes-to-tens-of-minutes as a synchronous, GIL-holding pybind11 call — `pytest-timeout`'s `--timeout=120` can't preempt it mid-call, so instead of failing cleanly it hangs and then hard-kills the whole session once the call finally returns. Now `pytest.skip`'d by default; opt in locally with `REGPOLY_RUN_HUGE_CROSSCHECK=1` plus a generous `--timeout` override.
- `rmt607` (`data/library/rmt_params.yaml`) crashes the interpreter outright with `Fatal Python error: Aborted` (SIGABRT) inside `test_me_notprimitive` — almost certainly an NTL-internal `abort()` (NTL's `Error()` handler aborts by default rather than throwing) on some degenerate characteristic-polynomial factorization, not a catchable C++ exception. `rmt521` (same family/method, same code path) runs fine immediately before it, so this isn't a blanket family failure, but since a SIGABRT can't be `xfail`'d or caught, the whole `rmt_params.yaml` family is deferred (`deferred: true`) pending an actual debugger session with a C++ build environment.

## Specs (read when relevant)

- `docs/theory/equidistribution-spec.md` — design of record for matricial equidistribution on non-full-period F₂-linear generators.
- `docs/theory/antithetic-check.md` — local-antitheticity test for a linear RNG point set.
