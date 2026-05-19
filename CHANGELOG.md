# Changelog

## Unreleased — documentation framework migration: MkDocs → Sphinx

Replaced the MkDocs Material site with a Sphinx + Doxygen + Breathe +
Exhale stack using the PyData Sphinx Theme — the canonical "scientific
Python library" framework (NumPy / SciPy / pandas / scikit-learn /
ITK / VTK all use it). Migration delivers:

- **Hierarchical C++ API**: Exhale auto-generates Class Hierarchy /
  File Hierarchy / Namespace Hierarchy / Full API Index pages from
  Doxygen XML. Click `regpoly::core::Generator` → see the full
  inheritance tree of every `*Gen` family.
- **Single visual language** across Python, C++, narrative pages, and
  notebooks (was: mkdocstrings Material widgets for Python +
  raw-Doxygen tables for C++).
- **Bidirectional Python ↔ C++ cross-references**: every Python class
  with a C++ counterpart links it via the Sphinx `:cpp:class:` domain
  role; every C++ class with a Python wrapper carries `@see :py:class:`
  in its Doxygen block. Click either side, jump to the other.
- **Hash-based notebook execution** (MyST-NB + jupyter-cache, with the
  cache key including the compiled pybind11 `.so` so notebook outputs
  invalidate on C++ algorithm changes).

Documentation conventions ([style guide](docs/dev/api-docs.md)):

- **Python**: NumPy-style docstrings via Napoleon. Required sections:
  `Parameters`, `Returns`, `Raises`. Recommended: `See Also`,
  `Examples`, `Notes`. Rewrote 28 docstrings across 8 modules from
  Google style; added bidirectional `:cpp:class:` cross-refs.
- **C++**: full Doxygen blocks (`@brief @param @return @throws
  @ingroup`) on every public symbol declared under
  `packages/regpoly-cpp/src/include/` (55 headers, ~830 blocks).
  `@code{.cpp}` examples on every user-constructible class.
  `@return` singular (not `@returns`). Class-level long descriptions
  preserve invariant prose.

New user-facing on-ramp pages:

- [Python / C++ bridge](docs/usage/python-cpp-bridge.md) — explains
  the layered split + the wrapper-to-C++ cross-reference table.
- [Equidistribution in 5 minutes](docs/theory/equidistribution-in-5-minutes.md)
  — beginner on-ramp before the design-of-record spec.
- [C++ usage](docs/usage/cpp.md) rewrite with a complete consumer
  `CMakeLists.txt` example and a runnable "hello equidistribution"
  C++ snippet that compiles against the pure-C++ install.

Build commands changed:

```
# Was:
cd docs && uv run mkdocs serve

# Now:
sudo apt install -y doxygen graphviz   # one-time
uv run --group docs sphinx-build -b html docs site
# or, for live-reload:
uv run --group docs sphinx-autobuild docs site
```

CI (`.github/workflows/docs.yml`):

- Builds the full uv workspace (including the pybind11 extension) so
  autodoc has live signatures and `:cpp:class:` cross-refs resolve.
- `apt install doxygen graphviz` step added.
- Jupyter-cache action-cached, keyed on a hash of every `.ipynb` +
  the compiled `.so`.
- Deploys via `actions/deploy-pages` (was: `mkdocs gh-deploy`).
- `sphinx-build -W` (warnings-as-errors) deferred until the
  Exhale-duplicate-declaration baseline drops below the current 62.

Removed: `mkdocs.yml`, `docs/gen_paper_pages.py` (logic now in a
`builder-inited` hook in `docs/conf.py`), 26 obsolete
`docs/api/python/*.md` mkdocstrings stubs, `docs/api/index.md`,
`docs/api/cpp/index.md`, `docs/dev/api-docs-style.md` (replaced by
`docs/dev/api-docs.md`), `docs/usage/notebooks.md` (replaced by
`docs/notebooks/index.md`).

Workspace notebook reorganisation: moved the 14 workspace-root
reference notebooks (`notebooks/*.ipynb`) into
`docs/notebooks/reference/` via `git mv` (preserves history). New
`docs/notebooks/index.md` hub renders both the 15 family notebooks
and the 14 reference notebooks. The 2 CA research notebooks under
`docs/notebooks/research/` are excluded from execution.

Canonical exemplars when documenting a new symbol:

- C++: `packages/regpoly-cpp/src/include/library/catalog.h` (full
  convention) and `core/generator.h` (abstract base + SIMD-override
  prose preservation).
- Python: `packages/regpoly/src/regpoly/library/__init__.py`
  (NumPy docstrings, bidirectional `:cpp:class:` cross-refs).
- Style guide: [docs/dev/api-docs.md](docs/dev/api-docs.md).

## Unreleased — legacy `.dat` reader extracted to add-on

**Breaking changes** (C++ public API + YAML schema + CLI surface):

- The C++ public API no longer exposes `<regpoly/legacy_reader.h>` or
  the `regpoly_legacy::` namespace. Both `legacy_reader.{cpp,h}` and
  `well_legacy_decode.{cpp,h}` were deleted from `packages/regpoly-cpp/`.
- `regpoly-cli` no longer offers `legacy-info` / `legacy-trans` subcommands.
- The YAML search schema rejects the `legacy_file:` component-source key
  (both `seek_config.cpp` and Python `Seek.from_yaml` raise a clear error
  pointing at the add-on).
- The `regpoly` Python CLI (`uv run regpoly ...`) no longer accepts the
  positional `nb_comp test_file gen_file1 [...]` form. `Seek.from_legacy()`
  is gone.
- The pybind11 module `_regpoly_cpp` no longer exposes
  `legacy_read_generator_specs` / `legacy_read_transformation_specs`.
- `from regpoly import LegacyReader` raises `ImportError`.

**Replacement**: install the new optional pure-Python `regpoly-legacy`
add-on. It ships:

- `regpoly_legacy.LegacyReader` (drop-in API replacement for the old
  `regpoly.io.legacy_reader.LegacyReader`).
- `regpoly_legacy.seek_from_legacy(nb_comp, test_file, gen_data_files)`.
- A `regpoly-legacy` console script with three subcommands: `info`,
  `trans`, `seek`.

The add-on parses `.dat` files in pure Python (no compiler needed) and
constructs the resulting `Generator` / `Transformation` objects through
the existing `regpoly.Generator.create()` factory — same C++ core under
the hood, no behavioural drift.

Historical legacy fixtures move from `shared/legacy_parameters/` to
`packages/regpoly-legacy/shared/legacy_parameters/`. Six of the seven
pre-v2 YAML configs that used `legacy_file:` move from
`shared/yaml/equidist/` to `packages/regpoly-legacy/shared/yaml/equidist/`
(the seventh, `carry32.yaml`, moves too — `test_seek_config` now uses a
hand-authored inline fixture).

### Pure-C++ build mode (no Python / pybind11 needed)

`packages/regpoly-cpp/CMakeLists.txt` adds a new option:

```
-DREGPOLY_BUILD_PYTHON_EXTENSION=ON   # default: builds _regpoly_cpp.so
-DREGPOLY_BUILD_PYTHON_EXTENSION=OFF  # pure-C++: skips find_package(pybind11)
                                      # and skips the _regpoly_cpp target
```

When OFF, the build produces only `libregpoly_core.a`, the `regpoly-cli`
binary, the public headers, and the `find_package(regpoly)` CMake
config package — no pybind11 dependency at configure, build, or install
time. The Python wheel build (scikit-build-core / `uv sync`) still
defaults to ON; nothing changes for Python consumers.

C++-only consumers can now do:

```
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF
cmake --build build
sudo cmake --install build
```

and link their own project via `find_package(regpoly REQUIRED)` →
`target_link_libraries(my_proj PRIVATE regpoly::regpoly)` without any
Python tooling on the system.

### `gen_primitive_factors_data.py` moved to `regpoly/tools/`

The dev-time codegen script that produces
`packages/regpoly-cpp/src/algebra/primitive_factors_data.cpp` from the
canonical JSON factor table moved from `packages/regpoly-cpp/scripts/`
to `packages/regpoly/src/regpoly/tools/`. Both regen scripts
(`generate_factors.py` produces the JSON; `gen_primitive_factors_data.py`
converts JSON → `.cpp`) now live side-by-side in `regpoly/tools/`.

Invocation: `python -m regpoly.tools.gen_primitive_factors_data`.
The generated `.cpp` is unchanged and remains committed to the repo, so
fresh builds don't depend on running either script.

## v2.0.0 — 2026-05-04

The v2.0 redesign concludes the 9-phase migration that pushed all
algorithmic logic to C++ (`regpoly_cpp`) and made Python (`regpoly`)
and the FastAPI web app (`regpoly_web`) thin wrappers over that
core. A C++-only user has full feature parity with a Python user.

### Architecture

- **Three-package `uv` workspace** with strict one-way deps enforced
  by `import-linter`: `regpoly_web → regpoly → regpoly_cpp`. The
  web app no longer imports `_regpoly_cpp` directly; everything
  goes through `regpoly` shims (`regpoly.introspection`,
  `regpoly.analyses.pis.analyze_single_generator`, …).
- **C++ surface**: `CombinedGenerator`, `Catalog`, all search
  drivers (PrimitiveSearch, Seek, TemperingOptimizer,
  TemperingSearch), the legacy `.dat` reader, and a vendored
  `yaml-cpp 0.8.0` for YAML config parsing.
- **`regpoly-cli`** standalone binary with subcommands: `catalog
  list/show/gen`, `legacy-info`, `legacy-trans`, `search`, `show`,
  `publish`. Equivalent C++ surface to the Python `regpoly` CLI.
- **Python shims**: `regpoly.library/__init__.py` shrunk from 606
  → 77 LOC, `regpoly.io.legacy_reader` 505 → 45 LOC; both delegate
  to the C++ implementations through pybind11.

### Web app

- **SQLite schema v2** with three typed result tables
  (`equidistribution_result`, `collision_free_result`,
  `tuplets_result`) mirroring the C++ analysis structs 1:1. Old v1
  databases auto-upgrade in place at startup
  (`regpoly_web.migrations.v2.migrate_v2_inplace`); the legacy
  `test_result` table is preserved for migrated rows. Standalone
  CLI `packages/regpoly-web/scripts/migrate_v2.py` for offline
  rehearsals.
- **FastAPI TestClient + 8 Playwright golden-path e2e tests**
  (opt-in via `@pytest.mark.e2e`) covering families page,
  tested-generator detail, primitive-search list, tempering cancel,
  publish/unpublish library_id surfacing, YAML import, SSE progress.

### Documentation

- **MkDocs Material site** with 17 per-family theory pages, an
  auto-stubbed paper page per catalog YAML
  (`docs/gen_paper_pages.py` + `mkdocs-gen-files`), 18 per-family
  notebooks generated from a parametric template
  (`docs/notebooks/families/_stamp.py` + `mkdocs-jupyter`), and an
  expanded `theory/f2-linear-generators.md` umbrella covering
  state evolution, characteristic polynomial, equidistribution +
  dimension gap, combined generators, and tempering.
- Continuous deployment to GitHub Pages on every push to `master`
  (`.github/workflows/docs.yml`).

### Testing

- **C++**: GoogleTest via CMake `FetchContent`, `ctest --test-dir
  build`. 125 tests at v2.0.0.
- **Python**: pytest split per-package; 84 default-lane tests
  (was 52 at Phase 0); MTToolBox cross-check + nbmake notebook
  exec gated by `@pytest.mark.slow`.
- **Web e2e**: 8 Playwright tests gated by `@pytest.mark.e2e`.

### Backwards-compatibility

- Legacy family aliases (`MersenneTwister` → `MTGen`, `MT` →
  `MTGen`, `MELG` → `MELGGen`, `TGFSR` → `TGFSRGen`, `PolyLCG` →
  `PolyLCGGen`, …) are **retained** in `regpoly.core.generator
  ::_FAMILY_ALIASES` and the C++ factory. Existing YAML configs
  under `shared/yaml/` and the test fixtures continue to work.
  These aliases will be removed in a future major release once
  downstream consumers have migrated.
- Legacy `.dat` parameter file format remains supported via
  `regpoly_cpp::io::legacy_reader`.

### Migration history

The 9-phase commit chain (all on `origin/master`):

- Phase 0: scaffolding (`952237e`, `710eb09`, `db07422`)
- Phase 1: CombinedGenerator + public C++ API (`dec556e`, `8e08fb5`)
- Phase 2: orchestration loops to C++ (10 commits)
- Phase 3: Catalog + YAML + legacy reader to C++ (`cfb0a07`,
  `b0e163c`, `8a0023b`)
- Phase 4: regpoly-cli feature-complete (`2e89ec7` → `ef77d7a`)
- Phase 5: web rewire (`c795c85` → `0974c97`)
- Phase 6: documentation content (`60db8c8` → `70fd446`)
- Phase 7: per-family notebooks (`56d6678`)
- Phase 8: cleanup + GH Pages deploy + v2.0.0 tag (THIS RELEASE)
