# Changelog

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
