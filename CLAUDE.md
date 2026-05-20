# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**regpoly_monorepo** is the canonical home of the REGPOLY codebase (v2.0+) — a research/academic toolkit for analysing and searching for high-quality combined pseudo-random number generators (PRNGs) based on modulo-2 linear recurrences (LFSRs over GF(2)). The codebase is split into three Python packages bound by a `uv` workspace:

| Package | Role | Build |
|---|---|---|
| `packages/regpoly-cpp/` | Native C++ core: GF(2) linear algebra, MT-family generators, lattice methods. Ships either as a Python wheel (`_regpoly_cpp` extension + thin `regpoly_cpp` package, via scikit-build-core) **or** as a pure-C++ install (`libregpoly_core.a` + `regpoly-cli` + `find_package(regpoly)` headers, via plain CMake with `-DREGPOLY_BUILD_PYTHON_EXTENSION=OFF` — no pybind11 needed). | scikit-build-core + CMake (optional pybind11) |
| `packages/regpoly/`     | Pure-Python library: YAML config, generator catalogue, search loop, equidistribution analyses, CLI (`regpoly`). Depends on `regpoly-cpp`. | setuptools |
| `packages/regpoly-web/` | FastAPI web shell + SQLite-backed result store. CLI (`regpoly-web`). Depends on `regpoly`. | setuptools |
| `packages/regpoly-legacy/` | **Optional** pure-Python add-on: reads pre-v2 `.dat` parameter files via `regpoly_legacy.LegacyReader` and exposes a `regpoly-legacy` CLI (`info` / `trans` / `seek`). Constructs `Generator` objects through the existing `regpoly` factory. Depends on `regpoly`. No C++ code of its own. | setuptools |

**Dependency arrows are strictly one-way.** Two layered contracts in import-linter:
- `regpoly-web → regpoly → regpoly-cpp` (catalog / YAML / web path)
- `regpoly-legacy → regpoly → regpoly-cpp` (pre-v2 `.dat` path)

No cross-imports in the other direction. `regpoly-legacy` and `regpoly-web` never import each other. If a future change tempts you to make `regpoly-cpp` import from `regpoly` (or `regpoly` import from `regpoly-legacy`), that's a sign the boundary is wrong — push the abstraction down instead.

## Repository Layout

```
regpoly_monorepo/
├── pyproject.toml                       workspace declaration + dependency-groups
│                                          (docs/lint/test/e2e), ruff, pytest,
│                                          import-linter contracts
├── README.md
├── CLAUDE.md                            this file
├── .github/workflows/                   ci.yml, ci-slow.yml, docs.yml
├── packages/
│   ├── regpoly-cpp/
│   │   ├── pyproject.toml
│   │   ├── CMakeLists.txt               C++ build via scikit-build-core (Python wheel mode)
│   │   │                                  OR plain cmake -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF
│   │   │                                  (pure-C++ mode, no pybind11 needed).
│   │   │                                  -DREGPOLY_BUILD_TESTS=ON enables GoogleTest via ctest.
│   │   ├── src/                         algebra/ analyses/ bindings/ cli/ core/
│   │   │   │                              generators/ lattice/ library/ search/
│   │   │   │                              transforms/ yaml_config/ + include/<subdir>/
│   │   │   └── regpoly_cpp/__init__.py  thin Python re-export of the .so
│   │   │                                  (only shipped in wheel mode)
│   │   └── tests/                       GoogleTest C++ tests (ctest) +
│   │       ├── CMakeLists.txt           FetchContent googletest
│   │       ├── test_smoke.cpp et al.
│   │       └── python/                  pytest tests for the pybind11 binding ABI
│   ├── regpoly/
│   │   ├── pyproject.toml
│   │   ├── src/regpoly/                 core/ analyses/ search/ io/ data/
│   │   │                                  library/ tools/ __init__.py
│   │   └── tests/                       pytest tests for the Python wrapper layer
│   │                                      (incl. MTToolBox cross-check, marked @pytest.mark.slow)
│   ├── regpoly-web/
│   │   ├── pyproject.toml
│   │   ├── src/regpoly_web/             FastAPI app + routes + tasks + templates
│   │   │                                  (PostgreSQL via psycopg3 since v2.x dockerize cutover)
│   │   ├── scripts/                     regpoly-web.sh, wipe-db.sh
│   │   └── tests/                       pytest + Playwright e2e (in the e2e marker)
│   └── regpoly-legacy/                  optional pure-Python add-on (since post-v2.0)
│       ├── pyproject.toml
│       ├── src/regpoly_legacy/          reader / seek_factory / cli / tools/
│       ├── shared/                      legacy_parameters/ + yaml/equidist/ for `.dat` flows
│       └── tests/python/                + tests/python/data/golden/ frozen snapshots
├── third_party/
│   ├── MTToolBox/                       vendored; used only by the slow-lane cross-check
│   └── dSMFT/                           vendored reference for dSFMT recurrence (comments only)
├── shared/                              runtime data only
│   └── yaml/                            YAML configs for the search CLI
│                                          (pre-v2 .dat fixtures live under
│                                           packages/regpoly-legacy/shared/legacy_parameters/)
├── docs/                                MkDocs Material site (workspace-root)
│   ├── mkdocs.yml
│   ├── index.md
│   ├── theory/                          F₂-linear theory + algorithm specs (was shared/docs and
│   │                                      shared/cpp_regpoly_docs/*.md)
│   ├── generators/                      one .md per family (was shared/cpp_regpoly_docs/generators/)
│   ├── library/                         catalog YAML + index (was shared/cpp_regpoly_docs/library/);
│   │                                      this is the **runtime catalog source** — the web app and
│   │                                      C++ Catalog both read from here. CMake will install it to
│   │                                      ${prefix}/share/regpoly/catalog/.
│   ├── papers/                          PDFs + bibliography (was shared/papers/)
│   ├── usage/                           {python,cpp,web,notebooks}.md
│   ├── dev/                             architecture, building, contributing, postgres
│   └── notebooks/                       per-family + research + utility notebooks
└── notebooks/                           pre-v2 notebooks kept for reference
```

## Build & Run

```bash
cd /home/frpan/projets/claude_projects/regpoly_monorepo

# Install every workspace member (editable) and build the C++ extension:
uv sync

# Run the search CLI on a workspace YAML:
uv run regpoly shared/yaml/equidist/mt19937.yaml

# Run the legacy `.dat` reader (add-on; takes positional args):
uv run regpoly-legacy seek 1 \
  packages/regpoly-legacy/shared/legacy_parameters/example1c \
  packages/regpoly-legacy/shared/legacy_parameters/96_1.dt

# Run the web app (Postgres-only; default DSN comes from $REGPOLY_DB_URL):
uv run regpoly-web                                       # uses $REGPOLY_DB_URL
uv run regpoly-web --db-url postgresql://user:pw@host/regpoly

# Run tests (default skips slow + e2e):
uv run pytest                         # all default-lane tests across packages
uv run pytest -m slow                 # MTToolBox cross-checks, nbmake notebooks
uv run pytest -m e2e                  # Playwright golden paths

# C++ tests (default Python-wheel build mode):
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure

# C++ tests (pure-C++ mode, no pybind11 needed):
cmake -S packages/regpoly-cpp -B build-cpp \
      -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF \
      -DREGPOLY_BUILD_TESTS=ON
cmake --build build-cpp -j
ctest --test-dir build-cpp --output-on-failure

# Docs site:
uv sync --group docs
cd docs && uv run mkdocs serve
```

The `tests/_mttoolbox_build.md` doc inside the regpoly package describes how to (re-)build `MTToolBox/samples/*/calc_equidist` binaries used by the cross-check tests; the `MTToolBox/Makefile`s have been refreshed to point at this monorepo's path. **Bypass autotools** when building those samples — the documented `g++` recipe works directly.

## Specs (read these when relevant)

- `docs/theory/equidistribution-spec.md` — **the design of record** for matricial equidistribution computation on F_2-linear generators that are not assumed full-period (Berlekamp–Massey → factor χ_f → invariant subspace → matricial DE core → guard → verify). Read whenever working on equidistribution where the characteristic polynomial may not be primitive.
- `docs/theory/antithetic-check.md` — algorithm for testing local antitheticity of a linear RNG point set.

## Past redesign — v2.0 plan (COMPLETE)

The 9-phase v2.0 redesign is **COMPLETE**. Tagged `v2.0.0` at commit `47c24e6` on 2026-05-04. Plan file: `~/.claude/plans/i-want-to-have-rippling-starlight.md` (preserved for reference). Final architecture: all algorithmic logic in C++; Python is a thin wrapper; the web app uses only Python (no `_cpp` imports); C++-only users have full feature parity via `regpoly-cli`, and can build the C++ core without any Python/pybind11 tooling via the CMake option `-DREGPOLY_BUILD_PYTHON_EXTENSION=OFF` (added v2.x post-release).

## Active plan — Simplify to three containers

Plan file: `~/.claude/plans/the-architecture-of-the-streamed-lollipop.md`. Goal: three-service compose stack (`db`, `web`, `worker`) deployable with `docker compose up -d`. No init container (migrations run in the FastAPI lifespan), no Caddy (no public exposure baked in — external access is handled out of stack by whatever fronts it, e.g. cloudflared on a NAS), no in-stack backups (named `pgdata` volume snapshotted externally), no host scripts. Image is multi-arch (`linux/amd64` + `linux/arm64`) so it runs on any common Docker host. Per-phase commit + push to `origin/master` is authorized **for the duration of this plan only** (overrides the general "no auto commit/push" rule).

## History (for context — don't act on it)

This monorepo was carved out of a larger research repo at `/home/frpan/projets/claude_projects/regpoly/` (still on disk, not part of this workspace). That parent repo contains:

- `working_code/` — original C library + LaTeX-driven build pipeline (`tcode` extracts headers from `.tex`). Predates this monorepo. Not maintained from here. Reference only.
- `MinimalCode/{c_regpoly,pkg_regpoly,pyregpoly}/` — sibling reduced re-implementations (C, Python). Independent of this monorepo. Each has its own README.
- `NEWPROJECT/` — a parallel C reorganisation building `POL`, `QMC`, `PROD`. Not used by this monorepo.

If you need to read those for cross-reference, do so via absolute paths under `/home/frpan/projets/claude_projects/regpoly/`. **Never modify them from this session** — they have their own owners and workflows.

The previous home of the C++ port was `MinimalCode/cpp_regpoly/`. Its contents have been moved into this monorepo wholesale: `src/cpp/` → `packages/regpoly-cpp/src/`, `src/python/web/` → `packages/regpoly-web/src/regpoly_web/`, the rest of `src/python/` → `packages/regpoly/src/regpoly/`, and the bundled `MTToolBox/` and `dSMFT/` into `third_party/`.

## Conventions

- **Math notation in docs.** Sphinx + MyST is configured with the `dollarmath` and `amsmath` extensions (see `docs/conf.py`'s `myst_enable_extensions`). Math symbols MUST be rendered as proper math, never as inline code. Use `$...$` for inline math (e.g. `$t(s) = m - d_{\max}(s)$`, not `` `t(s) = m - d_max(s)` ``) and `$$...$$` (or fenced ```` ```{math} ```` blocks) for display math. Bold vectors with `\mathbf{...}`, blackboard for fields (`\mathbb{F}_2`). Mirror the conventions documented in [`docs/theory/notation.md`](docs/theory/notation.md) and demonstrated throughout `docs/theory/equidistribution-spec.md`. Backtick-styled code is reserved for code identifiers (function names, file paths, parameter names) — never for variables, sets, formulas, or anything you'd typeset in LaTeX.
- **One-way deps.** Two layered import-linter contracts: `regpoly-web → regpoly → regpoly-cpp` (catalog / YAML / web path) and `regpoly-legacy → regpoly → regpoly-cpp` (pre-v2 `.dat` path). Nothing goes backwards. Defined in the workspace `pyproject.toml`.
- **`shared/` is workspace-shared runtime data only.** YAML search configs live under `shared/yaml/`. Pre-v2 `.dat` fixtures are NOT here anymore — they moved to `packages/regpoly-legacy/shared/legacy_parameters/`.
- **C++ extension stays in `regpoly-cpp`.** The `.so` is a build artefact of `regpoly-cpp` (when built in Python-wheel mode); `regpoly` consumes it via the `regpoly_cpp` package interface only.
- **Pure-C++ users.** `cmake -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF` builds `libregpoly_core.a` + `regpoly-cli` + headers + `find_package(regpoly)` config — no pybind11, no Python.
- **Git is per-monorepo, not per-package.** Run `git init` once at the workspace root.

## Deferred decisions (call out before changing)

- **Workspace lockfile.** `uv.lock` will appear after the first `uv sync`. Commit it.
- **`MTToolBox/Makefile`** absolute paths were re-baked via `config.status` after the move; a future relocation needs another `cd third_party/MTToolBox && ./config.status --recheck`.

## Auto-memory

User auto-memory lives at `~/.claude/projects/-home-frpan-projets-claude-projects-regpoly-monorepo/memory/`. Preserved feedback rules:

- *No auto commit/push* in general — but the active v2.0 redesign explicitly authorizes per-phase commit + push for its duration; see `MEMORY.md → project_redesign_v2.md`.
- *Safe masks only for the tempering optimiser*.
