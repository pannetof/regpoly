# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**regpoly_monorepo** is the canonical home of the REGPOLY codebase (v2.0+) — a research/academic toolkit for analysing and searching for high-quality combined pseudo-random number generators (PRNGs) based on modulo-2 linear recurrences (LFSRs over GF(2)). The codebase is split into three Python packages bound by a `uv` workspace:

| Package | Role | Build |
|---|---|---|
| `packages/regpoly-cpp/` | Native C++ core: GF(2) linear algebra, MT-family generators, lattice methods, pybind11 bindings. Publishes the `_regpoly_cpp` extension and a thin `regpoly_cpp` Python wrapper. | scikit-build-core + CMake |
| `packages/regpoly/`     | Pure-Python library: YAML config, generator catalogue, search loop, equidistribution analyses, CLI (`regpoly`). Depends on `regpoly-cpp`. | setuptools |
| `packages/regpoly-web/` | FastAPI web shell + SQLite-backed result store. CLI (`regpoly-web`). Depends on `regpoly`. | setuptools |

**Dependency arrows are strictly one-way:** `regpoly-web → regpoly → regpoly-cpp`. No cross-imports in the other direction. If a future change tempts you to make `regpoly-cpp` import from `regpoly`, that's a sign the boundary is wrong — push the abstraction down into `regpoly-cpp` instead.

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
│   │   ├── CMakeLists.txt               C++ build via scikit-build-core; -DREGPOLY_BUILD_TESTS=ON
│   │   │                                  enables the GoogleTest target via ctest
│   │   ├── src/                         {algebra,core,generators,transforms,lattice,bindings,include}/
│   │   │   └── regpoly_cpp/__init__.py  thin Python re-export of the .so
│   │   └── tests/                       GoogleTest C++ tests (ctest) +
│   │       ├── CMakeLists.txt           FetchContent googletest
│   │       ├── test_smoke.cpp
│   │       └── python/                  pytest tests for the pybind11 binding ABI
│   ├── regpoly/
│   │   ├── pyproject.toml
│   │   ├── src/regpoly/                 {core,library,analyses,search,io,data,tools,__init__.py}
│   │   └── tests/                       pytest tests for the Python wrapper layer
│   │                                      (incl. MTToolBox cross-check, marked @pytest.mark.slow)
│   └── regpoly-web/
│       ├── pyproject.toml
│       ├── src/regpoly_web/             FastAPI web app (was cpp_regpoly/src/python/web/)
│       ├── scripts/                     regpoly-web.sh, wipe-db.sh
│       ├── var/                         regpoly.db (user data, gitignored)
│       └── tests/                       pytest + (Phase 5) Playwright e2e
├── third_party/
│   ├── MTToolBox/                       vendored; used only by the slow-lane cross-check
│   └── dSMFT/                           vendored reference for dSFMT recurrence (comments only)
├── shared/                              runtime data only
│   ├── yaml/                            YAML configs for the search CLI
│   └── legacy_parameters/               legacy *.dat / example* parameter files (test fixtures)
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
│   ├── dev/                             architecture, building, contributing
│   └── notebooks/                       (Phase 7) per-family + research + utility notebooks
└── notebooks/                           legacy notebooks (Phase 7 will rewrite under docs/notebooks/)
```

## Build & Run

```bash
cd /home/frpan/projets/claude_projects/regpoly_monorepo

# Install everything in the workspace (editable), build the C++ extension:
uv sync

# Run the search CLI:
uv run regpoly shared/yaml/equidist/well19937a.yml

# Run the web app:
uv run regpoly-web --db packages/regpoly-web/var/regpoly.db

# Run tests (default skips slow + e2e):
uv run pytest                         # all default-lane tests across packages
uv run pytest -m slow                 # MTToolBox cross-checks, nbmake notebooks
uv run pytest -m e2e                  # Playwright golden paths (after Phase 5)

# C++ tests:
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure

# Docs site:
uv sync --group docs
cd docs && uv run mkdocs serve
```

The `tests/_mttoolbox_build.md` doc inside the regpoly package describes how to (re-)build `MTToolBox/samples/*/calc_equidist` binaries used by the cross-check tests; the `MTToolBox/Makefile`s have been refreshed to point at this monorepo's path. **Bypass autotools** when building those samples — the documented `g++` recipe works directly.

## Specs (read these when relevant)

- `docs/theory/equidistribution-spec.md` — **the design of record** for matricial equidistribution computation on F_2-linear generators that are not assumed full-period (Berlekamp–Massey → factor χ_f → invariant subspace → matricial DE core → guard → verify). Read whenever working on equidistribution where the characteristic polynomial may not be primitive.
- `docs/theory/antithetic-check.md` — algorithm for testing local antitheticity of a linear RNG point set.

## Past redesign — v2.0 plan (COMPLETE)

The 9-phase v2.0 redesign is **COMPLETE**. Tagged `v2.0.0` at commit `47c24e6` on 2026-05-04. Plan file: `~/.claude/plans/i-want-to-have-rippling-starlight.md` (preserved for reference). Final architecture: all algorithmic logic in C++; Python is a thin wrapper; the web app uses only Python (no `_cpp` imports); C++-only users have full feature parity via `regpoly-cli`.

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

- **One-way deps.** `regpoly-web` may import `regpoly`; `regpoly` may import `regpoly_cpp`; nothing goes backwards. Enforced by `import-linter` (contracts in workspace `pyproject.toml`).
- **`shared/` is runtime data only.** YAML config templates and legacy `.dat` test fixtures. Documentation, papers, and the catalog live under `docs/`.
- **`var/` and `regpoly.db` are user data.** Never commit DB files or move them into `src/`.
- **C++ extension stays in `regpoly-cpp`.** The `.so` is a build artefact of `regpoly-cpp`; `regpoly` consumes it via the `regpoly_cpp` package interface only.
- **Git is per-monorepo, not per-package.** Run `git init` once at the workspace root.

## Deferred decisions (call out before changing)

- **Workspace lockfile.** `uv.lock` will appear after the first `uv sync`. Commit it.
- **`MTToolBox/Makefile`** absolute paths were re-baked via `config.status` after the move; a future relocation needs another `cd third_party/MTToolBox && ./config.status --recheck`.

## Auto-memory

User auto-memory lives at `~/.claude/projects/-home-frpan-projets-claude-projects-regpoly-monorepo/memory/`. Preserved feedback rules:

- *No auto commit/push* in general — but the active v2.0 redesign explicitly authorizes per-phase commit + push for its duration; see `MEMORY.md → project_redesign_v2.md`.
- *Safe masks only for the tempering optimiser*.
