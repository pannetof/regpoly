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
├── pyproject.toml                       workspace declaration only (no [project])
├── README.md
├── CLAUDE.md                            this file
├── packages/
│   ├── regpoly-cpp/
│   │   ├── pyproject.toml
│   │   ├── CMakeLists.txt               drives the C++ build via scikit-build-core
│   │   ├── src/                         was cpp_regpoly/src/cpp/{algebra,core,generators,transforms,lattice,bindings,include}/
│   │   │   └── regpoly_cpp/__init__.py  thin Python re-export of the .so
│   │   └── tests/                       (empty — C++ tests via ctest go here)
│   ├── regpoly/
│   │   ├── pyproject.toml
│   │   ├── src/regpoly/                 was cpp_regpoly/src/python/{core,library,analyses,search,io,data,tools,__init__.py}
│   │   └── tests/                       was cpp_regpoly/tests/* (ALL pytest tests live here for now — split per-package later if useful)
│   └── regpoly-web/
│       ├── pyproject.toml
│       ├── src/regpoly_web/             was cpp_regpoly/src/python/web/
│       ├── scripts/                     regpoly-web.sh, wipe-db.sh
│       └── var/                         regpoly.db lives here (user data, gitignored)
├── third_party/
│   ├── MTToolBox/                       vendored MTToolBox; only used by regpoly/tests/test_mttoolbox_crosscheck.py
│   └── dSMFT/                           vendored reference for dSFMT recurrence (used as comment reference only)
├── shared/
│   ├── yaml/                            YAML configs for the search CLI (equidist, fullperiodsearch, generators, …)
│   ├── legacy_parameters/               legacy *.dat / example* parameter files (input fixtures)
│   ├── papers/                          papers served by the web app's /papers endpoint
│   ├── docs/                            cross-cutting algorithm specs — see "Specs" below
│   └── cpp_regpoly_docs/                ports of the old per-module documentation
└── notebooks/                           research notebooks — not packaged
```

`packages/regpoly/tests/` currently holds the entire pytest suite from `cpp_regpoly`. It's fine to leave it here until the boundaries naturally pull tests apart; do not pre-emptively split.

## Build & Run

```bash
cd /home/frpan/projets/claude_projects/regpoly_monorepo

# Install everything in the workspace (editable), build the C++ extension:
uv sync

# Run the search CLI:
uv run regpoly shared/yaml/equidist/well19937a.yml

# Run the web app:
uv run regpoly-web --db packages/regpoly-web/var/regpoly.db

# Run tests (default skips slow MTToolBox cross-checks):
uv run pytest packages/regpoly/tests
uv run pytest -m slow packages/regpoly/tests        # cross-checks too
```

The `tests/_mttoolbox_build.md` doc inside the regpoly package describes how to (re-)build `MTToolBox/samples/*/calc_equidist` binaries used by the cross-check tests; the `MTToolBox/Makefile`s have been refreshed to point at this monorepo's path. **Bypass autotools** when building those samples — the documented `g++` recipe works directly.

## Specs (read these when relevant)

- `shared/docs/equidistribution-spec.md` — **the design of record** for matricial equidistribution computation on F_2-linear generators that are not assumed full-period (Berlekamp–Massey → factor χ_f → invariant subspace → matricial DE core → guard → verify). Read whenever working on equidistribution where the characteristic polynomial may not be primitive.
- `shared/docs/antithetic-check.md` — algorithm for testing local antitheticity of a linear RNG point set.

## History (for context — don't act on it)

This monorepo was carved out of a larger research repo at `/home/frpan/projets/claude_projects/regpoly/` (still on disk, not part of this workspace). That parent repo contains:

- `working_code/` — original C library + LaTeX-driven build pipeline (`tcode` extracts headers from `.tex`). Predates this monorepo. Not maintained from here. Reference only.
- `MinimalCode/{c_regpoly,pkg_regpoly,pyregpoly}/` — sibling reduced re-implementations (C, Python). Independent of this monorepo. Each has its own README.
- `NEWPROJECT/` — a parallel C reorganisation building `POL`, `QMC`, `PROD`. Not used by this monorepo.

If you need to read those for cross-reference, do so via absolute paths under `/home/frpan/projets/claude_projects/regpoly/`. **Never modify them from this session** — they have their own owners and workflows.

The previous home of the C++ port was `MinimalCode/cpp_regpoly/`. Its contents have been moved into this monorepo wholesale: `src/cpp/` → `packages/regpoly-cpp/src/`, `src/python/web/` → `packages/regpoly-web/src/regpoly_web/`, the rest of `src/python/` → `packages/regpoly/src/regpoly/`, and the bundled `MTToolBox/` and `dSMFT/` into `third_party/`.

## Conventions

- **One-way deps.** `regpoly-web` may import `regpoly`; `regpoly` may import `regpoly_cpp`; nothing goes backwards. Add `import-linter` to CI if drift becomes a risk.
- **`shared/` is data, not code.** Anything code-shared between packages = a fourth package, not a `shared/` import.
- **`var/` and `regpoly.db` are user data.** Never commit DB files or move them into `src/`.
- **C++ extension stays in `regpoly-cpp`.** The `.so` is a build artefact of `regpoly-cpp`; `regpoly` consumes it via the `regpoly_cpp` package interface only.
- **Git is per-monorepo, not per-package.** Run `git init` once at the workspace root.

## Deferred decisions (call out before changing)

- **Test split.** Tests are unsplit on purpose. Move tests to a per-package `tests/` only when you have a concrete reason (e.g. CI job differentiation, or a test that genuinely belongs to web).
- **Workspace lockfile.** `uv.lock` will appear after the first `uv sync`. Commit it.
- **`MTToolBox/Makefile`** absolute paths were re-baked via `config.status` after the move; a future relocation needs another `cd third_party/MTToolBox && ./config.status --recheck`.

## Auto-memory

User auto-memory has been migrated into this project's memory directory (`~/.claude/projects/-home-frpan-projets-claude-projects-regpoly_monorepo/memory/`). The two preserved feedback rules — *no auto commit/push* and *safe masks only for the tempering optimiser* — apply unchanged here.
