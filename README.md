# regpoly --- C++/Python toolkit for combined PRNGs

Analyse the equidistribution properties of pseudo-random number
generators based on modulo-2 linear recurrences (LFSRs over GF(2)).
Uses NTL for polynomial arithmetic, with optional pybind11 bindings for
the Python wrapper layer.

The primary use-case is the full Python stack (`regpoly`, `regpoly-web`,
optional `regpoly-legacy` add-on) built on top of the C++ core. The
C++ core can also be built and installed **standalone**, without Python
or pybind11, for pure C++ consumers.

Search configuration uses YAML; the pre-v2 C-compatible `.dat`
parameter format is supported via the optional pure-Python
`regpoly-legacy` add-on package.

## Dependencies

**Required (always):**
- **libntl-dev** (NTL number theory library, must be installed system-wide)

**Required for the Python wheel build (the default `uv sync` path):**
- Python 3.10+
- pybind11 ≥ 2.12 (pulled in automatically by scikit-build-core)
- PyYAML (pulled in automatically)

**For pure-C++ builds:** none of the Python-side dependencies above are needed; see [docs/usage/cpp.md](docs/usage/cpp.md#pure-c-build-no-python-no-pybind11).

Install NTL on Debian/Ubuntu:

```bash
sudo apt install libntl-dev
```

## Installation

### Python (Python wheel, default)

```bash
cd regpoly_monorepo
uv sync           # builds the wheel via scikit-build-core
```

### Pure C++ (no Python, no pybind11)

```bash
cd regpoly_monorepo
cmake -S packages/regpoly-cpp -B build \
      -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF
cmake --build build -j
sudo cmake --install build
```

Produces `libregpoly_core.a`, the `regpoly-cli` binary, public headers, and the `find_package(regpoly)` CMake config package.

## Usage

### YAML mode (canonical)

```bash
uv run regpoly <config.yaml>
```

YAML configuration files live under `shared/yaml/` (workspace root).
They specify the generator family, parameters, tempering, and test
method in a single file.

**Mersenne Twister MT19937 with tempering:**

```bash
uv run regpoly shared/yaml/equidist/mt19937.yaml
```

**Two-component combined generator:**

```bash
uv run regpoly shared/yaml/equidist/example3.yaml
```

### Legacy `.dat` mode (optional add-on)

The pre-v2 C-compatible `.dat` parameter format is supported via the
`regpoly-legacy` add-on (installed by default in the workspace; ships
a `regpoly-legacy` console script). Fixtures live under
`packages/regpoly-legacy/shared/legacy_parameters/`.

**Single TGFSR generator:**

```bash
cd packages/regpoly-legacy/shared/legacy_parameters
uv run regpoly-legacy seek 1 example1c 96_1.dt
```

**Two-component combined generator:**

```bash
cd packages/regpoly-legacy/shared/legacy_parameters
uv run regpoly-legacy seek 2 example3 96_1.dt 155_2.dt
```

## YAML Configuration Format

A search configuration file ties together generators, tempering, and tests:

```yaml
search:
  seed: [1059674621, 1170794524]
  Lmax: 32

components:
  - generators:
      file: gen.well19937.yaml        # load from file
    tempering:
      - type: tempMK2
        w: 32
        eta: 7
        mu: 15
        u: 11
        l: 18
        b: 0x9d2c5680
        c: 0xefc60000

tests:
  - type: equidistribution
    method: lattice                    # or "matricial"
    max_gap_sum: 100000
```

Generators can also be defined inline:

```yaml
components:
  - generators:
      family: MTGen
      common:
        w: 32
        r: 624
        m: 397
        p: 31
      generators:
        - a: 0x9908b0df
```

## Equidistribution Methods

- **matricial** --- Gaussian elimination on the generator matrix. Fast for
  moderate state sizes.
- **lattice** --- Dual lattice basis reduction using Lenstra's algorithm with
  NTL for polynomial arithmetic. Required for large generators (k > 10000).

## Supported Generator Families

The 17 generator families are exposed under their canonical class names
(with the `-Gen` suffix). Legacy short tags (lowercase, French) and
pre-rename class names are still accepted via the Python `_FAMILY_ALIASES`
table and the C++ factory `reg_alias` entries.

- `TGFSRGen` --- Twisted GFSR
- `MTGen` --- Mersenne Twister
- `PolyLCGGen` --- Polynomial LCG
- `TauswortheGen` --- Tausworthe / LFSR
- `WELLGen` --- WELL / carry generators
- `F2wLFSRGen`, `F2wPolyLCGGen` --- GF(2^w) generators
- `MarsaXorshiftGen` --- Marsaglia Xor-shift
- `CellularAutomataGen` --- Two-rule (rule 90/150) cellular automata (Bhuvaneswari–Bhattacharjee 2026)
- `XoroshiroGen`, `XoshiroGen` --- Blackman–Vigna 2022 scrambled xor-shift
- `MELGGen` --- Maximally Equidistributed F₂-Linear (Harase-Kimoto)
- `SFMTGen`, `DSFMTGen` --- SIMD-Friendly MT (single / double precision)
- `MTGPGen` --- MT for Graphic Processors
- `TinyMT32Gen` --- Tiny MT (127-bit state)
- `RMT64Gen` --- Reducible MT (64-bit)

## Repository Layout

```
regpoly_monorepo/                       (uv workspace root)
├── pyproject.toml                      workspace + lint contracts + testpaths
├── packages/
│   ├── regpoly-cpp/                    C++ core (scikit-build-core wheel)
│   │   ├── pyproject.toml
│   │   ├── CMakeLists.txt              optional pybind11 (REGPOLY_BUILD_PYTHON_EXTENSION)
│   │   └── src/                        algebra/ analyses/ bindings/ cli/ core/
│   │                                     generators/ io/ lattice/ library/
│   │                                     search/ transforms/ yaml_config/
│   ├── regpoly/                        Python wrapper layer (setuptools)
│   │   └── src/regpoly/                core/ search/ analyses/ io/ tools/ library/
│   ├── regpoly-web/                    FastAPI + PostgreSQL (setuptools)
│   │   └── src/regpoly_web/            app + routes + templates + tasks
│   └── regpoly-legacy/                 optional pure-Python .dat reader (setuptools)
│       └── src/regpoly_legacy/         reader/ seek_factory/ cli/ tools/
├── shared/yaml/                        workspace-shared YAML search configs
├── third_party/                        MTToolBox (vendored), dSMFT (reference)
└── docs/                               Sphinx + Doxygen + Breathe + Exhale site
    ├── conf.py                         Sphinx configuration
    ├── Doxyfile                        Doxygen input → XML
    ├── api/                            auto-generated Python + C++ reference
    ├── generators/                     one .md per *Gen family (17 pages)
    ├── library/                        paper YAMLs + per-family params catalog
    ├── notebooks/                      29 worked notebooks (15 family + 14 reference)
    ├── papers/                         reference PDFs + auto-stubs from library/
    ├── theory/                         F₂-linear theory + algorithm specs
    ├── usage/                          python / cpp / python-cpp-bridge / web
    └── dev/                            architecture, building, contributing, api-docs
```

Builds are driven by `pyproject.toml` in each package (no top-level
`setup.py`). The C++ core can also be built and installed standalone
via plain CMake — see [docs/usage/cpp.md](docs/usage/cpp.md).

## Documentation

Full Sphinx site (PyData Sphinx Theme + autodoc + Doxygen + Breathe +
Exhale): **[docs/index.md](docs/index.md)** (deployed to GitHub Pages
on every push to `master` via
[`.github/workflows/docs.yml`](.github/workflows/docs.yml)).

Build locally:

```bash
uv sync --group docs
sudo apt install -y doxygen graphviz          # one-time
uv run --group docs sphinx-build -b html docs site
# or, for live-reload:
uv run --group docs sphinx-autobuild docs site
```

| Section | Where |
|---|---|
| **Theory** — F₂-linear framework, equidistribution, lattice methods, tempering, search format | [`docs/theory/`](docs/theory/index.md) |
| **Generators** — one page per family (17 total) | [`docs/generators/`](docs/generators/index.md) |
| **Usage** — Python, C++ CLI, Python/C++ bridge, web UI | [`docs/usage/python.md`](docs/usage/python.md), [`cpp.md`](docs/usage/cpp.md), [`python-cpp-bridge.md`](docs/usage/python-cpp-bridge.md), [`web.md`](docs/usage/web.md) |
| **Notebooks** — 15 per-family demos + 14 reference notebooks | [`docs/notebooks/`](docs/notebooks/index.md) |
| **Library** — published parameter catalogs | [`docs/library/`](docs/library/index.md) |
| **Papers** — bibliography + PDFs | [`docs/papers/`](docs/papers/index.md) |
| **API reference** — autodoc + Doxygen/Exhale per-symbol pages for the public Python + C++ surface | [`docs/api/python/index.md`](docs/api/python/index.md) |
| **Dev** — architecture, building, contributing, API doc style, Postgres dev-env, web design spec, user stories | [`docs/dev/`](docs/dev/architecture.md) |

### Package READMEs

- [`packages/regpoly-cpp/README.md`](packages/regpoly-cpp/README.md) — C++ core (both build modes).
- [`packages/regpoly/README.md`](packages/regpoly/README.md) — Python wrapper + `regpoly` CLI.
- [`packages/regpoly-web/README.md`](packages/regpoly-web/README.md) — FastAPI web UI.
- [`packages/regpoly-legacy/README.md`](packages/regpoly-legacy/README.md) — optional `.dat` add-on.
- [`packages/regpoly/tests/README.md`](packages/regpoly/tests/README.md) — test lanes + MTToolBox cross-check refresh recipe.
- [`packages/regpoly/tests/_mttoolbox_build.md`](packages/regpoly/tests/_mttoolbox_build.md) — building the upstream `calc_equidist` reference binaries.

### Deployment

- [`DEPLOY.md`](DEPLOY.md) — three-container Docker stack (db / web / worker).
- [`deploy/README.md`](deploy/README.md) — directory index for the `deploy/` artifacts.

### Other top-level docs

- [`CHANGELOG.md`](CHANGELOG.md) — release notes; v2.0.0 (2026-05-04) tagged the
  9-phase C++ rewrite; `Unreleased` covers the `regpoly-legacy` extraction and
  the optional pure-C++ build mode.
- [`CLAUDE.md`](CLAUDE.md) — project instructions for Claude Code sessions
  (workspace conventions, deferred decisions, history).
