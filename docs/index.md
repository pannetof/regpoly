---
sd_hide_title: true
---

# REGPOLY — C++/Python toolkit for combined PRNGs

::::{grid} 1 1 2 2
:gutter: 3
:padding: 4 4 0 0

:::{grid-item}
:columns: 12

# REGPOLY

A research toolkit for analysing and searching for high-quality combined
pseudo-random number generators based on modulo-2 linear recurrences (LFSRs
over GF(2)). Used in F₂-linear generator design and validation.
:::

::::

---

## Get started

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} 🐍 &nbsp; Get started in Python
:link: usage/python
:link-type: doc
:class-card: sd-shadow-sm

Install with `uv sync`, run a search on an MT19937 YAML config in one
command, or drive the library programmatically with `Catalog`,
`Combination`, `analyze_single_generator`.
:::

:::{grid-item-card} ⚙ &nbsp; Use REGPOLY from C++
:link: usage/cpp
:link-type: doc
:class-card: sd-shadow-sm

Pure-C++ build with no Python or pybind11 (`cmake -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF`).
Link via `find_package(regpoly)`. A runnable "hello equidistribution"
example is on the C++ usage page.
:::

:::{grid-item-card} 📓 &nbsp; Browse the notebooks
:link: notebooks/index
:link-type: doc
:class-card: sd-shadow-sm

29 worked examples — 15 per-family demos (MT, SFMT, WELL, MELG, …) plus
14 reference notebooks covering primitivity, equidistribution, charpoly
derivation, …
:::

::::

---

## Documentation

::::{grid} 1 2 2 4
:gutter: 3

:::{grid-item-card} 📐 &nbsp; Theory
:link: theory/index
:link-type: doc

F₂-linear framework, equidistribution, lattice methods, tempering,
LCP irreducibility, search format.
:::

:::{grid-item-card} 🧮 &nbsp; Generators
:link: generators/index
:link-type: doc

17 family pages — MT, SFMT, dSFMT, MTGP, TinyMT32, RMT64, WELL,
Tausworthe, TGFSR, PolyLCG, F2w-LFSR, F2w-PolyLCG, Marsaglia XOR,
Xoroshiro, Xoshiro, MELG, Cellular Automata.
:::

:::{grid-item-card} 📚 &nbsp; Library
:link: library/index
:link-type: doc

Published parameter sets indexed by paper.
:::

:::{grid-item-card} 📄 &nbsp; Papers
:link: papers/index
:link-type: doc

Bibliography + PDFs of the foundational papers behind each family.
:::

:::{grid-item-card} 🔍 &nbsp; API Reference
:link: api/python/index
:link-type: doc

Auto-generated reference for the public Python surface (`autodoc` +
NumPy-style docstrings) and the C++ surface (Doxygen + Exhale's
Class / File / Namespace hierarchy). Bidirectionally cross-linked.
:::

:::{grid-item-card} 🛠 &nbsp; Developer
:link: dev/architecture
:link-type: doc

Architecture, building, contributing, API doc style guide, web app
design spec, dev-env setup.
:::

:::{grid-item-card} 🚀 &nbsp; Deployment
:link: deployment
:link-type: doc

Three-container Docker stack for the web app (db / web / worker),
deployable on any Docker host.
:::

:::{grid-item-card} 🗒 &nbsp; Release notes
:link: changelog
:link-type: doc

Per-release changelog.
:::

::::

---

## Packages

REGPOLY is a `uv` workspace of four packages:

| Package | Role |
|---|---|
| [`regpoly-cpp`](packages/regpoly-cpp.md) | C++ core: GF(2) algebra, lattice methods, all generator families, search drivers, YAML, catalog. Publishes the `_regpoly_cpp` Python extension and a standalone `regpoly-cli` binary. Can also be built **without** Python/pybind11 (`-DREGPOLY_BUILD_PYTHON_EXTENSION=OFF`). |
| [`regpoly`](packages/regpoly.md) | Pure-Python wrapper around `regpoly_cpp`, plus result decoration (DataFrames, plots) and the Python `regpoly` CLI. |
| [`regpoly-web`](packages/regpoly-web.md) | FastAPI web UI for launching searches and browsing results (PostgreSQL via `psycopg3`). Uses `regpoly` only. |
| [`regpoly-legacy`](packages/regpoly-legacy.md) | Optional pure-Python add-on that reads the pre-v2 `.dat` parameter format and builds `Generator` objects through the same C++ factory. |

Dependency arrows are strictly one-way:

- `regpoly-web → regpoly → regpoly-cpp`
- `regpoly-legacy → regpoly → regpoly-cpp`

Both layered contracts are enforced by `import-linter` in CI.

```{toctree}
:hidden:
:caption: User Guide

theory/index
generators/index
library/index
papers/index
usage/python
usage/cpp
usage/python-cpp-bridge
usage/web
notebooks/index
deployment
packages/regpoly-cpp
packages/regpoly
packages/regpoly-web
packages/regpoly-legacy
```

```{toctree}
:hidden:
:caption: API Reference

api/python/index
api/cpp/library/library_root
```

```{toctree}
:hidden:
:caption: Developer

dev/architecture
dev/building
dev/contributing
dev/api-docs
dev/postgres
dev/testing-lanes
dev/mttoolbox-build
dev/regpoly-web-design-spec
dev/web-user-stories
```

```{toctree}
:hidden:
:caption: Reference

changelog
```
