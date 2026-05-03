# REGPOLY

A research toolkit for analysing and searching for high-quality combined pseudo-random number generators based on modulo-2 linear recurrences (LFSRs over GF(2)).

## Quick links

- **[F₂-linear generator theory](theory/f2-linear-generators.md)** — start here.
- **[Generator catalogue](generators/index.md)** — one page per family (MT, SFMT, WELL, MELG, …).
- **[Library](library/index.md)** — published parameter sets, indexed by paper.
- **[Papers](papers/index.md)** — bibliography and PDFs.
- **[Usage](usage/python.md)** — how to drive REGPOLY from Python, C++, the web UI, or notebooks.

## Layout

REGPOLY is a `uv` workspace of three packages:

| Package | Role |
|---|---|
| `regpoly-cpp` | C++ core: GF(2) algebra, lattice methods, all generator families, search drivers, YAML, catalog. Publishes the `_regpoly_cpp` Python extension and a standalone `regpoly-cli` binary. |
| `regpoly` | Pure-Python wrapper around `regpoly_cpp`, plus result decoration (DataFrames, plots) and the Python `regpoly` CLI. |
| `regpoly-web` | FastAPI web UI for launching searches and browsing results. Uses `regpoly` only. |

Dependency arrows are strictly one-way: `regpoly-web → regpoly → regpoly-cpp`.
