# REGPOLY

A research toolkit for analysing and searching for high-quality combined pseudo-random number generators based on modulo-2 linear recurrences (LFSRs over GF(2)).

## Quick links

- **[F₂-linear generator theory](theory/f2-linear-generators.md)** — start here.
- **[Theory index](theory/index.md)** — all the algorithm specs.
- **[Generator catalogue](generators/index.md)** — one page per family (MT, SFMT, WELL, MELG, …).
- **[Library](library/index.md)** — published parameter sets, indexed by paper.
- **[Papers](papers/index.md)** — bibliography and PDFs.
- **[Usage](usage/python.md)** — Python entry point. See also [C++ CLI](usage/cpp.md), [Web UI](usage/web.md), and [Notebooks](usage/notebooks.md).
- **[Dev](dev/architecture.md)** — architecture, [building](dev/building.md), [contributing](dev/contributing.md), [Postgres dev-env](dev/postgres.md), [web design spec](dev/regpoly-web-design-spec.md), [user stories](dev/web-user-stories.md).

## Layout

REGPOLY is a `uv` workspace of four packages:

| Package | Role |
|---|---|
| `regpoly-cpp` | C++ core: GF(2) algebra, lattice methods, all generator families, search drivers, YAML, catalog. Publishes the `_regpoly_cpp` Python extension and a standalone `regpoly-cli` binary. Can also be built **without** Python/pybind11 (`-DREGPOLY_BUILD_PYTHON_EXTENSION=OFF`). |
| `regpoly` | Pure-Python wrapper around `regpoly_cpp`, plus result decoration (DataFrames, plots) and the Python `regpoly` CLI. |
| `regpoly-web` | FastAPI web UI for launching searches and browsing results (PostgreSQL via `psycopg3`). Uses `regpoly` only. |
| `regpoly-legacy` | Optional pure-Python add-on that reads the pre-v2 `.dat` parameter format and builds `Generator` objects through the same C++ factory. |

Dependency arrows are strictly one-way:

- `regpoly-web → regpoly → regpoly-cpp`
- `regpoly-legacy → regpoly → regpoly-cpp`

Both layered contracts are enforced by `import-linter` in CI.
