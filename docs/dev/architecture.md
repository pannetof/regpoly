# Architecture

REGPOLY is a `uv` workspace of three packages with strict one-way dependencies:

```
regpoly-web      →  regpoly  →  regpoly-cpp
   FastAPI           Python      C++ core
```

`import-linter` enforces the layered arrow in CI. The C++ core carries every algorithm, the YAML search schema, the catalog reader/writer, and a standalone `regpoly-cli` binary so a user working exclusively in C++ has feature parity with a Python user for the *catalog-driven* path. The historical pre-v2 `.dat` parameter format is **no longer supported** — neither from C++ nor from Python.

## Build modes

`packages/regpoly-cpp/` supports two build modes:

- **Python wheel** (default, via scikit-build-core + `uv sync`): builds `libregpoly_core.a`, the `regpoly-cli` binary, **and** the `_regpoly_cpp` pybind11 extension. Required by every Python consumer (`regpoly`, `regpoly-web`).
- **Pure C++** (via `cmake -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF`): builds the static library + CLI + headers + `find_package(regpoly)` config package, with **no pybind11 or Python dependency at configure, build, or install time**. The pybind11 extension is skipped entirely.

Both modes share the same `regpoly_core` static library — the Python wheel just adds a thin pybind11 module on top.

## See also

- [Building](building.md) — full build / test recipes (Python + pure C++ + docs).
- [Contributing](contributing.md) — how to add a generator family + CI gates.
- [PostgreSQL for local dev](postgres.md) — spinning up a dev DB for `regpoly-web`.
- [Web — design spec](regpoly-web-design-spec.md) — the FastAPI app's internal design.
- [Web — user stories](web-user-stories.md) — UX scenarios driving the web UI.
