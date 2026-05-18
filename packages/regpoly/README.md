# regpoly

Pure-Python wrapper layer over the C++ core (`regpoly-cpp`). Provides the
high-level Python API, the `regpoly` console script, YAML config plumbing,
result decoration (DataFrames, plots), and the per-family catalog loader.

Depends on `regpoly-cpp`. Installed automatically as part of the workspace
`uv sync`.

## Run a search

```bash
uv run regpoly shared/yaml/equidist/mt19937.yaml
```

YAML search configs live under `shared/yaml/` at the workspace root.
See [docs/theory/search_format.md](../../docs/theory/search_format.md)
for the schema.

## Programmatic use

See [docs/usage/python.md](../../docs/usage/python.md) for the canonical
Python entry points (`regpoly.search.seek.Seek`, `regpoly.library.Catalog`,
`regpoly.analyses.pis.analyze_single_generator`, …).

## Source layout

```
src/regpoly/
├── analyses/         result decoration on top of regpoly_cpp analyses
├── core/             Generator, Combination, family aliases
├── data/             primitive-factor JSON + other static data
├── io/               YAML + CLI plumbing
├── library/          Catalog wrapper
├── search/           Seek, PrimitiveSearch, TemperingSearch Python shims
└── tools/            dev-time codegen (gen_primitive_factors_data, …)
```

## Tests

```bash
uv run pytest packages/regpoly/tests
uv run pytest -m slow packages/regpoly/tests   # MTToolBox cross-checks
```

See [`tests/README.md`](tests/README.md) for the slow-lane cross-check
workflow and [`tests/_mttoolbox_build.md`](tests/_mttoolbox_build.md) for
building the upstream `calc_equidist` reference binaries.
