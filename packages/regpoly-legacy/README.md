# regpoly-legacy

Optional pure-Python add-on for reading **pre-v2 `.dat` parameter files** into
`regpoly` `Generator` / `Transformation` objects. Installs as a workspace
member of `regpoly_monorepo`.

This package exists solely to support the historical, frozen `.dat` file
format produced by the pre-v2 C codebase. New code should use the YAML
catalog under `docs/library/` instead.

## What's inside

- `regpoly_legacy.LegacyReader` — main entry point (`read_generators`,
  `read_transformations`).
- `regpoly_legacy.seek_from_legacy(nb_comp, test_file, gen_data_files)` —
  builds a `regpoly.search.seek.Seek` from a positional legacy invocation.
- `regpoly_legacy.tools.transition_matrix` — utility for computing the
  GF(2) transition matrix of a legacy generator.
- `regpoly-legacy` console script with three subcommands:
  - `regpoly-legacy info <file.dat>`  — prints parsed generator specs.
  - `regpoly-legacy trans <file.dat>` — prints parsed transformation specs.
  - `regpoly-legacy seek <nb_comp> <test_file> <gen_file1> [...]` — runs a
    search driven by legacy positional inputs.

## What's *not* inside

No C++. No compiler needed. The .dat parser is pure Python; constructed
`Generator` and `Transformation` objects are built via `regpoly`'s existing
factory, which delegates to the `regpoly_cpp` C++ core.

## Fixtures

Historical fixtures used by the test suite live under `shared/legacy_parameters/`
inside this package. The pre-v2 `shared/legacy_parameters/` directory at the
workspace root has been removed.
