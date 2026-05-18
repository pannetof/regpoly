# Using REGPOLY from C++

REGPOLY ships a standalone C++ binary, `regpoly-cli`, plus a CMake-installable library so downstream C++ projects can link against the search drivers and the catalog without bringing in Python.

The pre-v2 `.dat` parameter format is **not** parseable from C++; reading those files requires the optional pure-Python `regpoly-legacy` add-on package (`uv run regpoly-legacy info FILE.dat`).

## Build

```bash
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build -j
```

The `regpoly-cli` binary is produced in `build/regpoly-cli`. Run the GoogleTest suite via `ctest`:

```bash
ctest --test-dir build --output-on-failure
```

### Pure-C++ build (no Python, no pybind11)

The default configure builds the `_regpoly_cpp` pybind11 extension alongside the static library and CLI, which requires `pybind11` to be `find_package`-able. C++-only consumers who don't want any Python tooling pass `-DREGPOLY_BUILD_PYTHON_EXTENSION=OFF`:

```bash
cmake -S packages/regpoly-cpp -B build \
      -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF \
      -DREGPOLY_BUILD_TESTS=ON
cmake --build build -j
```

This produces `libregpoly_core.a`, `regpoly-cli`, the public headers, and the `find_package(regpoly)` CMake config package — no `.so`, no pybind11 dependency. The GoogleTest suite (`ctest`) still runs in this mode.

## CLI subcommands

```text
regpoly-cli <command> [options]

  catalog list [--library DIR]
  catalog show PAPER_ID [--library DIR]
  catalog gen GEN_ID [--library DIR]
  search FILE.yaml
  show FILE.yaml
  publish FILE.yaml --paper PAPER_ID --gen-id GEN_ID
          [--display TEXT] [--target STR] [--starred] [--library DIR]
```

`search FILE.yaml` runs an equidistribution search against the same seek-style YAML format the Python CLI accepts (see [Search format](../theory/search_format.md)). `show FILE.yaml` prints a tested-generator YAML — components, tempering chain, and analysis results. `publish` appends a tested-generator entry to an existing paper YAML under the catalog and is the canonical way to grow `docs/library/` from C++.

## Linking from another CMake project

A future install step will publish a `regpoly-config.cmake` package so downstream consumers can write:

```cmake
find_package(regpoly REQUIRED)
target_link_libraries(my_project PRIVATE regpoly::regpoly)
```

Until then, link against the `regpoly_core` static library and add the public headers from `packages/regpoly-cpp/src/include/` to your include path.

## C++/Python parity

A C++-only user has feature parity with a Python user for the **catalog- and YAML-driven** path: the same YAML configs, the same catalog, the same search drivers. The Python wrappers add result decoration (DataFrames, plots) but do not gate any algorithm.

The one piece of feature drift is the pre-v2 `.dat` parameter format. That reader was retired from the C++ core in v2.x — it now lives only in the pure-Python `regpoly-legacy` add-on, so a C++-only user must either convert their `.dat` fixture to the YAML catalog format up-front, or fall back to invoking `regpoly-legacy seek` from the command line.

## See also

- [Python usage](python.md) — wrapper layer over `regpoly_cpp`.
- [Search format](../theory/search_format.md) — YAML schema accepted by `search` and `regpoly`.
- [Architecture](../dev/architecture.md) — package layout + dependency arrows.
