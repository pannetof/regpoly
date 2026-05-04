# Using REGPOLY from C++

REGPOLY ships a standalone C++ binary, `regpoly-cli`, plus a CMake-installable library so downstream C++ projects can link against the search drivers, the catalog, and the legacy `.dat` reader without bringing in Python.

## Build

```bash
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build -j
```

The `regpoly-cli` binary is produced in `build/regpoly-cli`. Run the GoogleTest suite via `ctest`:

```bash
ctest --test-dir build --output-on-failure
```

## CLI subcommands

```text
regpoly-cli <command> [options]

  catalog list [--library DIR]
  catalog show PAPER_ID [--library DIR]
  catalog gen GEN_ID [--library DIR]
  legacy-info FILE.dat [-L N]
  legacy-trans FILE.dat
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

A C++-only user has feature parity with a Python user: the same YAML configs, the same catalog, the same legacy `.dat` reader. The Python wrappers add result decoration (DataFrames, plots) but do not gate any algorithm.

## See also

- [Python usage](python.md) — wrapper layer over `regpoly_cpp`.
- [Search format](../theory/search_format.md) — YAML schema accepted by `search` and `regpoly`.
- [Architecture](../dev/architecture.md) — package layout + dependency arrows.
