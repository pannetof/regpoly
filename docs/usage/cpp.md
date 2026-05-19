# Using REGPOLY from C++

REGPOLY ships a standalone C++ binary, `regpoly-cli`, plus a
CMake-installable library so downstream C++ projects can link against
the search drivers and the catalog without bringing in Python.

The pre-v2 `.dat` parameter format is **not** parseable from C++;
reading those files requires the optional pure-Python `regpoly-legacy`
add-on (`uv run regpoly-legacy info FILE.dat`).

## Build

```bash
cmake -S packages/regpoly-cpp -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The `regpoly-cli` binary lands in `build/regpoly-cli`.

### Pure-C++ build (no Python, no pybind11)

The default configure builds the `_regpoly_cpp` pybind11 extension
alongside the static library and CLI, which requires pybind11 to be
`find_package`-able. C++-only consumers who don't want any Python
tooling pass `-DREGPOLY_BUILD_PYTHON_EXTENSION=OFF`:

```bash
cmake -S packages/regpoly-cpp -B build-cpp \
      -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF \
      -DREGPOLY_BUILD_TESTS=ON \
      -DCMAKE_INSTALL_PREFIX=$PWD/build-cpp/install
cmake --build build-cpp --target install -j
```

This produces `libregpoly_core.a`, `regpoly-cli`, the public headers
under `${prefix}/include/regpoly/`, and the `find_package(regpoly)`
CMake config package — no `.so`, no pybind11 dependency. The
GoogleTest suite still runs via `ctest`.

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

`search FILE.yaml` runs an equidistribution search against the same
seek-style YAML format the Python CLI accepts (see
[Search format](../theory/search_format.md)). `show FILE.yaml` prints
a tested-generator YAML — components, tempering chain, and analysis
results. `publish` appends a tested-generator entry to an existing
paper YAML under the catalog and is the canonical way to grow
`docs/library/` from C++.

## Linking from a downstream CMake project

After the pure-C++ install above, a consumer project finds REGPOLY via
the standard `find_package` mechanism:

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_app CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# yaml-cpp is a private link dependency of regpoly_core. Bring it in
# either via apt (libyaml-cpp-dev) or by vendoring the same submodule
# REGPOLY uses (third_party/yaml-cpp). For the vendored path:
add_subdirectory(third_party/yaml-cpp yaml-cpp_build EXCLUDE_FROM_ALL)

find_package(regpoly REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE
    regpoly::regpoly_core
    yaml-cpp::yaml-cpp)
```

Configure your consumer with:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/regpoly/build-cpp/install \
              -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
cmake --build build -j
```

The public C++ surface is the umbrella header
`<regpoly/regpoly.h>` (abstract `Generator`, `Transformation`,
`Combination`, `CombinedGenerator`, `BitVect`, the lattice / equidist
runners) plus the catalog (`<regpoly/catalog.h>`) and the YAML config
loader (`<regpoly/seek_config.h>`). Concrete generator families are
not part of the public ABI — build them via the catalog or the search
config loader rather than directly.

## "Hello equidistribution" — a runnable example

The snippet below loads a workspace YAML config, builds the combined
generator it describes, runs the matricial equidistribution analysis,
and prints the total dimension defect (SE). It is exactly what
{func}`regpoly.analyses.pis.analyze_single_generator` does on the
Python side, in ~40 lines of straight C++.

`main.cpp`:

```cpp
// Build and run matricial equidistribution on whatever generator the
// supplied YAML search config describes. Print SE.
#include <regpoly/regpoly.h>
#include <regpoly/seek_config.h>
#include <regpoly/combined.h>
#include <regpoly/equidistribution_runner.h>
#include <climits>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
    using namespace regpoly::core;
    namespace yc = regpoly::yaml_config;

    const std::string yaml_path =
        argc > 1 ? argv[1] : "shared/yaml/equidist/mt19937.yaml";
    constexpr int L = 32;

    // 1. Load the YAML search config and build the runtime Combination.
    auto cfg   = yc::load_seek_config(yaml_path);
    auto built = yc::build_search(cfg);
    Combination& comb = *built.combination;

    // 2. Snapshot the active combination — materialise the tempering
    //    chain into a single CombinedGenerator we can iterate.
    auto combined = build_combined_from_combination(comb);

    // 3. Matricial equidistribution. delta + mse are capped so nothing
    //    is rejected — we just want the score.
    std::vector<int> delta(L + 1, INT_MAX);
    auto res = run_matricial_equidistribution(
        *combined, comb.k_g(), L, /*Lmax=*/L, delta, /*mse=*/INT_MAX);

    std::cout << "Equidistribution from " << yaml_path << '\n'
              << "  k_g      = " << comb.k_g() << '\n'
              << "  L        = " << L << '\n'
              << "  SE       = " << res.se << '\n'
              << "  verified = " << (res.verified ? "yes" : "no") << '\n';
    return 0;
}
```

Build + run:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/regpoly/build-cpp/install \
              -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
cmake --build build -j
./build/my_app /path/to/regpoly/shared/yaml/equidist/mt19937.yaml
```

Expected output (the SE is the *initial* score for the YAML's
specified tempering; a real search would minimise it):

```text
Equidistribution from .../mt19937.yaml
  k_g      = 19937
  L        = 32
  SE       = 6750
  verified = yes
```

## C++/Python parity

A C++-only user has feature parity with a Python user for the
**catalog- and YAML-driven** path: same YAML configs, same catalog,
same search drivers. The Python wrappers add result decoration
(DataFrames, plots) but do not gate any algorithm.

One piece of feature drift: the pre-v2 `.dat` parameter format was
retired from the C++ core in v2.x — it now lives only in the
pure-Python `regpoly-legacy` add-on. A C++-only user must either
convert their `.dat` fixture to the YAML catalog format up-front, or
fall back to `uv run regpoly-legacy seek` from the command line.

See [The Python / C++ bridge](python-cpp-bridge.md) for the layered
view of which concerns live where, including the wrapper-class
cross-reference table.

## See also

- [Python usage](python.md) — wrapper layer over `regpoly_cpp`.
- [The Python / C++ bridge](python-cpp-bridge.md) — what each layer
  owns + the wrapper-to-C++ map.
- [Search format](../theory/search_format.md) — YAML schema accepted by
  `search` / `regpoly` / the snippet above.
- [Architecture](../dev/architecture.md) — package layout + dependency
  arrows.
