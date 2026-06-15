# regpoly — C++ core

Self-contained C++ library + CLI for analysing and searching combined PRNGs
based on GF(2) linear recurrences. **No Python required.** The Python wrapper
lives in `../python/`; the catalog data in `../data/`.

## What you get

- `libregpoly_core.a` — the static library (GF(2) algebra, MT-family generators,
  lattice/equidistribution methods, catalog + YAML config readers).
- `regpoly-cli` — command-line tool (`catalog`, `search`, `show`, `publish`).
- `find_package(regpoly)` support: the imported target `regpoly::regpoly`, the
  umbrella header `<regpoly/regpoly.h>`, and a CMake config that pulls in NTL and
  yaml-cpp transitively.

## Dependencies

- CMake ≥ 3.20, a C++17 compiler.
- **NTL** (+ GMP) — `apt install libntl-dev` (pulls GMP). Located via the bundled
  `cmake/FindNTL.cmake` (honours `NTL_ROOT` / `CMAKE_PREFIX_PATH`).
- **yaml-cpp ≥ 0.8** — used if found (`apt install libyaml-cpp-dev`); otherwise
  fetched and installed automatically (needs network at configure time; set
  `FETCHCONTENT_SOURCE_DIR_YAML-CPP` for offline builds).

## Build, test, install

```bash
# from the repo root
cmake --preset cpp-only          # = -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF -DREGPOLY_BUILD_TESTS=ON
cmake --build --preset cpp-only -j
ctest --preset cpp-only          # GoogleTest suite

cmake --install build-cpp --prefix /opt/regpoly
# installs: lib/libregpoly_core.a, bin/regpoly-cli, include/regpoly/*.h,
#           lib/cmake/regpoly/* (config + FindNTL), share/regpoly/catalog + papers
```

`regpoly-cli` finds its catalog automatically after install (it looks at
`<bindir>/../share/regpoly/catalog`); override with `--library DIR` or
`$REGPOLY_CATALOG_DIR`.

## Use it from your own CMake project

```cmake
find_package(regpoly REQUIRED)          # add /opt/regpoly to CMAKE_PREFIX_PATH
add_executable(app main.cpp)
target_link_libraries(app PRIVATE regpoly::regpoly)
```
```cpp
#include <regpoly/regpoly.h>
#include <regpoly/catalog.h>
```

A complete, build-tested consumer is in [`examples/`](examples/) — it's the CI
acceptance gate (configured against the install prefix, no Python involved).

## Key CMake options

| Option | Default | Meaning |
|---|---|---|
| `REGPOLY_BUILD_PYTHON_EXTENSION` | ON | Build the pybind11 `_regpoly_cpp` module. **OFF** = pure C++. |
| `REGPOLY_BUILD_CLI` | ON | Build `regpoly-cli`. |
| `REGPOLY_BUILD_TESTS` | OFF | Build the GoogleTest suite (`ctest`). |
| `REGPOLY_INSTALL_CMAKE_PACKAGE` | ON | Install `find_package(regpoly)` config. |
| `REGPOLY_NATIVE_ARCH` | ON | `-march=native` (dev). OFF + `REGPOLY_TUNE_ARCH` for portable/tuned builds. |
