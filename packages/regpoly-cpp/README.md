# regpoly-cpp

C++ core of REGPOLY: GF(2) linear algebra, F₂-linear generator families,
equidistribution analyses (matricial + lattice), YAML search-config parser,
catalog reader, the `regpoly-cli` standalone binary, and an optional
pybind11 extension (`_regpoly_cpp`) consumed by the `regpoly` Python wrapper.

## Two build modes

**Python wheel (default)** — built by `uv sync` from the workspace root via
scikit-build-core. Produces `libregpoly_core.a`, `regpoly-cli`, *and*
the `_regpoly_cpp` pybind11 extension.

```bash
cd ../..
uv sync
```

**Pure C++ (no Python, no pybind11)** — plain CMake. Produces the static
library, the CLI, public headers, and a `find_package(regpoly)` config
package, with no Python dependency at configure/build/install time.

```bash
cmake -S . -B build -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF
cmake --build build -j
sudo cmake --install build
```

## Tests

```bash
cmake -S . -B build -DREGPOLY_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Python tests for the pybind11 binding ABI live under `tests/python/`
and are picked up by the workspace `pytest`.

## Source layout

```
src/
├── algebra/          GF(2) polynomial / matrix primitives, primitivity
├── analyses/         equidistribution + tuple-test cores
├── bindings/         pybind11 module (Python wheel mode only)
├── cli/              regpoly-cli entry point
├── core/             Generator, CombinedGenerator, factory, parameter specs
├── generators/       per-family implementations (one .cpp per family)
├── include/          public headers (mirrors src/ subdir layout)
├── io/               catalog reader/writer, YAML helpers
├── lattice/          dual-lattice equidistribution machinery
├── library/          Catalog
├── search/           PrimitiveSearch, Seek, TemperingSearch, TemperingOptimizer
├── transforms/       tempering transformations (TempMK, etc.)
└── yaml_config/      YAML search-config parser
```

## Documentation

- [docs/dev/architecture.md](../../docs/dev/architecture.md) — layered package architecture.
- [docs/dev/building.md](../../docs/dev/building.md) — full build recipes.
- [docs/usage/cpp.md](../../docs/usage/cpp.md) — using `regpoly-cli` and linking against the installed CMake package.
- [docs/generators/](../../docs/generators/) — one page per generator family.
