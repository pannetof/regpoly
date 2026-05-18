# regpoly --- C++/Python toolkit for combined PRNGs

Analyse the equidistribution properties of pseudo-random number
generators based on modulo-2 linear recurrences (LFSRs over GF(2)).
Uses NTL for polynomial arithmetic, with optional pybind11 bindings for
the Python wrapper layer.

The primary use-case is the full Python stack (`regpoly`, `regpoly-web`,
optional `regpoly-legacy` add-on) built on top of the C++ core. The
C++ core can also be built and installed **standalone**, without Python
or pybind11, for pure C++ consumers.

Search configuration uses YAML; the pre-v2 C-compatible `.dat`
parameter format is supported via the optional pure-Python
`regpoly-legacy` add-on package.

## Dependencies

**Required (always):**
- **libntl-dev** (NTL number theory library, must be installed system-wide)

**Required for the Python wheel build (the default `uv sync` path):**
- Python 3.10+
- pybind11 ≥ 2.12 (pulled in automatically by scikit-build-core)
- PyYAML (pulled in automatically)

**For pure-C++ builds:** none of the Python-side dependencies above are needed; see [docs/usage/cpp.md](docs/usage/cpp.md#pure-c-build-no-python-no-pybind11).

Install NTL on Debian/Ubuntu:

```bash
sudo apt install libntl-dev
```

## Installation

### Python (Python wheel, default)

```bash
cd regpoly_monorepo
uv sync           # builds the wheel via scikit-build-core
```

### Pure C++ (no Python, no pybind11)

```bash
cd regpoly_monorepo
cmake -S packages/regpoly-cpp -B build \
      -DREGPOLY_BUILD_PYTHON_EXTENSION=OFF
cmake --build build -j
sudo cmake --install build
```

Produces `libregpoly_core.a`, the `regpoly-cli` binary, public headers, and the `find_package(regpoly)` CMake config package.

## Usage

### YAML mode

```bash
python -m regpoly.io.cli <config.yaml>
```

YAML configuration files are in the `yaml/` directory. They specify the
generator family, parameters, tempering, and test method in a single file.

**Single generator, lattice method (WELL19937a):**

```bash
python -m regpoly.io.cli yaml/equidist/well19937a.yml
```

**Mersenne Twister MT19937 with tempering:**

```bash
python -m regpoly.io.cli yaml/equidist/mt19937.yaml
```

**Two-component combined generator:**

```bash
python -m regpoly.io.cli yaml/equidist/example3.yaml
```

### Legacy mode

```bash
python -m regpoly.io.cli <nb_comp> <test_file> <gen_file1> [gen_file2 ...]
```

This mode reads the same input files as the C `POL` program and now
lives in the optional `regpoly-legacy` add-on package (installed by
default in the workspace). The fixtures themselves live under
`packages/regpoly-legacy/shared/legacy_parameters/`.

**Single TGFSR generator:**

```bash
cd packages/regpoly-legacy/shared/legacy_parameters
uv run regpoly-legacy seek 1 example1c 96_1.dt
```

**Two-component combined generator:**

```bash
cd packages/regpoly-legacy/shared/legacy_parameters
uv run regpoly-legacy seek 2 example3 96_1.dt 155_2.dt
```

## YAML Configuration Format

A search configuration file ties together generators, tempering, and tests:

```yaml
search:
  seed: [1059674621, 1170794524]
  Lmax: 32

components:
  - generators:
      file: gen.well19937.yaml        # load from file
    tempering:
      - type: tempMK2
        w: 32
        eta: 7
        mu: 15
        u: 11
        l: 18
        b: 0x9d2c5680
        c: 0xefc60000

tests:
  - type: equidistribution
    method: lattice                    # or "matricial"
    max_gap_sum: 100000
```

Generators can also be defined inline:

```yaml
components:
  - generators:
      family: MTGen
      common:
        w: 32
        r: 624
        m: 397
        p: 31
      generators:
        - a: 0x9908b0df
```

## Equidistribution Methods

- **matricial** --- Gaussian elimination on the generator matrix. Fast for
  moderate state sizes.
- **lattice** --- Dual lattice basis reduction using Lenstra's algorithm with
  NTL for polynomial arithmetic. Required for large generators (k > 10000).

## Supported Generator Families

All families are exposed under their new canonical class name (with the
`-Gen` suffix). Legacy short tags (lowercase, French) and pre-rename
class names are still accepted via `_FAMILY_ALIASES` and the
backwards-compat aliases in `factory.cpp` / `bindings.cpp`.

- `TGFSRGen` --- Twisted GFSR
- `MTGen` --- Mersenne Twister
- `PolyLCGGen` --- Polynomial LCG
- `TauswortheGen` --- Tausworthe / LFSR
- `WELLGen` --- WELL / carry generators
- `F2wLFSRGen`, `F2wPolyLCGGen` --- GF(2^w) generators
- `MarsaXorshiftGen` --- Marsaglia Xor-shift
- `MELGGen` --- Maximally Equidistributed F₂-Linear (Harase-Kimoto)
- `SFMTGen`, `DSFMTGen` --- SIMD-Friendly MT (single / double precision)
- `MTGPGen` --- MT for Graphic Processors
- `TinyMT32Gen` --- Tiny MT (127-bit state)
- `RMT64Gen` --- Reducible MT (64-bit)

## Source Organization

```
cpp_regpoly/
  src/
    cpp/                  C++ implementation
      include/<bucket>/   headers, grouped by responsibility
      core/               base classes (generator, params, factory, enumerator)
      algebra/            BitVect, GaussMatrix, GF(2^w), Berlekamp-Massey
      generators/         18 generator families (one .cpp per family)
      transforms/         output tempering: permutation, temper_mk, lag_mask
      lattice/            ME analysis: dual_lattice, me_harase, me_notprimitive,
                          me_notprimitive_simd, me_helpers, temper_optimizer
      bindings/           pybind11 module
    python/               Python package (imported as `regpoly`)
      core/               wrappers: generator, combination, transformation, ...
      search/             seek, search_primitive, tempering_search, ...
      io/                 cli, tested_generator
      (pre-v2 .dat reader lives in the optional regpoly-legacy add-on:
       packages/regpoly-legacy/src/regpoly_legacy/reader.py)
      analyses/           equidistribution / collision-free / tuplets tests
      library/            published-generators paper catalog
      web/                FastAPI app + routes + templates + tasks
  docs/
    generators/           per-family markdown (one file per *Gen class)
    library/              paper YAMLs + cross-check fixture YAMLs
  papers/                 reference PDFs for paper YAMLs
  tests/                  pytest suite (smoke, library, MTToolBox crosscheck, …)
  setup.py                pybind11 build configuration
  pyproject.toml          package metadata (package_dir maps regpoly → src/python)
```
