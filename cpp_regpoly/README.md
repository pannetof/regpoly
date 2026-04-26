# cpp_regpoly --- C++/Python implementation of REGPOLY

C++ accelerated Python package for analyzing the equidistribution properties of
pseudo-random number generators based on modulo-2 linear recurrences.
Uses pybind11 for Python bindings and NTL for polynomial arithmetic.

This is the primary implementation, supporting both YAML configuration files
and the legacy C-compatible input format.

## Dependencies

- Python 3.10+
- pybind11 (installed automatically)
- PyYAML (installed automatically)
- **libntl-dev** (NTL number theory library, must be installed system-wide)

Install NTL on Debian/Ubuntu:

```bash
sudo apt install libntl-dev
```

## Installation

```bash
cd cpp_regpoly
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Usage

### YAML mode

```bash
python -m regpoly.io.cli <config.yaml>
```

YAML configuration files are in the `yaml/` directory. They specify the
generator family, parameters, tempering, and test method in a single file.

**Single generator, lattice method (WELL19937a):**

```bash
python -m regpoly.io.cli ../yaml/equidist/well19937a.yml
```

**Mersenne Twister MT19937 with tempering:**

```bash
python -m regpoly.io.cli ../yaml/equidist/mt19937.yaml
```

**Two-component combined generator:**

```bash
python -m regpoly.io.cli ../yaml/equidist/example3.yaml
```

### Legacy mode

```bash
python -m regpoly.io.cli <nb_comp> <test_file> <gen_file1> [gen_file2 ...]
```

This mode reads the same input files as the C `POL` program.

**Single TGFSR generator:**

```bash
cd ../legacy_parameters
python -m regpoly.io.cli 1 example1c 96_1.dt
```

**Two-component combined generator:**

```bash
cd ../legacy_parameters
python -m regpoly.io.cli 2 example3 96_1.dt 155_2.dt
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
- `MatsumotoGen` --- Matsumoto generators
- `AC1DGen` --- AC1D generators
- `MELGGen` --- Maximally Equidistributed F₂-Linear (Harase-Kimoto)
- `SFMTGen`, `DSFMTGen` --- SIMD-Friendly MT (single / double precision)
- `MTGPGen` --- MT for Graphic Processors
- `TinyMT32Gen` --- Tiny MT (127-bit state)
- `RMT64Gen` --- Reducible MT (64-bit)
- `XorShift128Gen` --- Marsaglia 128-bit Xor-shift

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
      io/                 cli, legacy_reader, tested_generator
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
