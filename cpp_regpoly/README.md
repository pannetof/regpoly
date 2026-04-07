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
python -m regpoly.cli <config.yaml>
```

YAML configuration files are in the `yaml/` directory. They specify the
generator family, parameters, tempering, and test method in a single file.

**Single generator, lattice method (WELL19937a):**

```bash
python -m regpoly.cli ../yaml/search.well19937a.yaml
```

**Mersenne Twister MT19937 with tempering:**

```bash
python -m regpoly.cli ../yaml/search.mt19937.yaml
```

**Two-component combined generator:**

```bash
python -m regpoly.cli ../yaml/search.example3.yaml
```

### Legacy mode

```bash
python -m regpoly.cli <nb_comp> <test_file> <gen_file1> [gen_file2 ...]
```

This mode reads the same input files as the C `POL` program.

**Single TGFSR generator:**

```bash
cd ../legacy_parameters
python -m regpoly.cli 1 example1c 96_1.dt
```

**Two-component combined generator:**

```bash
cd ../legacy_parameters
python -m regpoly.cli 2 example3 96_1.dt 155_2.dt
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
      family: MT
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

- `tgfsr` --- Twisted GFSR
- `MT` --- Mersenne Twister
- `polylcg` --- Polynomial LCG
- `tausworthe` --- Tausworthe
- `carry` --- Carry / WELL generators
- `genf2w` --- GF(2^w) generators (polynomial LCG and LFSR variants)
- `marsaxorshift` --- Marsaglia Xor-shift
- `matsumoto` --- Matsumoto generators
- `AC1D` --- AC1D generators

## Source Organization

```
cpp_regpoly/
  src/
    cpp/
      include/       C++ headers
      src/           C++ source files (generators, gauss, lattice, bindings)
    regpoly/
      cli.py         command-line entry point
      seek.py        search orchestration (YAML + legacy)
      generateur.py  Python generator wrapper
      tests/         equidistribution, collision-free, tuplets tests
  setup.py           pybind11 build configuration
  pyproject.toml     package metadata
```
