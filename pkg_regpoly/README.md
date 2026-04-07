# pkg_regpoly --- Cython-accelerated Python implementation of REGPOLY

Python package for analyzing the equidistribution properties of pseudo-random
number generators based on modulo-2 linear recurrences.
Performance-critical modules are optionally compiled with Cython.

This package has no external dependencies beyond Python itself.
It supports only the legacy C-compatible input format (no YAML).

## Dependencies

- Python 3.10+
- Cython 3.0+ (optional, for compilation; falls back to pure Python)

## Installation

```bash
cd pkg_regpoly
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

To install with Cython acceleration:

```bash
pip install cython>=3.0
pip install -e .
```

## Usage

```bash
python -m regpoly.cli <nb_comp> <test_file> <gen_file1> [gen_file2 ...]
```

- `nb_comp` --- number of generator components (1 or 2)
- `test_file` --- test configuration file (legacy format)
- `gen_file_i` --- generator data file(s), one per component

### Examples

Run from the `legacy_parameters/` directory where the data files are located:

```bash
cd ../legacy_parameters
```

**Single TGFSR generator:**

```bash
python -m regpoly.cli 1 example1c 96_1.dt
```

**Single TGFSR generator (k=155):**

```bash
python -m regpoly.cli 1 example1c 155_2.dt
```

**Two-component combined generator:**

```bash
python -m regpoly.cli 2 example3 96_1.dt 155_2.dt
```

## Supported Generator Families

- `tgfsr` --- Twisted GFSR
- `polylcg` --- Polynomial LCG
- `tausworthe` --- Tausworthe
- `genf2w` --- GF(2^w) generators
- `MT` --- Mersenne Twister

## Cython Modules

The following modules are compiled with Cython when available:

- `bitvect` --- bit vector operations
- `gauss_matrix` --- Gaussian elimination
- `generateur` --- base generator class
- `polylcg`, `tausworthe`, `tgfsr`, `genf2w` --- generator implementations

Without Cython, all modules run as pure Python (slower but fully functional).

## Source Organization

```
pkg_regpoly/
  src/regpoly/
    cli.py             command-line entry point
    seek.py            search orchestration
    generateur.py      base generator class
    bitvect.py         bit vector implementation
    gauss_matrix.py    Gaussian elimination
    generators/        generator family implementations
    tests/             equidistribution and collision-free tests
  setup.py             Cython build configuration
  pyproject.toml       package metadata
```
