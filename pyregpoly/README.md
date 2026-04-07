# pyregpoly --- Pure Python implementation of REGPOLY

Pure Python reference implementation for analyzing the equidistribution
properties of pseudo-random number generators based on modulo-2 linear
recurrences.

No compilation or external dependencies required. This is the simplest
implementation, suitable for understanding the algorithms and for development
and testing. It supports only the legacy C-compatible input format.

## Dependencies

- Python 3.10+

No external packages needed.

## Usage

This is not an installable package. Run directly from the directory:

```bash
cd pyregpoly
python main.py <nb_comp> <test_file> <gen_file1> [gen_file2 ...]
```

- `nb_comp` --- number of generator components (1 or 2)
- `test_file` --- test configuration file (legacy format)
- `gen_file_i` --- generator data file(s), one per component

### Examples

Run from the `pyregpoly/` directory, pointing to data files in
`legacy_parameters/`:

**Single TGFSR generator:**

```bash
python main.py 1 ../legacy_parameters/example1c ../legacy_parameters/96_1.dt
```

**Single TGFSR generator (k=155):**

```bash
python main.py 1 ../legacy_parameters/example1c ../legacy_parameters/155_2.dt
```

**Two-component combined generator:**

```bash
python main.py 2 ../legacy_parameters/example3 ../legacy_parameters/96_1.dt ../legacy_parameters/155_2.dt
```

## Supported Generator Families

- `tgfsr` --- Twisted GFSR
- `polylcg` --- Polynomial LCG
- `tausworthe` --- Tausworthe

## Running Tests

```bash
pip install pytest
pytest
```

## Source Organization

```
pyregpoly/
  main.py              command-line entry point
  seek.py              search orchestration
  generateur.py        base generator class
  bitvect.py           bit vector implementation
  matrix.py            Gaussian elimination
  combinaison.py       generator combinations
  polylcg.py           Polynomial LCG generator
  tausworthe.py        Tausworthe generator
  tgfsr.py             Twisted GFSR generator
  temper_mk.py         Matsumoto-Kurita tempering
  permutation.py       permutation tempering
  equidistribution_test.py   ME test
  collision_free_test.py     CF test
  test_*.py            unit tests
```
