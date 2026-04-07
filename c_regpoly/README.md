# c_regpoly --- C implementation of REGPOLY

C library and command-line tool (`POL`) for analyzing the equidistribution
properties of pseudo-random number generators based on modulo-2 linear
recurrences (LFSRs over GF(2)).

## Dependencies

- **ZEN** --- finite-field arithmetic library (included separately in the project)
- **TestU01** --- statistical testing library (optional, required only for distance tests)
- GCC with C99 support

## Building

Edit the `Makefile` to set the paths for ZEN and TestU01 to match your installation:

```makefile
DIRZEN      = /path/to/ZEN
DIRTESTU01LIB = /path/to/testu01/install/lib
DIRTESTU01INC = /path/to/testu01/install/include
```

Then build:

```bash
cd c_regpoly
make all        # builds the POL executable
make clean      # removes all build artifacts
```

## Usage

```
./POL <nb_components> <test_config> <gen_file_1> [gen_file_2 ...]
```

- `nb_components` --- number of generator components in the combination (1 or 2)
- `test_config` --- test configuration file specifying the method and thresholds
- `gen_file_i` --- generator data file(s) for each component

Generator and configuration files are expected in the working directory
(typically `legacy_parameters/`).

### Examples

Run from the repository root (`MinimalCode/`), where the `legacy_parameters/`
directory contains the data files:

```bash
cd legacy_parameters
```

**Single TGFSR generator, matricial method:**

```bash
../c_regpoly/POL 1 example1c 96_1.dt
```

**Single TGFSR generator (k=155), matricial method:**

```bash
../c_regpoly/POL 1 example1c 155_2.dt
```

**Two-component combined generator:**

```bash
../c_regpoly/POL 2 example3 96_1.dt 155_2.dt
```

**WELL19937a generator, lattice method:**

```bash
../c_regpoly/POL 1 example_well44497b carry32_624_31_final.dat
```

**WELL44497b generator, lattice method:**

```bash
../c_regpoly/POL 1 example_well44497b carry32_1391_15_final.dat
```

## Test Configuration File Format

Each line of the configuration file specifies one parameter:

```
1059674621  1170794524     # Seeds for tempering parameter RNG
32                         # Lmax (maximum resolution)
0 trans32.dat              # Tempering: 0=disabled, 1=enabled + tempering file
1                          # Number of tempering tries per generator
10000 matrix               # Max sum of dimension gaps + method (matrix|lattice)
0                          # Per-resolution gap restriction (0=none)
0                          # Non-successive output test (0=disabled)
```

The method can be `matrix` (Gaussian elimination) or `lattice`
(dual lattice basis reduction with Lenstra's algorithm).

## Generator Data File Format

Each generator type has its own file format. The first line is the type tag.

**TGFSR** (`tgfsr`):

```
tgfsr
32 3            # w r
5               # number of generators
1 9965523b      # m a  (one line per generator)
1 d0ca64b1
...
```

**Mersenne Twister** (`MT`):

```
MT
1               # number of generators
32 624 397 31 9908b0df    # w r m p a
```

**Carry / WELL** (`carry`):

```
carry
32 624 31       # w r p
1               # number of generators
70 179 449  0 -25 0 0 00000000 00000000 00000000  ...  # m1 m2 m3 + 8 matrices
```

Each matrix entry is: `type pi0 pi1 pi2 pu0 pu1 pu2` (type + 3 integer params
+ 3 hex params).

## Output

POL prints for each generator:
- The characteristic polynomial hamming weight
- Generator parameters
- Dimension gaps table for resolutions l = 1 to Lmax
- The sum of dimension gaps over Psi_12

Example output:

```
hammingweigth = 8585
Carry Generator:
19937
 w=  32  r=624  p=  31  m1= 70  m2=179  m3=449

  Dimension gaps for every resolution
=======+=====+=====+=====+=====+
RESOL  |    1|    2|    3|    4|
-------+-----+-----+-----+-----|
ECART  |     |     |     |    1|
-------+-----+-----+-----+-----|
DIM    |19937| 9968| 6645| 4983|
=======+=====+=====+=====+=====+
```

A generator with all zeros in the ECART row is maximally equidistributed (ME)
at those resolutions.

## Source Organization

```
c_regpoly/
  src/           C source files
  include/       C header files
  lib/           compiled objects and libREGPOLY.a
  Makefile
```
