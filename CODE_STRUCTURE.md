# MinimalCode — Code Structure

## Purpose

MinimalCode is a reduced version of the REGPOLY v2.0 library.  Its single
executable, `POL`, searches for high-quality combined pseudo-random number
generators (PRNGs) built from GF(2) linear recurrences.  It tests each
candidate against equidistribution criteria and prints the ones that pass.

Invocation:

```
./POL <NbComp> <test_file> <gen_data_file1> ... <gen_data_file(NbComp)>
```

---

## Directory Layout

```
MinimalCode/
├── Makefile               Build script — compiles all modules into libREGPOLY.a, links POL
├── POL                    Resulting executable
├── src/                   All C source files (19 files)
├── include/               Header files (18 files)
└── lib/                   Object files and libREGPOLY.a (generated)
```

---

## Module Map

The 18 compiled modules (plus `rechercheguide.c` which is the `main`) fall
into five logical layers:

```
┌─────────────────────────────────────────────────────┐
│                   rechercheguide.c  (main / driver)  │
├─────────────────────────────────────────────────────┤
│   readfiles.c     (parse config & data files)        │
├────────────────────────┬────────────────────────────┤
│  mecf.c                │  tuplets.c                  │
│  (equidistribution)    │  (tuple equidistrib.)       │
├────────────────────────┴────────────────────────────┤
│  combinaisons.c    (combined generator framework)    │
│  tempering.c       (output transformation dispatch)  │
├──────────────┬──────────────┬───────────────────────┤
│  Generator   │  Tempering   │  Math / support        │
│  modules     │  modules     │  modules               │
│              │              │                        │
│  tausworthe  │  temperMK    │  vectorsF2             │
│  polylcg     │  permut      │  basisreduc            │
│  tgfsr       │              │  polynomials           │
│  genf2w      │              │  regpoly               │
│  MT          │              │  myzen                 │
│              │              │  timer                 │
└──────────────┴──────────────┴───────────────────────┘
```

---

## File-by-File Description

### Entry Point

#### `rechercheguide.c`
The `main()` function.  Parses command-line arguments, opens the test file,
then calls `Seek()`.  `Seek()` drives the entire search loop:

1. Calls `ReadSearch` and `ReadGenDataFiles` to initialize the `Combinaison`.
2. Iterates over all candidate combined generators with `FirstGen` / `NextGen`.
3. For each candidate, calls `UpdateAllTrans`, then `TestEquid` (or
   `TestMETemperMK` if MK-opt mode).
4. If the ME criterion passes, calls `TestTuplets`.
5. Prints generators that pass both criteria.
6. Prints a summary at the end (total tested, ME count, retained count, CPU time).

---

### Configuration & I/O

#### `readfiles.c` / `readfiles.h`
Parses the two kinds of input files:

- **`ReadSearch`** — reads the main test-parameter file: RNG seed, output
  width `Lmax`, tempering descriptors per component, ME gap bounds, and
  tuple-test parameters.  Calls `SetMRG32k3a` / `MRG32k3a` (from `myzen`)
  to initialize the internal RNG used for random tempering parameter selection.

- **`ReadGenDataFiles`** — reads one data file per component; dispatches to
  the appropriate `ReadDataXXX` function based on a generator-type tag
  (`taus`, `polylcg`, `tgfsr`, `genf2w`, `MT`).

- **`ReadTempering`** — reads a tempering-transformation descriptor file;
  supports `permut` (bit permutation) and `tempMK` / `tempMKopt` /
  `tempMK2` / `tempMK2opt` (Matsumoto-Kurita I and II, optionally with
  optimization).

---

### Generator Framework

#### `combinaisons.c` / `combinaisons.h`
Defines the three core types:

| Type | Description |
|------|-------------|
| `Generateur` | A single PRNG component (function pointers for iteration, init, characteristic polynomial, display) |
| `Component`  | A `Generateur` array plus an array of `Transformation` objects; holds one active generator at a time |
| `Combinaison`| An array of `Component`s; the combined generator is their XOR output |

Key functions:
- `AllocComponentsInCombinaison` / `FreeComponentsInCombinaison` — memory management.
- `FirstGen` / `NextGen` — advance the search through all generator parameter combinations.
- `DispCurrentComb` — prints the parameters of the current combined generator.

#### `tempering.c` / `tempering.h`
Manages the `Transformation` type and applies output transformations:

- `Transform` — applies a sequence of `Transformation` objects to a `BitVect`.
- `UpdateAllTrans` / `UpdateTrans` — randomly re-selects the tempering
  parameters for all components (called once per candidate in the search loop).
- `AllocTransInComponent` / `AddTransInComponent` — build up the transformation
  pipeline for a component.

---

### Equidistribution Analysis

#### `mecf.c` / `mecf.h`
The main equidistribution (ME = Maximally Equidistributed) analysis engine.

- **`paramMecf`** — stores per-resolution gap tables (`ecart`, `Lambda`), set
  indicators (`phi12`, `psi12`), and cumulative gap sums.
- **`SetMethod`** — selects the analysis algorithm:
  - `METHOD_MATRICIAL` — direct matrix method (`TestME`).
  - `METHOD_DUALLATTICE` — dual-lattice method.
  - `METHOD_NOTHING` — skip ME verification.
- **`TestME`** — builds the transition matrix of the combined generator and
  computes `t(l)` (resolution of equidistribution) for each output resolution
  `l` from 1 to `Lmax`.  Uses Gaussian elimination on GF(2) matrices.
- **`PrepareMat` / `PrepareMat2`** — fills a `Matrix` with the output
  bit-vectors of the combined generator over one full period (used by both
  `mecf` and `tuplets`).
- **`ResolutionEquid`** — given a matrix and a set of output indices, computes
  the equidistribution resolution for those indices (uses
  `SpecialGaussianElimination`).
- **`isPresqueME` / `isME`** — predicate functions used by the search loop to
  decide whether to proceed to the tuple test.
- **`DispTable` / `DispME`** — display gap tables and ME status to stdout.

#### `tuplets.c` / `tuplets.h`
Tests the **non-successive** (tuple) equidistribution criterion
Delta(t_1,...,t_d).

- **`paramTuplets`** — stores depth `d`, per-depth target dimensions `tuph[]`,
  computed gaps `gap[]` and `DELTA[]`, and pass/fail thresholds.
- **`TestTuplets`** — reuses `PrepareMat` output and calls `ResolutionEquid`
  for every combination of `dim` non-successive output indices (via the
  internal `nextCombination` iterator).  Only `EQUIDISTRIBUTION` test type is
  supported; other types cause an immediate error exit.
- **`isTuplets`** — returns TRUE if the generator passes the tuple criterion
  (or if tuple testing is disabled).
- **`DispTuplets`** — prints gap tables for successive and non-successive
  dimensions.

---

### Generator Implementations

Each generator module defines `ReadDataXXX`, `AllocXXX`, `InitXXX`,
`IterationXXX`, `PolyCharXXX`, and `DispXXX`, all wired into a `Generateur`
struct via function pointers.

#### `tausworthe.c` / `tausworthe.h`
Tausworthe (LFSR) generator.  The recurrence is a linear feedback shift
register over GF(2).  Parameterized by degree `k`, feedback polynomial, and
step size `s`.  `ReadDataTaus` reads lists of `(k, s, poly)` triples.

#### `polylcg.c` / `polylcg.h`
Polynomial LCG — a linear congruential generator whose multiplier is a
polynomial over GF(2^k).  Parameterized by degree `k`, the modulus polynomial,
and the multiplier polynomial.  `ReadDataPolyLCG` reads the parameter list.

#### `tgfsr.c` / `tgfsr.h`
Twisted Generalized Feedback Shift Register (TGFSR).  A matrix-based GF(2)
recurrence with a twist transformation.  `ReadDataTGFSR` reads degree `k`,
the twist matrix, and the characteristic polynomial.

#### `genf2w.c` / `genf2w.h`
Generator over F_{2^w} — arithmetic is done in the extension field GF(2^w)
rather than GF(2) directly.  `ReadDataGenF2w` reads the field degree `w`,
the recurrence degree `k`, and the characteristic polynomial coefficients.

#### `MT.c` / `MT.h`
Mersenne Twister.  A specific parameterization of TGFSR with a very large
period (2^(km) - 1 for word size m and degree k).  `ReadDataMT` reads the
standard MT parameters (degree `k`, word size `w`, twist matrix B, tempering
parameters if pre-defined).

---

### Tempering Implementations

#### `temperMK.c` / `temperMK.h`
Matsumoto-Kurita (MK) tempering — a systematic output transformation that
improves equidistribution one resolution at a time.

- Supports two variants: type I (`v XOR (v << s) AND b`) and type II (adds
  an upper-triangular step).
- In `opt` mode (`tempMKopt`), the search optimizes `b` to maximize
  equidistribution up to resolution `limitv`.
- **`TestMETemperMK`** — called instead of `TestEquid` when MK optimization
  is requested; iterates the MK optimization loop inside the search.
- **`UpdateTemperMK`** — randomly re-samples the `b` and `c` bit-vectors when
  `aleatoire > 0`.

#### `permut.c` / `permut.h`
Bit permutation tempering — reorders the bits of the output word according to a
fixed permutation.  `InitPermut` generates a random permutation of `precision`
bits.  `UpdatePermut` re-randomizes the permutation at each `UpdateTrans` call.

---

### Mathematics & Support

#### `vectorsF2.c` / `vectorsF2.h`
The GF(2) linear algebra primitives used everywhere:

- **`BitVect`** — a bit-vector of length `n` stored as an array of `ulong`
  words.  Supports XOR, AND, shifts (left/right, rotative), copy, compare,
  and display.
- **`Matrix`** — a 2-D array of `BitVect` rows.  Supports Gaussian elimination
  (`GaussianElimination`, `CompleteElimination`, `SpecialGaussianElimination`),
  matrix-vector multiply, transpose, inverse, and display.
- These are the lowest-level building blocks; every other module depends on them.

#### `polynomials.c` / `polynomials.h`
Polynomial arithmetic over GF(2)[x]:

- `GCDPoly` — GCD of two polynomials.
- `InvModPoly` — modular inverse.
- `PrimitivePoly` — primality test for a polynomial.
- `ProductPoly`, `SumPoly`, `ModPoly` — standard arithmetic.
- Used by generator modules to verify and manipulate characteristic polynomials.

#### `basisreduc.c` / `basisreduc.h`
Lattice basis reduction over GF(2):

- `PermuteCoord` — finds a permutation of the output bits that maximizes the
  equidistribution figures of merit.
- `BasisReduction` — reduces a lattice basis (used in the dual-lattice method
  of `mecf`).
- Only called when `METHOD_DUALLATTICE` is selected.

#### `regpoly.c` / `regpoly.h`
Core numeric utilities:

- `intmin`, `intmax`, `floatmax` — min/max helpers.
- `MersennePrime` — checks whether a degree `k` corresponds to a Mersenne prime.
- `ReadLn` — skips to the next line in a FILE (used pervasively by `readfiles`).
- `MRG32k3a` / `SetMRG32k3a` — the internal MRG32k3a random number generator
  used for random tempering parameter selection during the search.

#### `myzen.c` / `myzen.h`
Wrapper around the external **ZEN** finite-field library.  Provides
`AllocZen`, `FreeZen`, and related calls needed by `genf2w` and `polynomials`
for arithmetic in extension fields GF(2^w).

#### `timer.c` / `timer.h`
Portable CPU timer:

- `timer_Init` — records the start time.
- `timer_Val` — returns elapsed CPU seconds since `timer_Init`.
- Used by `rechercheguide.c` to report search duration.

---

## Data Flow

```
Command line
    │
    ▼
rechercheguide.c: main()
    │
    ├─► ReadSearch()          ← test_file (ME bounds, tuplet params, seeds)
    │       └─► ReadTempering()  ← tempering descriptor files
    │
    ├─► ReadGenDataFiles()    ← gen_data_file(s) (generator parameters)
    │
    ├─► FirstGen()            (initialize Combinaison to first candidate)
    │
    └─► [search loop]
            │
            ├─ UpdateAllTrans()         re-sample tempering parameters
            ├─ TestEquid() / TestMETemperMK()   equidistrib. criterion
            │       └─ PrepareMat() → GaussianElimination()
            │
            ├─ isPresqueME() ─(pass)─► TestTuplets()
            │                               └─ ResolutionEquid() per combination
            │
            ├─ isTuplets() ──(pass)─► DispCurrentComb() + DispTuplets()
            │
            └─ NextGen()              advance to next candidate
```

---

## Supported Generator Types

| Tag in data file | Module | Generator family |
|-----------------|--------|-----------------|
| `taus`    | tausworthe.c | Tausworthe (LFSR) |
| `polylcg` | polylcg.c    | Polynomial LCG |
| `tgfsr`   | tgfsr.c      | Twisted GFSR |
| `genf2w`  | genf2w.c     | Generator over F_{2^w} |
| `MT`      | MT.c         | Mersenne Twister |

## Supported Tempering Types

| Tag in tempering file | Module | Description |
|-----------------------|--------|-------------|
| `permut`      | permut.c   | Random bit permutation |
| `tempMK`      | temperMK.c | Matsumoto-Kurita type I |
| `tempMKopt`   | temperMK.c | MK type I with optimization |
| `tempMK2`     | temperMK.c | Matsumoto-Kurita type II |
| `tempMK2opt`  | temperMK.c | MK type II with optimization |
