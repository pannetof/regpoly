# Building MTToolBox `calc_equidist`-style binaries (no autotools)

These binaries are required by `tests/_gen_mttoolbox_reference.py`
to materialize the per-family `*_mttoolbox_reference.py` modules
that the cross-check test asserts against.

The user environment may lack autotools (`aclocal-1.16`, etc.), so
the upstream `make` recipe in each `samples/<family>/Makefile.am`
is bypassed — we invoke `g++` directly.

## Prerequisites (one-time)

Runtime libraries (apt names assume Debian/Ubuntu):
- `libntl` (C++ bindings to NTL)
- `libgmp10`
- `libgf2x3`

If the `-dev` packages are not installed, link directly against the
runtime sonames using `-l:libgmp.so.10 -l:libgf2x.so.3`.

MTToolBox's static lib must exist at
`/home/frpan/projets/claude_projects/regpoly/MTToolBox/lib/libMTToolBox.a`.
Build it once with:

```bash
cd /home/frpan/projets/claude_projects/regpoly/MTToolBox/lib
g++ -O2 -std=c++11 -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS \
    -I../include -c period.cpp -o period.o
g++ -O2 -std=c++11 -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS \
    -I../include -c AlgorithmPrimitivity.cpp -o AlgorithmPrimitivity.o
ar rcs libMTToolBox.a period.o AlgorithmPrimitivity.o
```

(`version.c` requires a missing `config.h` and is not used by any of
the binaries we build.)

## Per-binary recipes

All binaries link against `MTToolBox/lib/libMTToolBox.a` plus NTL,
GMP, GF2X. Common flags:

```bash
CXXFLAGS="-O2 -std=c++11 -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS"
INCLUDE="-I/home/frpan/projets/claude_projects/regpoly/MTToolBox/include"
LIBS="-L/home/frpan/projets/claude_projects/regpoly/MTToolBox/lib \
      -lMTToolBox -lntl -l:libgmp.so.10 -l:libgf2x.so.3"
```

### SFMT — `samples/sfmtdc/calc_equidist`

Needs `Annihilate.o` first.

```bash
cd MTToolBox/samples/sfmtdc
g++ $CXXFLAGS $INCLUDE -c Annihilate.cpp -o Annihilate.o
g++ $CXXFLAGS $INCLUDE -o calc_equidist calc_equidist.cpp Annihilate.o $LIBS
```

CLI: `[-v] [-r] [-s SEED] "mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,parity1,parity2,parity3,parity4"`

Verbose output: `k(NN) = NNNN d(NN) = NNNN` per v=1..32 (32-bit).

### dSFMT — `samples/dSFMTdc/calc_equidist`

```bash
cd MTToolBox/samples/dSFMTdc
g++ $CXXFLAGS $INCLUDE -c Annihilate.cpp -o Annihilate.o
g++ $CXXFLAGS $INCLUDE -o calc_equidist calc_equidist.cpp Annihilate.o $LIBS
```

CLI: `[-v] [-s SEED] "mexp,pos1,sl1,msk1,msk2,fix1,fix2,parity1,parity2"`

Verbose output: `k(NN) = NNNN \t d(NN) = NNNN` per v=1..52 (52-bit).

### RMT64 — `samples/RMT/calc_equidist`

```bash
cd MTToolBox/samples/RMT
g++ $CXXFLAGS $INCLUDE -o calc_equidist calc_equidist.cpp $LIBS
```

CLI: `"mexp,pos,mata,parity,mskb,mskc"` (mexp/pos decimal; rest hex).

Output: `k(N) = NN \t d(N) = NN` per v=1..64.

### MT19937 — `samples/MTDC/mt1`

```bash
cd MTToolBox/samples/MTDC
g++ $CXXFLAGS $INCLUDE -o mt1 mt1.cpp $LIBS
```

No CLI args. Output: `k(NN):NNNN  d(NN):NNNN` (note `:` instead of `=`,
and tightly-packed) per v=1..32.

### XORSHIFT128 — `samples/XORSHIFT/xorshift-{2,5}`

Each binary has its (a,b,c) hardcoded:
- `xorshift-2`: (5, 14, 1)
- `xorshift-5`: (20, 11, 7)

```bash
cd MTToolBox/samples/XORSHIFT
for n in 2 5; do
  g++ $CXXFLAGS $INCLUDE -o xorshift-$n xorshift-$n.cpp $LIBS
done
```

No CLI args. Output: `k(NN):NNN  d(NN):NNN` per v=1..32.

## Verification

After building, smoke-test each binary against a known parameter
set committed in `docs/library/*_params.yaml`. Then run
`python tests/_gen_mttoolbox_reference.py` to emit the per-family
reference modules.
