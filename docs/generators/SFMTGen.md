# SFMT — SIMD-oriented Fast Mersenne Twister

**C++ class:** `SFMTGen`
**Name in code:** `"SFMTGen"` (alias: `"SFMT"`)

`SFMTGen` is the SIMD-Friendly Mersenne Twister of Saito and
Matsumoto (2008). It produces a 128-bit block per recurrence step
using parallel-prefix shifts and a parameterized lane mask, giving
high output throughput on processors with 128-bit SIMD instructions
while preserving an MT-style Mersenne-prime period of
$2^{\mathrm{mexp}} - 1$. Output is consumed 32 bits at a time
(four 32-bit lanes per block) by downstream code.

## Mathematical recurrence

The state is a circular array of $N = \lfloor \mathrm{mexp}/128
\rfloor + 1$ blocks of 128 bits, each block exposing four 32-bit
lanes `u[0..3]` (`u[0]` is the low 32 bits). For block index $i$,
with $a = \mathrm{state}[i]$,
$b = \mathrm{state}[(i + \mathrm{pos1}) \bmod N]$,
$c = \mathrm{state}[(i + N - 2) \bmod N]$,
$d = \mathrm{state}[(i + N - 1) \bmod N]$:

```
x = lshift128(a, sl2)        # 128-bit-wide left  byte-shift by sl2 bytes
y = rshift128(c, sr2)        # 128-bit-wide right byte-shift by sr2 bytes

for l in 0..3:
    state[i].u[l] = a.u[l]
                  ^ x.u[l]
                  ^ ((b.u[l] >> sr1) & msk[l])
                  ^ y.u[l]
                  ^ (d.u[l] << sl1)
```

`sl2`/`sr2` are byte-granular shifts (the SIMD `pslldq`/`psrldq`
operations); `sl1`/`sr1` are bit-granular per-lane shifts; the
4-lane mask vector `(msk1, msk2, msk3, msk4)` parametrises the
linear feedback. After `generate_all()` runs the recurrence over
all $N$ blocks, `next()` consumes successive 32-bit lanes from the
buffer, refilling on exhaustion.

### Characteristic polynomial

Computed at runtime via Berlekamp–Massey on the bit-stream produced
by `next()`; the SIMD layout precludes a clean closed form. The
matricial equidistribution analysis uses the non-primitive path
([`docs/theory/equidistribution-spec.md`](../theory/equidistribution-spec.md))
on the recovered minimal polynomial.

## Parameters

| Name   | Type       | Role       | Rand type | Rand args | Description |
|--------|------------|------------|-----------|-----------|-------------|
| `mexp` | `int`      | structural | --        | --        | Mersenne exponent (607, 1279, 2281, 4253, 11213, 19937, 44497, 86243, 132049, 216091). |
| `pos1` | `int`      | structural | --        | --        | Linear-feedback offset in 128-bit blocks. |
| `sl1`  | `int`      | structural | --        | --        | Per-lane left-shift bit count. |
| `sl2`  | `int`      | structural | --        | --        | 128-bit-wide left byte-shift count. |
| `sr1`  | `int`      | structural | --        | --        | Per-lane right-shift bit count. |
| `sr2`  | `int`      | structural | --        | --        | 128-bit-wide right byte-shift count. |
| `msk`  | `uint_vec` | structural | `none`    | --        | 4-element lane mask `(msk1, msk2, msk3, msk4)`. |

Roles:

- **structural** — fixed by the user; defines the shape of $A$.
- **search** — randomized by the search loop. SFMT exposes no
  search-randomized parameters in the regpoly catalog; per-mexp
  tuples come from the upstream sfmt-1.5.1 `SFMT-paramsNNN.h`
  headers. Per-mexp parity vectors used by `period_certify()`
  are looked up internally and need not be supplied.

## State size (period exponent)

```
k = 128 * N                      with N = mexp / 128 + 1
```

Classical examples:

- **SFMT-607:**   mexp = 607   ⇒ $N = 5$,   $k = 128 \cdot 5 = 640$.
- **SFMT-19937:** mexp = 19937 ⇒ $N = 156$, $k = 128 \cdot 156 = 19968$.
- **SFMT-44497:** mexp = 44497 ⇒ $N = 348$, $k = 128 \cdot 348 = 44544$.

## Search examples

### Verify a published SFMT parameter set

```yaml
search:
  family: SFMTGen
  L: 32
  limit:
    max_tries: 1

structural_params:
  mexp: 607
  pos1: 2
  sl1:  15
  sl2:  3
  sr1:  13
  sr2:  3
  msk:  [0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f]
```

### Verify a high-period SFMT19937

```yaml
search:
  family: SFMTGen
  L: 32
  output: sfmt19937_results.yaml
  limit:
    max_tries: 1

structural_params:
  mexp: 19937
  pos1: 122
  sl1:  18
  sl2:  1
  sr1:  11
  sr2:  1
  msk:  [0xdfffffef, 0xddfecb7f, 0xbffaffff, 0xbffffff6]
```

## Tempering

SFMT folds its tempering into the recurrence itself — the per-lane
mask `msk` and the shifts `sr1`/`sl1` together play the role of the
output transformation. There is no separate post-state tempering
pass (contrast with MT's two shift-mask steps and dSFMT's lung-XOR
step). See
[`theory/tempering_optimization.md`](../theory/tempering_optimization.md)
for the search-side details that apply when treating `msk` as a
tunable.

## Implementation notes

- **Lane / byte-shift conventions.** `lshift128` and `rshift128`
  treat the four 32-bit lanes as a single 128-bit value with `u[0]`
  the low 32 bits, mirroring `_mm_slli_si128` / `_mm_srli_si128`
  semantics. `sl2` and `sr2` count bytes (multiplied by 8 internally
  to get bit shifts).
- **Buffer refill.** `next()` lazy-runs `generate_all()` only when
  the 4-lane buffer is exhausted (`buf_idx_ >= 4 * N`); the
  matricial test forces the refill explicitly via
  `simd_advance_one_word()`.
- **Parity certification.** `period_certify()` adjusts a single
  bit of `state[0]` to land in the main orbit when the parity
  inner-product is zero, using a per-mexp parity table baked into
  the source (no user-supplied parity required).

---

## References

- M. Saito and M. Matsumoto. *SIMD-oriented Fast Mersenne Twister:
  a 128-bit Pseudorandom Number Generator.* Monte Carlo and
  Quasi-Monte Carlo Methods 2006, Springer (2008), 607–622.
  [DOI](https://doi.org/10.1007/978-3-540-74496-2_36).
  Published as
  [`saito-matsumoto-2008`](/library/saito-matsumoto-2008) in the
  library.
