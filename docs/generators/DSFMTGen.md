# DSFMTGen — Double-precision SIMD-oriented Fast Mersenne Twister

**C++ class:** `DSFMTGen`
**Registered name:** `"DSFMTGen"`
**Legacy aliases:** `dSFMTGen`

`DSFMTGen` is the dSFMT generator of Saito (2009): a SIMD MT variant
that emits 52-bit IEEE-754 mantissa fractions in a 2-lane SIMD packing,
giving direct double-precision floats in $[1, 2)$ without rejection.
The recurrence threads a 128-bit "lung" carry through a circular
buffer of 128-bit blocks and is the standard choice when the consumer
needs a stream of doubles rather than 32-bit integers.

## Mathematical recurrence

The state is an array of $N + 1$ 128-bit blocks, where
$N = \lfloor (\mathrm{mexp} - 128)/104 \rfloor + 1$. The first $N$
blocks hold the active circular buffer; the trailing block is the
"lung" carry $L$. Each 128-bit block exposes two 64-bit lanes
$u[0], u[1]$, of which only the low 52 bits per lane are output bits.

For block index $i$, with $a = \mathrm{state}[i]$,
$b = \mathrm{state}[(i + \mathrm{pos1}) \bmod N]$, and the lung $L$:

```
new_L.u[0] = (a.u[0] << SL1) ^ (L.u[1] >> 32) ^ (L.u[1] << 32) ^ b.u[0]
new_L.u[1] = (a.u[1] << SL1) ^ (L.u[0] >> 32) ^ (L.u[0] << 32) ^ b.u[1]
state[i].u[0] = (new_L.u[0] >> SR) ^ (new_L.u[0] & MSK1) ^ a.u[0]
state[i].u[1] = (new_L.u[1] >> SR) ^ (new_L.u[1] & MSK2) ^ a.u[1]
L = new_L
```

`SR = 12` is hard-coded (matches the upstream dSFMT 2.2.5 reference);
`SL1`, `MSK1`, `MSK2` are per-mexp parameters; the masks restrict the
linear feedback to the 52-bit mantissa region.

### Characteristic polynomial

Computed at runtime via Berlekamp–Massey on the bit-stream produced
by `next()`; the 52-bit mantissa packing and lung carry preclude a
clean closed form. The matricial equidistribution analysis uses the
non-primitive path
([`docs/theory/equidistribution-spec.md`](../theory/equidistribution-spec.md))
on the recovered minimal polynomial.

## Parameters

| Name   | Type  | Role       | Rand type | Rand args | Description |
|--------|-------|------------|-----------|-----------|-------------|
| `mexp` | `int` | structural | --        | --        | Mersenne exponent (521, 1279, 2203, 4253, 11213, 19937, 44497, 86243, 132049, 216091). |
| `pos1` | `int` | structural | --        | --        | Recurrence offset in 128-bit blocks. |
| `sl1`  | `int` | structural | --        | --        | Left-shift count `SL1` for the SIMD recurrence. |
| `msk1` | `int` | structural | --        | --        | 64-bit mask for lane `u[0]` (restricts feedback to mantissa bits). |
| `msk2` | `int` | structural | --        | --        | 64-bit mask for lane `u[1]`. |

Roles:

- **structural** — fixed by the user; defines the shape of $A$.
- **search** — randomized by the search loop. dSFMT exposes no
  search-randomized parameters in the regpoly catalog; per-mexp
  $(\mathrm{pos1}, \mathrm{sl1}, \mathrm{msk1}, \mathrm{msk2})$
  tuples come from the upstream dSFMT 2.2.5 distribution.

## State size (period exponent)

```
k = 128 * (N + 1)        with N = floor((mexp - 128) / 104) + 1
```

The period of the output stream is $2^{\mathrm{mexp}} - 1$; the
F$_2$-linear state space is larger because the lung adds a 128-bit
trailing block.

Classical examples:

- **dSFMT-521:**   mexp = 521   ⇒ $N = 4$,  $k = 128 \cdot 5 = 640$.
- **dSFMT-19937:** mexp = 19937 ⇒ $N = 191$, $k = 128 \cdot 192 = 24576$.
- **dSFMT-44497:** mexp = 44497 ⇒ $N = 426$, $k = 128 \cdot 427 = 54656$.

## Search examples

### Verify a published dSFMT parameter set

```yaml
search:
  family: dSFMTGen
  L: 52
  limit:
    max_tries: 1

structural_params:
  mexp: 521
  pos1: 3
  sl1:  25
  msk1: 0x000fbfefff77efff
  msk2: 0x000ffeebfbdfbfdf
```

## Implementation notes

- **52-bit output width.** Each `next()` call publishes the low 52
  bits of one lane (MSB-first into the BitVect). The Combinaison
  pipeline takes the first $L$ bits; for the natural setting
  `L = 52`, the entire mantissa is consumed.
- **Lung-carry positioning.** The lung is held at trailing index
  $N$ in the state buffer, not interleaved with the active blocks —
  this matches `MTToolBox/samples/dSFMTdc/dSFMTsearch.hpp` and the
  standalone dSFMT 2.2.5 reference at `dSMFT/dSFMT-common.h`.
- **All-zero guard.** `init()` injects a single `1` bit into the
  lung if the seed is degenerate (mirrors dSFMT's
  `period_certification` fallback). `set_raw_state()` skips this
  guard so SIMD-PIS canonical basis vectors stay F$_2$-linear.

---

## References

- M. Saito. *An application of finite field: design and
  implementation of 128-bit instruction-based fast pseudorandom
  number generator.* Master's thesis, Hiroshima University, 2009.
  [DOI](https://doi.org/10.15027/29071).
  Published as [`saito-2009-dsfmt`](../papers/saito-2009-dsfmt.md) in
  the library.
