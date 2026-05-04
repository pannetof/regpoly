# TinyMT32 — Tiny Mersenne Twister (32-bit)

**C++ class:** `TinyMT32Gen`
**Name in code:** `"TinyMT32Gen"` (alias: `"TinyMT32"`)

`TinyMT32Gen` is the 127-bit-state Tiny MT variant of Saito and
Matsumoto (2011), later codified as IETF RFC 8682. It is a compact
MT-style F$_2$-linear generator sized for embedded contexts and
parallel-stream applications, yet retaining the prime period
$2^{127} - 1$. Each TinyMT instance is parametrised by a 96-bit
per-stream tuple `(mat1, mat2, tmat)`; the structural shift counts
are fixed constants of the family.

## Mathematical recurrence

The state is four 32-bit words `s[0..3]`, treated as a 127-bit
linear-feedback shift register. At each step:

```
y = s[3]
x = (s[0] & 0x7fffffff) ^ s[1] ^ s[2]
x ^= x << 1
y ^= (y >> 1) ^ x

ns[0] = s[1]
ns[1] = s[2]
ns[2] = x ^ (y << 10)
ns[3] = y
if (y & 1):
    ns[1] ^= mat1
    ns[2] ^= mat2
s = ns

# Tempering:
t1  = s[0] ^ (s[2] >> 8)
out = s[3] ^ t1
if (t1 & 1):
    out ^= tmat
```

The mask `0x7fffffff` on `s[0]` removes the top bit of the first
word so the linear state spans 127 bits (one bit shy of $4 \times
32$). `mat1`, `mat2`, `tmat` are the three per-stream parameters
that distinguish parallel TinyMT instances.

### Characteristic polynomial

For a `(mat1, mat2, tmat)` tuple chosen by the upstream `tinymt32dc`
search tool, $\chi_f$ is the primitive polynomial of degree 127
over $\mathrm{GF}(2)$ guaranteeing the maximum period $2^{127} - 1$.
regpoly recovers the minimal polynomial at runtime via Berlekamp–
Massey on the bit-stream produced by `next()`.

## Parameters

| Name   | Type  | Role   | Rand type | Rand args | Description |
|--------|-------|--------|-----------|-----------|-------------|
| `mat1` | `int` | search | --        | --        | First per-stream matrix entry (XOR'd into `s[1]` when `y` is odd). |
| `mat2` | `int` | search | --        | --        | Second per-stream matrix entry (XOR'd into `s[2]` when `y` is odd). |
| `tmat` | `int` | search | --        | --        | Tempering matrix entry (XOR'd into the output when `t1` is odd). |

Roles:

- **structural** — fixed by the user; defines the shape of $A$.
  TinyMT32 has no user-tunable structural parameters: the state size
  (4 words), word width (32 bits), shifts (1, 8, 10), and top-bit
  mask are all hard-coded constants of the family.
- **search** — randomized by the search loop. The full
  `(mat1, mat2, tmat)` triple is the search axis; the upstream
  `tinymt32dc` parameter-search tool emits triples that yield
  primitive $\chi_f$ of degree 127.

## State size (period exponent)

```
k = 127           (4 * 32 - 1, masking the top bit of s[0])
```

There is one canonical TinyMT32 size; the family does not vary
$k$ across published parametrisations, only the per-stream
$(\mathrm{mat1}, \mathrm{mat2}, \mathrm{tmat})$ tuple.

## Search examples

### Verify the canonical TinyMT32 default tuple

```yaml
search:
  family: TinyMT32Gen
  L: 32
  limit:
    max_tries: 1

# No structural_params — TinyMT32 has none.

fixed_params:
  mat1: 0x8f7011ee
  mat2: 0xfc78ff1f
  tmat: 0x3793fdff
```

## Tempering

The output transformation is the lagged XOR
`t1 = s[0] ^ (s[2] >> 8)` combined with a parameterized
shift-mask via `tmat`: `out = s[3] ^ t1 ^ (tmat if t1 is odd
else 0)`. Unlike MT's two-stage `(>> u) & d`, `(<< s) & b` chain,
the tempering bit width is fixed and `tmat` is the only freely-
chosen tempering coefficient. See
[`theory/tempering_optimization.md`](../theory/tempering_optimization.md)
for the search-side details.

## Implementation notes

- **Top-bit mask.** The `& 0x7fffffff` on `s[0]` is essential to
  the 127-bit linear state interpretation; dropping it would give a
  128-bit reducible state.
- **Seed degeneracy guard.** `init()` sets `s[3] = 1` when the
  supplied seed is the all-zero fixed point of the recurrence,
  mirroring `MTToolBox::tinymt32::seed(1)`. This only fires from
  external `init` calls; the matricial test never feeds an
  all-zero perturbation.

---

## References

- M. Saito and M. Matsumoto. *Tiny Mersenne Twister (TinyMT): a
  small-sized variant of Mersenne Twister.* Hiroshima University
  technical material, released alongside the MTGP paper, 2011.
- M. Saito, M. Matsumoto, V. Roca, and E. Baccelli. *TinyMT32
  Pseudorandom Number Generator (PRNG).* IETF RFC 8682, 2020.
  [DOI](https://doi.org/10.17487/RFC8682).
  Published as
  [`saito-matsumoto-2011-tinymt`](/library/saito-matsumoto-2011-tinymt)
  in the library.
