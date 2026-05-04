# MTGP — Mersenne Twister for Graphic Processors

**C++ class:** `MTGPGen`
**Name in code:** `"MTGPGen"` (alias: `"MTGP"`)

`MTGPGen` is the GPU-oriented MT variant of Saito and Matsumoto
(2013). The state is partitioned across small "blocks" so the
recurrence can be evaluated cooperatively by GPU threads, while
preserving an MT-style Mersenne-prime period of
$2^{\mathrm{mexp}} - 1$. The recurrence pairs each step with a
constant-time table lookup that subsumes the standard MT twist.

## Mathematical recurrence

The state is $N = \lfloor \mathrm{mexp}/32 \rfloor + 1$ words of 32
bits, indexed by `idx`. At each step:

```
x1   = state[idx]
x2   = state[(idx + 1) % N]
y_in = state[(idx + pos) % N]
t    = state[(idx + pos - 1) % N]

y    = (x1 & mask) ^ x2 ^ (y_in << sh1)
mat  = tbl[y & 0x0f]
new  = y ^ (y >> sh2) ^ mat
state[idx] = new

# Tempering:
tt   = t ^ (t >> 16)
out  = new ^ tmp_tbl[tt & 0x0f]
idx  = (idx + 1) % N
```

The two 16-entry tables `tbl` and `tmp_tbl` encode both the twist
matrix and the output tempering as 4-bit-indexed XOR masks; a GPU
thread computes them by table lookup rather than by branch on the
low bit.

### Characteristic polynomial

Computed at runtime via Berlekamp–Massey on the bit-stream produced
by `next()`. The combination of the partial-word mask and the
table-driven twist makes a closed form impractical; the matricial
equidistribution analysis uses the non-primitive path
([`docs/theory/equidistribution-spec.md`](../theory/equidistribution-spec.md))
on the recovered minimal polynomial.

## Parameters

| Name      | Type       | Role       | Rand type | Rand args | Description |
|-----------|------------|------------|-----------|-----------|-------------|
| `mexp`    | `int`      | structural | --        | --        | Mersenne exponent (e.g. 11213, 19937, 44497, 86243, 132049, 216091). |
| `pos`     | `int`      | structural | --        | --        | Recurrence offset (in 32-bit words). |
| `sh1`     | `int`      | structural | --        | --        | Pre-twist left shift on `state[(idx+pos) % N]`. |
| `sh2`     | `int`      | structural | --        | --        | Post-twist right shift. |
| `mask`    | `int`      | structural | --        | --        | 32-bit element mask applied to `state[idx]`. |
| `tbl`     | `uint_vec` | structural | `none`    | --        | 16-entry twist-output table (4-bit-indexed XOR masks). |
| `tmp_tbl` | `uint_vec` | structural | `none`    | --        | 16-entry tempering table (4-bit-indexed XOR masks). |

Roles:

- **structural** — fixed by the user; defines the shape of $A$.
- **search** — randomized by the search loop. MTGP exposes no
  search-randomized parameters in the regpoly catalog; published
  per-mexp tuples come from the upstream MTGPDC parameter-search
  tool.

## State size (period exponent)

```
k = 32 * N                       with N = mexp / 32 + 1
```

Classical examples:

- **MTGP-11213:**  mexp = 11213  ⇒ $N = 351$,  $k = 32 \cdot 351 = 11232$.
- **MTGP-23209:**  mexp = 23209  ⇒ $N = 726$,  $k = 32 \cdot 726 = 23232$.
- **MTGP-44497:**  mexp = 44497  ⇒ $N = 1391$, $k = 32 \cdot 1391 = 44512$.

## Search examples

### Verify a published MTGP parameter set

```yaml
search:
  family: MTGPGen
  L: 32
  limit:
    max_tries: 1

structural_params:
  mexp: 11213
  pos:  84
  sh1:  12
  sh2:  4
  mask: 0xfff80000
  tbl:     [0x00000000, 0x71588353, 0xdfa887c1, 0xaef00492,
            0x4ba66c6e, 0x3afeef3d, 0x944eebaf, 0xe51668fc,
            0xa53da0e1, 0xd46523b2, 0x7a952720, 0x0bcda473,
            0xee9bcc8f, 0x9fc34fdc, 0x31734b4e, 0x402bc81d]
  tmp_tbl: [0x00000000, 0x200040bb, 0x1082c61e, 0x308286a5,
            0x10940c5a, 0x30944ce1, 0x0016ca44, 0x20168aff,
            0x00000000, 0x200040bb, 0x1082c61e, 0x308286a5,
            0x10940c5a, 0x30944ce1, 0x0016ca44, 0x20168aff]
```

## Tempering

The output tempering is the second 4-bit table lookup
`tmp_tbl[(t ^ (t >> 16)) & 0x0f]`, where `t` is the lagged state
word `state[(idx + pos - 1) % N]`. The twist table `tbl` plays the
role of the standard MT twist coefficient `a` but indexed by 4 bits
instead of 1, so the linear feedback path itself is part of the
tempering search rather than a single bit-mask. See
[`theory/tempering_optimization.md`](../theory/tempering_optimization.md)
for the search-side details.

## Implementation notes

- **No SIMD-aware overrides.** Unlike SFMT/dSFMT, MTGPGen runs the
  recurrence one word at a time through the same `next()` path used
  by all scalar generators; the "GPU block" partitioning is a
  parameter-search property of the family, not a runtime layout
  choice in regpoly.
- **Cross-check status.** Upstream MTToolBox does not ship a
  per-mexp `calc_equidist` for MTGP — only an aggregate
  parameter-search tool — so the cross-check fixture
  (`docs/library/mtgp_params.yaml`) is currently marked deferred.

---

## References

- M. Saito and M. Matsumoto. *Variants of Mersenne Twister suitable
  for graphic processors.* ACM Trans. Math. Softw. **39** (2013),
  1–20.
  [DOI](https://doi.org/10.1145/2427023.2427030).
  Published as
  [`saito-matsumoto-2013-mtgp`](/library/saito-matsumoto-2013-mtgp)
  in the library.
