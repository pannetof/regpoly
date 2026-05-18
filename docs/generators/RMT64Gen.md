# RMT64 — Reducible Mersenne Twister (64-bit)

**C++ class:** `RMT64Gen`
**Name in code:** `"RMT64Gen"` (alias: `"RMT64"`)

`RMT64Gen` is a 64-bit Mersenne-Twister variant with a *reducible*
characteristic polynomial — the period is composite rather than
Mersenne-prime. It is shipped as an MTToolBox sample
(`samples/RMT/rmt64dc`) and exercises the non-primitive
equidistribution analysis path
([`docs/theory/equidistribution-spec.md`](../theory/equidistribution-spec.md))
end-to-end.

## Mathematical recurrence

The state is $N = \lfloor \mathrm{mexp}/64 \rfloor + 1$ words of 64
bits, indexed by `index`. With fixed inner mask
`maska = 0x5555555555555555`:

```
index = (index + 1) % N
x     = state[index]
y     = state[(index + pos) % N]
y    ^= y << 17
mat   = (x & 1) ? mata : 0
x     = y ^ (x >> 1) ^ mat
state[index] = x

# Tempering:
x ^= (x >> 29) & maska
x ^= (x << 17) & maskb
x ^= (x << 37) & maskc
x ^= (x >> 43)
out = x
```

The inner shift-XOR `y ^= y << 17` and the four-stage tempering chain
match `MTToolBox/samples/RMT/rmt64.hpp` line-for-line.

### Characteristic polynomial

Reducible by construction (the "R" in RMT). Recovered at runtime via
Berlekamp–Massey on the bit-stream produced by `next()`; the
matricial equidistribution analysis then factors $\chi_f$ and works
on the chosen invariant subspace via the non-primitive path.

## Parameters

| Name    | Type  | Role       | Rand type | Rand args | Description |
|---------|-------|------------|-----------|-----------|-------------|
| `mexp`  | `int` | structural | --        | --        | Total state degree (period exponent of the dominant irreducible factor). |
| `pos`   | `int` | structural | --        | --        | Recurrence offset in 64-bit words. |
| `mata`  | `int` | structural | --        | --        | 64-bit twist coefficient (matrix-A entry). |
| `maskb` | `int` | structural | --        | --        | 64-bit tempering mask (paired with `<< 17`). |
| `maskc` | `int` | structural | --        | --        | 64-bit tempering mask (paired with `<< 37`). |

Roles:

- **structural** — fixed by the user; defines the shape of $A$.
- **search** — randomized by the search loop. RMT64Gen exposes no
  search-randomized parameters in the regpoly catalog; per-mexp
  tuples are produced by upstream `samples/RMT/rmt64dc` and
  consumed verbatim.

## State size (period exponent)

```
k = 64 * N                       with N = mexp / 64 + 1
```

The `mexp` parameter names the exponent of the dominant irreducible
factor of $\chi_f$ (the "useful" part of the period); the F$_2$-linear
state space is $k = 64 \cdot N$ bits and the full period is the
product of the irreducible factor periods.

Classical examples:

- **RMT-521:**   mexp = 521   ⇒ $N = 9$,   $k = 64 \cdot 9 = 576$.
- **RMT-1279:**  mexp = 1279  ⇒ $N = 20$,  $k = 64 \cdot 20 = 1280$.
- **RMT-19937:** mexp = 19937 ⇒ $N = 312$, $k = 64 \cdot 312 = 19968$.

## Search examples

### Verify a published RMT64 parameter set

```yaml
search:
  family: RMT64Gen
  L: 64
  limit:
    max_tries: 1

structural_params:
  mexp:  1279
  pos:   1
  mata:  0x121b78491deedbc8
  maskb: 0xe3beff7ffbae0000
  maskc: 0xfaf9fae000000000
```

## Tempering

The four-stage tempering chain `(>> 29) & maska`, `(<< 17) & maskb`,
`(<< 37) & maskc`, `>> 43` is part of the recurrence as shipped:
`maskb` and `maskc` are the only freely-chosen masks, the other two
are fixed (`maska = 0x5555…5555` and the trailing right-shift by 43
has no mask). See
[`theory/tempering_optimization.md`](../theory/tempering_optimization.md)
for the search-side details.

## Implementation notes

- **Non-primitive path is mandatory.** RMT64 is the canonical
  exerciser of regpoly's non-primitive matricial equidistribution
  computation; running it with the primitive path will give
  garbage (or fail the spec's primitivity guard).
- **Cross-check delta.** For a subset of mexps, regpoly's
  notprimitive method selects an invariant subspace that differs
  by $\le 1$ bit per $v$ from MTToolBox's annihilate-then-equidist
  convention; the `cross_check_xfail` markers in
  `docs/library/rmt_params.yaml` document the affected entries.
- **Seed degeneracy guard.** `init()` sets `state[0] = 1` when the
  supplied seed is all-zero, mirroring `MTToolBox::seed(1)`; the
  matricial test never feeds an all-zero perturbation, so this only
  triggers from an external `init` call.

---

## References

- S. Harase. *On the F2-linear relations of Mersenne Twister
  pseudorandom number generators.* Math. Comput. Simulation **100**
  (2014), 103–113.
  [DOI](https://doi.org/10.1016/j.matcom.2014.02.002).
  Published as [`mttoolbox-rmt`](../papers/mttoolbox-rmt.md) in the
  library (the practical source of truth for RMT64 parameters
  remains the upstream MTToolBox `samples/RMT/` distribution).
