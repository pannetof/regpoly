# WELL

**C++ class:** `WELLGen`
**Name in code:** `"Carry Generator"` (historical; returned by
`Generator::name()`)
**Aliases accepted by the factory:** `WELLRNG`

Well Equidistributed Long-period Linear (WELL) generators of Panneton,
L'Ecuyer, & Matsumoto (2006). The state is `r` words of `w` bits
stored in a circular buffer. The recurrence uses 8 algorithm slots
T0..T7. Each slot is filled with a transformation class Mi (i=0..6)
from Table I of the paper; the M-class table further down lists each
class's formula and computational cost.

## Mathematical recurrence

Slots T0..T7 are positions in the WELL recurrence; each is filled with
one of the M0..M6 classes per the parameter table. The state consists
of r words V[0], V[1], ..., V[r-1] in a circular buffer indexed by i.
At each step:

```
z0 = (V[(i+r-1) % r] & upper_mask) | (V[(i+r-2) % r] & lower_mask)
z1 = T0(V[i]) XOR T1(V[(i+m1) % r])
z2 = T2(V[(i+m2) % r]) XOR T3(V[(i+m3) % r])
z3 = z1 XOR z2
z4 = T4(z0) XOR T5(z1) XOR T6(z2) XOR T7(z3)

V[(i+r-1) % r] = z4
V[i]           = z3
i              = (i+r-1) % r
```

The masks split the boundary word to achieve the desired state size:

```
upper_mask = bits w-1 down to p    (top w-p bits)
lower_mask = bits p-1 down to 0    (bottom p bits)
```

### Transformation matrix classes

Each of the 8 algorithm slots T0..T7 is filled with one of the seven
M-class transformations from paper Table I. Each Mi takes a small set
of named arguments matching the paper:

| M | Args | Formula |
|---|------|---------|
| 0 | (none) | `y = 0` |
| 1 | (none) | `y = x` (identity) |
| 2 | `t` | `y = x ≫ t` if `t ≥ 0`, else `x ≪ −t` |
| 3 | `t` | `y = x ⊕ shift(x, t)` (same sign rule as M2) |
| 4 | `a` | `y = (x ≫ 1) ⊕ a` if LSB(x), else `x ≫ 1` |
| 5 | `t`, `b` | `y = x ⊕ (shift(x, t) & b)` |
| 6 | `q`, `t`, `s`, `a` | `y = (rotate_left(x, q) & d_s) ⊕ (a if x_t = 1 else 0)` |

For M6, `d_s` is the 32-bit mask with the (s+1)-th bit zero and all
others one (paper Table I's `d_s`); the `x_t = 1` test selects the
bit at MSB-position t. Both are synthesised at runtime from the
small-integer fields, matching the paper exactly.

**Bit-indexing convention.** M4's `x_{w-1}` resolves to the LSB under
the paper's column-vector convention (§1). M6's `t` and `s` are
0-indexed from the MSB.

**Range constraints** (paper Table I):
- M2/M3/M5: `−w ≤ t ≤ w` (signed shift; sign chooses direction)
- M6: `0 ≤ q, t, s < w`
- M4/M5/M6: `a`, `b` ∈ {0, …, 2^w − 1}

### Computational cost

Each Mi class has an associated cost (number of operations):
`M0=0, M1=1, M2=2, M3=3, M4=5, M5=4, M6=8`. Costs are not monotonic
in M-index (M4 > M5) because conditional XOR is more expensive than
masked shift on most architectures. The total cost of a generator is
the sum over its 8 algorithm slots, repeating Mi costs as needed.

## Parameters

| Name       | Type         | Role       | Rand type | Description |
|------------|--------------|------------|-----------|-------------|
| `w`        | `int`        | structural | --        | Word size in bits (default: 32) |
| `r`        | `int`        | structural | --        | Number of words in the circular buffer |
| `p`        | `int`        | structural | --        | Separation point for boundary word masking |
| `m1`       | `int`        | search     | `range` (1 to r-1) | Offset to first tap word |
| `m2`       | `int`        | search     | `range` (1 to r-1) | Offset to second tap word |
| `m3`       | `int`        | search     | `range` (1 to r-1) | Offset to third tap word |
| `matrices` | `struct_map` | search     | `none` — pin every slot, **or** drive the search with a cost cap (see [Cost-bounded search](#cost-bounded-search)) | Map keyed by `T0..T7`; each value is a dict with `M` plus that class's named args |

### `matrices` shape

`matrices` is a YAML / Python map keyed by the slot label (`T0` … `T7`).
Each value is itself a map: an `M` key naming the class (integer 0..6
or string `"M0".."M6"`) plus the per-class arg names from the table
above. Reading the YAML transcribes paper Table II directly:

```yaml
matrices:
  T0: {M: 3, t: -25}      # T0 = M3(-25)
  T1: {M: 3, t:  27}      # T1 = M3(27)
  T2: {M: 2, t:   9}      # T2 = M2(9)
  T3: {M: 3, t:   1}
  T4: {M: 1}              # M1, no args
  T5: {M: 3, t:  -9}
  T6: {M: 3, t: -21}
  T7: {M: 3, t:  21}
```

## State size (period exponent)

```
k = w * r - p
```

Classical examples:
- **WELL19937:** w=32, r=624, p=31 => k = 32*624 - 31 = 19937
- **WELL44497:** w=32, r=1391, p=15 => k = 32*1391 - 15 = 44497
- **WELL512:**  w=32, r=16,  p=0  => k = 32*16 = 512
- **WELL1024:** w=32, r=32,  p=0  => k = 32*32 = 1024

## Search examples

### WELL512a (paper Table II)

```yaml
search:
  family: WELLRNG
  L: 32
  limit:
    max_tries: 1

structural_params:
  w: 32
  r: 16
  p: 0

fixed_params:
  m1: 13
  m2: 9
  m3: 5
  matrices:
    T0: {M: 3, t: -16}
    T1: {M: 3, t: -15}
    T2: {M: 3, t:  11}
    T3: {M: 0}
    T4: {M: 3, t:  -2}
    T5: {M: 3, t: -18}
    T6: {M: 2, t: -28}
    T7: {M: 5, t:  -5, b: 0xda442d24}      # paper a1
```

### WELL19937a (paper Table II)

```yaml
search:
  family: WELLRNG
  L: 32
  output: well19937_results.yaml
  limit:
    max_tries: 1

structural_params:
  w: 32
  r: 624
  p: 31

fixed_params:
  m1: 70
  m2: 179
  m3: 449
  matrices:
    T0: {M: 3, t: -25}
    T1: {M: 3, t:  27}
    T2: {M: 2, t:   9}
    T3: {M: 3, t:   1}
    T4: {M: 1}
    T5: {M: 3, t:  -9}
    T6: {M: 3, t: -21}
    T7: {M: 3, t:  21}
```

### WELL44497a (paper Table II) — uses M6

```yaml
search:
  family: WELLRNG
  L: 32
  limit:
    max_tries: 1

structural_params:
  w: 32
  r: 1391
  p: 15

fixed_params:
  m1: 23
  m2: 481
  m3: 229
  matrices:
    T0: {M: 3, t: -24}
    T1: {M: 3, t:  30}
    T2: {M: 3, t: -10}
    T3: {M: 2, t: -26}
    T4: {M: 1}
    T5: {M: 3, t:  20}
    T6: {M: 6, q: 9, t: 14, s: 5, a: 0xb729fcec}     # paper a7
    T7: {M: 1}
```

## Cost-bounded search

The full-period search has two modes for WELL families:

1. **Pinned matrices** — every slot's `M`-class (and per-class args)
   is fixed in `fixed_params.matrices`; the search only varies the
   integer offsets (`m1`, `m2`, `m3`).
2. **Cost-bounded** — the user sets a `max_cost` budget; the search
   *samples* a fresh `matrices` map under the cap on every iteration,
   alongside the offsets. This is the default workflow on the web UI.

Pinning `matrices` and setting `max_cost` are **mutually exclusive**:
pinning fixes the cost, and `max_cost` only does something when the
search varies matrices. The driver rejects the combination at start.

### YAML example

```yaml
search:
  family: WELLRNG
  L: 32
  limit:
    max_tries: 100000
    max_cost: 12          # 0 < max_cost ≤ 64

structural_params:
  w: 32
  r: 32
  p: 0

fixed_params:
  m1:                     # randomized
  m2:                     # randomized
  m3:                     # randomized
  # `matrices` MUST be omitted under a cost cap.
```

### Sampling algorithm

The C++ helper `WELLGen::random_matrices(w, max_cost, rng)` (also
exposed to Python as `regpoly.well.random_matrices(w, max_cost, seed)`)
draws an admissible `matrices` map:

1. Up to 64 attempts of pure rejection sampling — 8 i.i.d. uniform Mi
   draws from `{M0..M6}`; accept if total cost ≤ `max_cost`.
2. On rejection, fall back to a greedy-budgeted draw: shuffle slot
   indices, then for each slot pick `Mi` uniformly from
   `{Mi : cost(Mi) ≤ remaining_budget}`. This always succeeds because
   `cost(M0) = 0`.

Per-class arg sampling uses the paper's stated ranges (with degenerate
values that collapse the class excluded — e.g. `t = 0` is excluded from
M2/M3/M5).

### Web UI

The `/primitive-search?family=WELLGen` form is **max_cost-only** for
WELL: the `matrices` editor is hidden and a required `max_cost` field
gates submission. To pin matrices, use the YAML/CLI path or POST
directly to `/api/primitive-searches`.

---

## References

- F. Panneton, P. L'Ecuyer, and M. Matsumoto. *Improved long-period
  generators based on linear recurrences modulo 2.* ACM Trans. Math.
  Softw. **32** (2006), 1–16.
  [DOI](https://doi.org/10.1145/1132973.1132974)
