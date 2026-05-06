# WELLRNG

**C++ class:** `WELLRNG`
**Name in code:** `"Carry Generator"`
**Legacy aliases:** `carry`, `Carry2Gen`

Well Equidistributed Long-period Linear (WELL) generators. The state
is r words of w bits stored in a circular buffer. The recurrence uses
8 algorithm slots T0..T7. Each slot is filled with a transformation
class Mi (i=0..6) from Table I of Panneton et al. (2006); the M-class
table further down lists each class's formula.

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

### Transformation matrix types

Each of the 8 algorithm slots T0..T7 is filled with one of the seven
M-class transformations from paper Table I:

| Type | Paper Mi         | Formula | Parameters |
|------|------------------|---------|------------|
| 0    | M0               | `y = 0` | none |
| 1    | M1               | `y = x` (identity) | none |
| 2    | M2(t)            | `y = x ≫ t` if `t ≥ 0`, else `x ≪ −t` | `t` ∈ `paramsint[0]` |
| 3    | M3(t)            | `y = x ⊕ shift(x, t)` (same sign rule as M2) | `t` ∈ `paramsint[0]` |
| 4    | M4(a)            | `y = (x ≫ 1) ⊕ a` if LSB(x), else `x ≫ 1` | `a` ∈ `paramsulong[0]` |
| 5    | M5(t, b)         | `y = x ⊕ (shift(x, t) & b)` | `t` ∈ `paramsint[0]`, `b` ∈ `paramsulong[0]` |
| 6    | M6\*(q, a, b, c) | `y = (rotate_left(x, q) & b) ⊕ (a if (x & c) else 0)` | `q` ∈ `paramsint[0]`, `a, b, c` ∈ `paramsulong[0..2]` |

**M6\*** is a strict generalisation of paper Table I's M6: arguments
`b` and `c` are full 32-bit masks, supersetting the paper's
parametric `d_s` and `x_t = 1` test. Setting `b = ~0`,
`c = (1 << t)`, and `a` as published reproduces the exact paper
formula.

**Bit-indexing convention.** `LSB(x)` corresponds to the paper's
`x_{w-1}` under the column-vector convention of §1.

Each matrix entry has:
- `type`: integer 0-6 selecting the Mi class
- `paramsint[0..2]`: up to 3 integer parameters (shift amounts)
- `paramsulong[0..2]`: up to 3 unsigned parameters (bitmasks)

### Computational cost

Each Mi class has an associated cost (number of operations):
M0=0, M1=1, M2=2, M3=3, M4=5, M5=4, M6=8. Costs are not monotonic in
M-index (M4 > M5) because conditional XOR is more expensive than
masked shift on most architectures. The total cost of a generator is
the sum over its 8 algorithm slots, repeating Mi costs as needed.

## Parameters

| Name        | Type       | Role       | Rand type | Description |
|-------------|------------|------------|-----------|-------------|
| `w`         | `int`      | structural | --        | Word size in bits (default: 32) |
| `r`         | `int`      | structural | --        | Number of words in the circular buffer |
| `p`         | `int`      | structural | --        | Separation point for boundary word masking |
| `m1`        | `int`      | search     | `range` (1 to r-1) | Offset to first tap word |
| `m2`        | `int`      | search     | `range` (1 to r-1) | Offset to second tap word |
| `m3`        | `int`      | search     | `range` (1 to r-1) | Offset to third tap word |
| `mat_types` | `int_vec`  | search     | `none` (must be provided) | List of 8 integers: the Mi class index (0–6) for each algorithm slot T0..T7 |
| `mat_pi`    | `int_vec`  | search     | `none` (must be provided) | Flat list of 24 integers: 3 integer params per matrix (8 x 3) |
| `mat_pu`    | `uint_vec` | search     | `none` (must be provided) | Flat list of 24 unsigned integers: 3 unsigned params per matrix (8 x 3) |

### Matrix parameter encoding

The `mat_pi` and `mat_pu` vectors are flat arrays of length 24 (8
matrices x 3 parameters each). For matrix j (0-indexed), the
parameters are at indices `j*3`, `j*3+1`, `j*3+2`. Unused parameters
should be set to 0.

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
  # T0=M3(-16), T1=M3(-15), T2=M3(11), T3=M0,
  # T4=M3(-2),  T5=M3(-18), T6=M2(-28), T7=M5(-5, a1)
  mat_types: [3, 3, 3, 0, 3, 3, 2, 5]
  mat_pi:    [-16, 0, 0,  -15, 0, 0,  11, 0, 0,    0, 0, 0,
               -2, 0, 0,  -18, 0, 0,  -28, 0, 0,  -5, 0, 0]
  # Only slot T7 (an M5 instance) carries a mask at index 7*3=21.
  mat_pu:    [0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,
              0, 0, 0,    0, 0, 0,    0, 0, 0,    0xda442d24, 0, 0]
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
  # T0=M3(-25), T1=M3(27), T2=M2(9), T3=M3(1),
  # T4=M1,      T5=M3(-9), T6=M3(-21), T7=M3(21)
  mat_types: [3, 3, 2, 3, 1, 3, 3, 3]
  mat_pi:    [-25, 0, 0,  27, 0, 0,   9, 0, 0,   1, 0, 0,
                0, 0, 0,  -9, 0, 0, -21, 0, 0,  21, 0, 0]
  mat_pu:    [0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,
              0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0]
```

### Minimal-cost WELL1024 search shape (synthetic)

This is not a published generator; it shows the parameter shape for a
search that explores only the cheap M-classes (M0, M1, M2, M3).

```yaml
search:
  family: WELLRNG
  L: 32
  limit:
    max_tries: 100000

structural_params:
  w: 32
  r: 32
  p: 0

fixed_params:
  m1:                # randomized
  m2:                # randomized
  m3:                # randomized
  # Slot classes: M3, M1, M2, M2, M3, M1, M3, M0 — no masks needed.
  mat_types: [3, 1, 2, 2, 3, 1, 3, 0]
  mat_pi:    [16, 0, 0,   0, 0, 0,  15, 0, 0,  11, 0, 0,
              19, 0, 0,   0, 0, 0,   4, 0, 0,   0, 0, 0]
  mat_pu:    [0, 0, 0,    0, 0, 0,   0, 0, 0,   0, 0, 0,
              0, 0, 0,    0, 0, 0,   0, 0, 0,   0, 0, 0]
```

---

## References

- F. Panneton, P. L'Ecuyer, and M. Matsumoto. *Improved long-period
  generators based on linear recurrences modulo 2.* ACM Trans. Math.
  Softw. **32** (2006), 1–16.
  [DOI](https://doi.org/10.1145/1132973.1132974)
