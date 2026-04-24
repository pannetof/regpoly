# WELLRNG

**C++ class:** `WELLRNG`
**Name in code:** `"Carry Generator"`
**Legacy aliases:** `carry`, `Carry2Gen`

Well Equidistributed Long-period Linear (WELL) generators. The state
is r words of w bits stored in a circular buffer. The recurrence uses
8 transformation matrices (T0 through T7) applied to combinations of
state words at configurable offsets, providing high flexibility in
tuning equidistribution properties.

## Mathematical recurrence

The state consists of r words V[0], V[1], ..., V[r-1] in a circular
buffer indexed by i. At each step:

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

Each of the 8 matrices T0-T7 is one of the following types:

| Type | Name       | Formula | Parameters |
|------|------------|---------|------------|
| 0    | T0         | `v XOR (v >> s)` | `s` (shift amount, in `paramsint[0]`) |
| 1    | Identity   | `v` | none |
| 2    | T2         | `if LSB(v) then (v >> 1) XOR a else (v >> 1)` | `a` (mask, in `paramsulong[0]`) |
| 3    | T3         | `v >> s` | `s` (shift amount, in `paramsint[0]`) |
| 4    | T4         | `v XOR ((v >> s) & a)` | `s` (in `paramsint[0]`), `a` (mask, in `paramsulong[0]`) |
| 5    | T5         | `rotate(v, s) & b XOR (if v & c then a else 0)` | `s` (in `paramsint[0]`), `a` (in `paramsulong[0]`), `b` (in `paramsulong[1]`), `c` (in `paramsulong[2]`) |
| 6    | T6         | `v XOR (v >> s1) XOR (v >> s2) XOR (v >> s3)` | `s1, s2, s3` (in `paramsint[0..2]`) |
| 7    | ZERO       | `0` | none |

Each matrix entry has:
- `type`: integer 0-7 selecting the formula
- `paramsint[0..2]`: up to 3 integer parameters (shift amounts)
- `paramsulong[0..2]`: up to 3 unsigned parameters (bitmasks)

### Computational cost

Each type has an associated cost (number of operations):
T0=3, T1=1, T2=5, T3=2, T4=4, T5=8, T6=7, T7=0. The total cost of a
generator is the sum of the costs of its 8 matrices.

## Parameters

| Name        | Type       | Role       | Rand type | Description |
|-------------|------------|------------|-----------|-------------|
| `w`         | `int`      | structural | --        | Word size in bits (default: 32) |
| `r`         | `int`      | structural | --        | Number of words in the circular buffer |
| `p`         | `int`      | structural | --        | Separation point for boundary word masking |
| `m1`        | `int`      | search     | `range` (1 to r-1) | Offset to first tap word |
| `m2`        | `int`      | search     | `range` (1 to r-1) | Offset to second tap word |
| `m3`        | `int`      | search     | `range` (1 to r-1) | Offset to third tap word |
| `mat_types` | `int_vec`  | search     | `none` (must be provided) | List of 8 integers: the type (0-7) for each matrix T0-T7 |
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

### WELL512-class with known matrix configuration

```yaml
search:
  family: WELLRNG
  L: 32
  limit:
    max_tries: 50000

structural_params:
  w: 32
  r: 16
  p: 0

fixed_params:
  m1:                # randomized in [1, 15]
  m2:                # randomized in [1, 15]
  m3:                # randomized in [1, 15]
  mat_types: [0, 1, 0, 3, 4, 1, 6, 7]
  mat_pi:
    # T0(s=13), T1=Id, T2=T0(s=9), T3=T3(s=5), T4=T4(s=18), T5=Id, T6=T6(s1=1,s2=7,s3=15), T7=Zero
    [13, 0, 0,   0, 0, 0,   9, 0, 0,   5, 0, 0,   18, 0, 0,   0, 0, 0,   1, 7, 15,   0, 0, 0]
  mat_pu:
    # Unsigned params: only T4 needs a mask at index 4*3=12
    [0, 0, 0,   0, 0, 0,   0, 0, 0,   0, 0, 0,   0xD3E43FCD, 0, 0,   0, 0, 0,   0, 0, 0,   0, 0, 0]
```

### WELL19937 with fixed offsets and matrix types

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
  mat_types: [0, 1, 0, 3, 4, 1, 6, 7]
  mat_pi:  [25, 0, 0,   0, 0, 0,   27, 0, 0,   7, 0, 0,   22, 0, 0,   0, 0, 0,   3, 11, 19,   0, 0, 0]
  mat_pu:  [0, 0, 0,    0, 0, 0,   0, 0, 0,    0, 0, 0,    0xE73C5628, 0, 0,   0, 0, 0,   0, 0, 0,   0, 0, 0]
```

### Minimal cost WELL1024 (only identity and zero matrices)

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
  mat_types: [0, 1, 3, 3, 0, 1, 0, 7]
  mat_pi:  [16, 0, 0,   0, 0, 0,   15, 0, 0,   11, 0, 0,   19, 0, 0,   0, 0, 0,   4, 0, 0,   0, 0, 0]
  mat_pu:  [0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,   0, 0, 0,   0, 0, 0]
```

---

## References

- F. Panneton, P. L'Ecuyer, and M. Matsumoto. *Improved long-period
  generators based on linear recurrences modulo 2.* ACM Trans. Math.
  Softw. **32** (2006), 1–16.
  [DOI](https://doi.org/10.1145/1132973.1132974)
