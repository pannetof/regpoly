# MarsaXorshiftGen — Marsaglia Xor-shift

**C++ class:** `MarsaXorshiftGen`
**Registered name:** `"MarsaXorshiftGen"` (C++ factory)
**Python-only aliases:** `marsaxorshift` (via `regpoly.core.generator._FAMILY_ALIASES`).
The C++ factory has no `reg_alias` entry for this family.

Marsaglia's Xor-shift generators. A family of generators based on
XOR and shift operations. Multiple types are supported, ranging from
a single-word generator to multi-component generators with arbitrary
shift combinations.

## Mathematical recurrence

The shift function `ShiftR(v, s)` is defined as: if s > 0 then v >> s,
else v << (-s). Positive values give right shifts, negative values give
left shifts.

### Type 1 -- single word

State: a single w-bit word v.

```
v ^= ShiftR(v, a)
v ^= ShiftR(v, b)
v ^= ShiftR(v, c)
```

### Type 2 -- two-component

(Historical note: legacy `.dat` files tag this recurrence with sub-codes
21..25 to drive sibling-variant expansion in the parser; the runtime
collapses all five back to `type: 2`. See
`packages/regpoly-legacy/src/regpoly_legacy/reader.py::_read_marsaxorshift`.)

State: r words in a circular buffer. Two words are selected (x from
position r-1 and y from position m-1) and combined:

```
y = V[m-1]
x = V[r-1]
Right-shift the state by one w-bit word
For j = 0, 1, 2:
    if p[j] != 0: x ^= ShiftR(x, p[j])
For j = 0, 1, 2:
    if q[j] != 0: y ^= ShiftR(y, q[j])
V[0] = x XOR y
```

### Type 3 -- multi-tap

State: r words in a circular buffer. Multiple taps are combined using
64-bit arithmetic to match the legacy C implementation:

```
t = 0
For each tap (position, shift):
    v64 = V[position-1] | (V[position] << 32)   (when position < r)
    v32 = V[position-1]
    t ^= v32 XOR ShiftR(v64, shift)
Right-shift the state by one w-bit word
V[0] = t
```

### Type 4 -- two-component, 2 shifts each

State: r words in a circular buffer. Like type 2 but with exactly
2 shifts per component:

```
y = V[m-1]
x = V[r-1]
Right-shift the state by one w-bit word
x ^= ShiftR(x, p[0])
x ^= ShiftR(x, p[1])
y ^= ShiftR(y, q[0])
y ^= ShiftR(y, q[1])
V[0] = x XOR y
```

### Type 100 -- general

State: r words in a circular buffer. Arbitrary number of components,
each with an arbitrary number of shifts:

```
t = 0
For each entry (position, shifts[]):
    temp = V[position-1]
    For each shift in shifts[]:
        temp ^= ShiftR(temp, shift)
    t ^= temp
Right-shift the state by one w-bit word
V[0] = t
```

## Parameters

| Name             | Type      | Role       | Rand type | Description |
|------------------|-----------|------------|-----------|-------------|
| `type`           | `int`     | structural | --        | Runtime generator type: 1, 2, 3, 4, or 100. (Legacy `.dat` files may write 21..25 to denote sibling-variant rows; the parser normalises those to `type: 2`.) |
| `w`              | `int`     | structural | --        | Word size in bits, ≥ 2. Default: `32`. The runtime uses a fast `uint64_t` kernel for `w ≤ 64`; widths past 64 transparently route through a BitVect-backed wide path (slower but unbounded). Type 3's cross-word read is a w=32-only Marsaglia C artifact; for any other w the tap collapses to a regular one-shift xorshift. |
| `r`              | `int`     | structural | --        | Number of state words. Default: `1`. |
| `m`              | `int`     | structural | --        | Secondary tap position (types 2 and 4). Default: `0`. |
| `shifts`         | `int_vec` | search     | `none`    | Type 1 only: three shift values [a, b, c] |
| `p`              | `int_vec` | search     | `none`    | Types 2, 4: shift values for the first component |
| `q`              | `int_vec` | search     | `none`    | Types 2, 4: shift values for the second component |
| `tap_positions`  | `int_vec` | search     | `none`    | Type 3: positions of taps (1-indexed) |
| `tap_shifts`     | `int_vec` | search     | `none`    | Type 3: shift value for each tap |
| `mi_positions`   | `int_vec` | search     | `none`    | Type 100: positions of components (1-indexed) |
| `mi_shifts`      | `int_vec` | search     | `none`    | Type 100: all shift values, concatenated across components |
| `mi_counts`      | `int_vec` | search     | `none`    | Type 100: number of shifts per component |
| `nb_taps`        | `int`     | enumerator | `none`    | Type 100 search only: number of taps to enumerate. Default: `3`. Ignored at instantiation. |
| `shifts_per_tap` | `int_vec` | enumerator | `none`    | Type 100 search only: per-tap shift count (length must equal `nb_taps`). Default: a vector of 1s. Ignored at instantiation. |

### Type-specific parameter usage

- **Type 1:** uses `shifts` only
- **Type 2:** uses `p` (3 values, padded with 0), `q` (3 values, padded with 0). Zero entries skip the corresponding `(I + S^k)` factor.
- **Type 3:** uses `tap_positions`, `tap_shifts` (parallel arrays, same length). Note: a tap at position `r` falls back to a non-cross-word read at runtime (the `idx + 1 < r` guard in `next()`).
- **Type 4:** uses `p` (2 values), `q` (2 values)
- **Type 100:** uses `mi_positions`, `mi_shifts`, `mi_counts`. Search-time tuning via `nb_taps` and `shifts_per_tap` (see above).

**Note:** Random generation is not supported for any type-specific
parameters (`rand_type=none` for all). All parameters must be provided
explicitly.

## State size (period exponent)

```
k = w * r
```

For type 1, r = 1, so k = w. The period of the generator is 2^k - 1
when the characteristic polynomial is primitive over GF(2).

## Search examples

Since random generation is not supported for this family, the examples
below show verification of known generators.

### Verify a type 1 generator (w=32)

```yaml
search:
  family: MarsaXorshiftGen
  L: 32
  limit:
    max_tries: 1

structural_params:
  type: 1
  w: 32
  r: 1

fixed_params:
  shifts: [13, -17, 5]   # fixed: v ^= (v<<13); v ^= (v>>17); v ^= (v<<5)
```

### Verify a type 4 two-component generator

```yaml
search:
  family: MarsaXorshiftGen
  L: 32
  limit:
    max_tries: 1

structural_params:
  type: 4
  w: 32
  r: 8
  m: 2

fixed_params:
  p: [-2, 17]
  q: [3, -18]
```

### Verify a type 100 general generator

```yaml
search:
  family: MarsaXorshiftGen
  L: 32
  limit:
    max_tries: 1

structural_params:
  type: 100
  w: 32
  r: 5

fixed_params:
  mi_positions: [1, 2, 5]
  mi_shifts: [-28, -17, 3, -18, 13, -12]
  mi_counts: [1, 2, 3]
```
