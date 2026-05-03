# MatsumotoGen

**C++ class:** `MatsumotoGen`
**Name in code:** `"Matsumoto Generator"`
**Legacy aliases:** `matsumoto`

Matsumoto-style generators using 32-bit words with word-level shifts
and XOR operations. Three types are supported, each with a different
recurrence structure. The state is a circular buffer of 32-bit words.

## Mathematical recurrence

All types use 32-bit words. The shift function `SHIFT(v, i)` is
defined as: if i > 0 then v >> i, else v << (-i). This means positive
values give right shifts and negative values give left shifts.

### Type 1

State: V[0], V[1], ..., V[n]. Recurrence:

```
z0 = V[n] XOR SHIFT(V[0], p[0]) XOR SHIFT(V[0], p[1])
z1 = z0 XOR V[m]
z2 = z1 XOR SHIFT(z1, p[2])
z3 = V[0] XOR z2 XOR SHIFT(z2, p[3])
Left-shift the state by one 32-bit word
V[n-1] = z3
V[n]   = z2
```

### Type 2

State: V[0], V[1], ..., V[n]. Like type 1, but with a conditional
twist coefficient `a` (stored in `paramsunsigned[0]`) applied when
`V[n]` is odd:

```
z0 = V[m]
if V[n] is odd:
    z1 = (V[n] >> 1) XOR SHIFT(z0, p[0]) XOR SHIFT(z0, p[1]) XOR a
else:
    z1 = (V[n] >> 1) XOR SHIFT(z0, p[0]) XOR SHIFT(z0, p[1])
z2 = V[0]
z3 = z2 XOR SHIFT(z2, p[2])
z4 = z3 XOR SHIFT(z3, p[3])
Left-shift the state by one 32-bit word
V[n-1] = z1
V[n]   = z4
```

### Type 3

State: V[0], V[1], ..., V[n], V[n+1], V[n+2]. Uses 3 extra words for
a 4-step recurrence:

```
v0 = V[0]
z0 = V[n] XOR v0 XOR SHIFT(V[m], p[0])
z1 = V[n+1] XOR z0 XOR SHIFT(z0, p[1])
z2 = V[n+2] XOR z1 XOR SHIFT(z1, p[2])
Left-shift the state by one 32-bit word
V[n-1] = v0 XOR z0 XOR z1 XOR z2
V[n]   = z0
V[n+1] = z1
V[n+2] = z2
```

## Parameters

| Name              | Type       | Role       | Rand type | Description |
|-------------------|------------|------------|-----------|-------------|
| `type`            | `int`      | structural | --        | Generator type: 1, 2, or 3 |
| `n`               | `int`      | structural | --        | Number of main state words minus 1 (state indices run from 0 to n for types 1,2; from 0 to n+2 for type 3) |
| `m`               | `int`      | search     | `range`   | Feedback tap position, in the range [1, n-1] |
| `paramsint`       | `int_vec`  | search     | `none`    | Shift parameters p[0], p[1], p[2], ... (positive = right shift, negative = left shift). The number of parameters depends on the type. |
| `paramsunsigned`  | `uint_vec` | search     | `none`    | Unsigned parameters (e.g., twist coefficient `a` for type 2). May be empty for types 1 and 3. |

### Shift parameter conventions

- **Type 1:** 4 shift values: p[0], p[1], p[2], p[3]
- **Type 2:** 4 shift values: p[0], p[1], p[2], p[3]; plus `paramsunsigned[0]` = twist coefficient `a`
- **Type 3:** 3 shift values: p[0], p[1], p[2]

**Note:** Random generation for `paramsint` and `paramsunsigned` is not
supported (`rand_type=none`). These parameters must be provided
explicitly for each search attempt, or the search must fix them.

## State size (period exponent)

```
k = (n + 1) * 32     for types 1 and 2
k = (n + 3) * 32     for type 3
```

The period of the generator is 2^k - 1 when the characteristic
polynomial is primitive over GF(2).

## Search examples

### Search over m with fixed shifts -- type 1

```yaml
search:
  family: MatsumotoGen
  L: 32
  limit:
    max_tries: 1000

structural_params:
  type: 1
  n: 24

fixed_params:
  m:                                # randomized: uniform in [1, n-1]
  paramsint: [-15, -17, 13, -7]    # fixed shift parameters
```

### Verify a known type 2 generator (TT800-like)

```yaml
search:
  family: MatsumotoGen
  L: 32
  limit:
    max_tries: 1

structural_params:
  type: 2
  n: 24

fixed_params:
  m: 7                                # fixed
  paramsint: [1, -15, 11, -17]        # fixed shift parameters
  paramsunsigned: [0x8ebfd028]        # fixed twist coefficient a
```

### Type 3 with output file

```yaml
search:
  family: MatsumotoGen
  L: 32
  output: matsumoto_type3_results.yaml
  limit:
    max_tries: 5000
    max_seconds: 60

structural_params:
  type: 3
  n: 15

fixed_params:
  m:                                # randomized: uniform in [1, n-1]
  paramsint: [-17, 15, -22]        # fixed shift parameters
```
