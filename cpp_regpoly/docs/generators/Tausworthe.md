# Tausworthe

**C++ class:** `Tausworthe`
**Name in code:** `"Tausworthe Generator"`
**Legacy aliases:** `taus`, `taus2`

Decimated LFSR (Linear Feedback Shift Register) generator. The
underlying LFSR is defined by a primitive polynomial over GF(2), and
the output sequence is formed by decimating the LFSR bit stream with
step s to produce L-bit output words.

## Mathematical recurrence

The LFSR is defined by a characteristic polynomial over GF(2):

```
P(z) = z^k + z^{q_1} + z^{q_2} + ... + z^{q_{t-1}} + 1
```

where k > q_1 > q_2 > ... > q_{t-1} > 0. The LFSR bit sequence
{x_n} satisfies:

```
x_n = x_{n-k+q_1} XOR x_{n-k+q_2} XOR ... XOR x_{n-k+q_{t-1}}    (for n >= k)
```

The output is formed by decimation with step s. The n-th output word
u_n of L bits is:

```
u_n = (x_{ns}, x_{ns+1}, ..., x_{ns+L-1})
```

When the polynomial is primitive, the LFSR has maximal period 2^k - 1.

### Quick Tausworthe mode

When `quicktaus=true` (the default), the decimation step is
automatically derived from the polynomial to enable an efficient
implementation:

```
s = k - Q[-2]
```

where Q is the sorted list of polynomial exponents and Q[-2] is the
second-largest exponent. This choice allows the recurrence to be
computed using only whole-word XOR and shift operations.

## Parameters

| Name        | Type      | Role       | Rand type | Description |
|-------------|-----------|------------|-----------|-------------|
| `poly`      | `int_vec` | structural | --        | Polynomial exponents in any order. Must include the degree k as the largest element, plus all intermediate exponents and 0. Example: `[k, q_1, q_2, ..., 0]`. |
| `s`         | `int`     | has default | --       | Decimation step. When omitted and `quicktaus=true`, automatically set to `k - Q[-2]`. When omitted and `quicktaus=false`, defaults to 1. |
| `quicktaus` | `bool`    | structural | --        | Enable quick Tausworthe mode (default: `true`). When true, uses an optimized recurrence that requires s = k - Q[-2]. |

### How k is derived

The state size k is not given directly. It is the largest element of
the `poly` vector:

```
k = max(poly)
```

## Search examples

### Simple trinomial-based Tausworthe generator

```yaml
search:
  family: Tausworthe
  L: 32
  limit:
    max_tries: 1

structural_params:
  poly: [31, 3, 0]
  quicktaus: true

fixed_params:
  # No search parameters -- all are structural or auto-derived.
  # s is automatically set to 31 - 3 = 28.
```

### Pentanomial with explicit decimation step

```yaml
search:
  family: Tausworthe
  L: 32
  limit:
    max_tries: 1

structural_params:
  poly: [89, 38, 19, 3, 0]
  quicktaus: false

fixed_params:
  s: 13   # fixed decimation step
```

### Quick Tausworthe with a known Mersenne exponent

```yaml
search:
  family: Tausworthe
  L: 64
  limit:
    max_tries: 1

structural_params:
  poly: [521, 32, 0]
  quicktaus: true

fixed_params:
  # s auto-derived: 521 - 32 = 489
```
