# PolyLCG

**C++ class:** `PolyLCG`
**Name in code:** `"Polynomial LCG"`
**Legacy aliases:** `polylcg`

Polynomial Linear Congruential Generator over GF(2). The state is a
k-bit register that evolves by multiplication by z in the quotient ring
GF(2)[z] / P(z).

## Mathematical recurrence

The state is a polynomial S(z) of degree less than k over GF(2).
At each step:

```
S(z) <- z * S(z)  mod  P(z)
```

where P(z) is the characteristic (feedback) polynomial of degree k.

In terms of bit operations on the k-bit state register:

1. Read the most-significant bit (MSB) of the state.
2. Left-shift the state by 1 position.
3. If the MSB was 1, XOR the state with the feedback polynomial.
4. Mask to k bits.

The characteristic polynomial has the form:

```
P(z) = z^k + z^{e_1} + z^{e_2} + ... + z^{e_t} + 1
```

The generator has maximal period 2^k - 1 when P(z) is primitive over
GF(2).

## Parameters

| Name   | Type      | Role       | Rand type        | Description |
|--------|-----------|------------|------------------|-------------|
| `k`    | `int`     | structural | --               | State size in bits (degree of P(z)) |
| `poly` | `int_vec` | search     | `poly_exponents` | Exponent positions of P(z), excluding the leading term z^k. Specified as a list `[e_1, e_2, ..., 0]` where each entry is an exponent with a nonzero coefficient. The constant term (exponent 0) must always be included. |

### Polynomial encoding

The polynomial P(z) = z^k + z^{e_1} + z^{e_2} + ... + 1 is specified
by listing the exponents of the non-leading terms in descending order:

```
poly: [e_1, e_2, ..., 0]
```

For example, the trinomial z^31 + z^3 + 1 is written as `poly: [3, 0]`.

## State size (period exponent)

```
k = k   (given directly as parameter)
```

The period of the generator is 2^k - 1 when P(z) is primitive.

## Search examples

### Simple search -- find primitive trinomials of degree 89

```yaml
search:
  family: PolyLCG
  L: 32
  limit:
    max_tries: 50000

structural:
  k: 89

search_params:
  poly:          # randomized: random primitive polynomial exponents
```

### Fixed polynomial -- verify a known generator

```yaml
search:
  family: PolyLCG
  L: 32
  limit:
    max_tries: 1

structural:
  k: 31

search_params:
  poly: [3, 0]  # fixed: z^31 + z^3 + 1
```

### Larger degree with output file

```yaml
search:
  family: PolyLCG
  L: 64
  output: polylcg521_results.yaml
  limit:
    max_tries: 100000
    max_seconds: 60

structural:
  k: 521

search_params:
  poly:          # randomized
```
