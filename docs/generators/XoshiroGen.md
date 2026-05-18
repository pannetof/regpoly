# XoshiroGen

**C++ class:** `XoshiroGen`
**Name in code:** `"Xoshiro Generator"`
**Legacy aliases:** `Xoshiro`

The xoshiro (xor/shift/rotate) linear engine introduced by Blackman &
Vigna (2022). Parametrised by a word width `w`, a state word count
`r ∈ {4, 8}`, and a pair `(A, B)` of shift and rotation amounts.

> Note: the paper writes "k" for the word count; this codebase reserves
> `k` for the state size in bits, so the parameter is named `r`.

## Mathematical recurrence

### r = 4 (S_{4w})

```
t       = s[1] << A
s[2]   ^= s[0]
s[3]   ^= s[1]
s[1]   ^= s[2]
s[0]   ^= s[3]
s[2]   ^= t
s[3]    = rotl(s[3], B)
```

### r = 8 (S_{8w})

```
t       = s[1] << A
s[2]   ^= s[0]
s[5]   ^= s[1]
s[1]   ^= s[2]
s[7]   ^= s[3]
s[3]   ^= s[4]
s[4]   ^= s[5]
s[0]   ^= s[6]
s[6]   ^= s[7]
s[6]   ^= t
s[7]    = rotl(s[7], B)
```

The corresponding `4w × 4w` and `8w × 8w` transition matrices are given
in §3.2 of the paper.

## Parameter sets from the literature

See [`docs/library/blackman-vigna-2022-scrambled.yaml`](../library/blackman-vigna-2022-scrambled.yaml)
for the Table 2 (64-bit) and Table 5 (32-bit) parameter sets.

## Output scramblers

The +, ++ and ** scramblers from §4 of the paper are *not* implemented
here. This class exposes only the F₂-linear engine — the raw
state-array output — which is what equidistribution and primitivity
analyses operate on.

## References

- Blackman & Vigna 2022 — *Scrambled Linear Pseudorandom Number
  Generators*, ACM TOMS 47(4):1-32.
  [`docs/papers/xoroshiro.pdf`](../papers/xoroshiro.pdf)
