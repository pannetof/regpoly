# XoroshiroGen

**C++ class:** `XoroshiroGen`
**Name in code:** `"Xoroshiro Generator"`
**Legacy aliases:** `Xoroshiro`

The xoroshiro (xor/rotate/shift/rotate) linear engine introduced by
Blackman & Vigna (2022). Parametrised by a word width `w` and a state
word count `r ≥ 2`, plus a triple `(A, B, C)` of intra-word rotations
and shifts.

> Note: the paper writes "k" for the word count; this codebase reserves
> `k` for the state size in bits (a regpoly-wide convention), so the
> parameter is named `r` instead.

## Mathematical recurrence

State: `r` words `s_0, s_1, ..., s_{r-1}` of `w` bits each. One step of
the recurrence is

```
t = s_0 XOR s_{r-1}
new s_0      = s_1                 (shifted up)
...
new s_{r-3}  = s_{r-2}             (shifted up)
new s_{r-2}  = rotl(s_0, A) XOR t XOR (t << B)
new s_{r-1}  = rotl(t, C)
```

The corresponding `rw × rw` transition matrix `X_{rw}` is given in §3.1
of the paper.

## Parameter sets from the literature

See [`docs/library/blackman-vigna-2022-scrambled.yaml`](../library/blackman-vigna-2022-scrambled.yaml)
for the Table 2 (64-bit) and Table 5 (32-bit) parameter sets.

## Output scramblers

The +, ++, * and ** scramblers from §4 of the paper are *not*
implemented here. This class exposes only the F₂-linear engine — the
raw state-array output — which is what equidistribution and
primitivity analyses operate on.

## References

- Blackman & Vigna 2022 — *Scrambled Linear Pseudorandom Number
  Generators*, ACM TOMS 47(4):1-32.
  [`docs/papers/xoroshiro.pdf`](../papers/xoroshiro.pdf)
