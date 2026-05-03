# DSFMTGen — Double-precision SIMD-oriented Fast Mersenne Twister

`DSFMTGen` is the dSFMT generator of Saito (2009): a SIMD MT variant
that emits 52-bit IEEE-754 mantissa fractions in a 2-lane SIMD packing,
giving direct double-precision floats in [1, 2) without rejection.

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `mexp` | int | Mersenne exponent (521, 1279, 2203, 4253, 11213, 19937, 44497, 86243, 132049, 216091). |
| `pos1` | int | Recurrence offset in 128-bit blocks. |
| `sl1` | int | Left-shift count for the SIMD recurrence. |
| `msk1`, `msk2` | uint64 | Per-lane mask vectors restricting bits to the 52-bit mantissa region. |
| `fix1`, `fix2` | uint64 | Lung-carry fix-up vectors. |
| `pcv1`, `pcv2` | uint64 | Period certification vector. |

## Recurrence

The dSFMT recurrence threads a 128-bit "lung" carry through a circular
buffer of 128-bit blocks. Output extraction packs two 52-bit fractions
per block into IEEE-754 doubles in [1, 2).

## Provenance

- **Reference**: Saito, M., "An application of finite field: design and
  implementation of 128-bit instruction-based fast pseudorandom number
  generator", Master's thesis (Hiroshima University), 2009.
- **Cross-check fixtures**: `docs/library/dsfmt_params.yaml`.
- **Library entry** (planned): `docs/library/saito-2009-dsfmt.yaml`.
