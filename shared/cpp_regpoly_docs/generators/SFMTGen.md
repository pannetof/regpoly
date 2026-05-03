# SFMTGen — SIMD-oriented Fast Mersenne Twister

`SFMTGen` is the SIMD-Friendly Mersenne Twister of Saito and Matsumoto
(2008). It produces a 128-bit block per recurrence step using
parallel-prefix shifts and a parameterized mask vector, giving high
output throughput on processors with 128-bit SIMD instructions while
preserving an MT-style Mersenne-prime period of $2^{\text{mexp}}-1$.

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `mexp` | int | Mersenne exponent (e.g. 607, 1279, 2203, 11213, 19937, 44497, 86243, 132049, 216091). |
| `pos1` | int | Linear feedback offset (in 128-bit blocks). |
| `sl1`, `sl2` | int | Left-shift counts for SIMD recurrence (sl2 is byte-shift). |
| `sr1`, `sr2` | int | Right-shift counts (sr2 is byte-shift). |
| `msk1`..`msk4` | uint32 | Per-lane mask values. |
| `parity1`..`parity4` | uint32 | Period certification vector (period adjustment). |

## Recurrence

The SFMT state is a circular array of $N$ 128-bit blocks
(`N = mexp/128 + 1`). Each step computes
$x[i] \gets x[i] \oplus \mathrm{sl_{128}}(x[i]) \oplus \mathrm{sr_{128}}(x[i+\text{pos}_1]) \oplus \mathrm{simd}(x[i+N-2])$
where the SIMD term mixes per-lane shifts and the mask vector.

## Provenance

- **Reference**: Saito, M., Matsumoto, M., "SIMD-oriented Fast Mersenne
  Twister: a 128-bit Pseudorandom Number Generator", Monte Carlo and
  Quasi-Monte Carlo Methods 2006, Springer, 2008.
- **Cross-check fixtures**: `docs/library/sfmt_params.yaml`.
- **Library entry**: `docs/library/saito-matsumoto-2008.yaml`.
