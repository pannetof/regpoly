# Equidistribution Report — Bhuvaneswari & Bhattacharjee 2026

All equidistribution values below come from the **C++ production kernel** (`EquidistributionTest.run()` with the matricial DE method, which is the default for `CellularAutomataGen`).  The per-resolution gaps `Δ_ℓ` are converted to per-dimension gaps `Λ_t` via the production `EquidistributionResults._conv_ecarts(C)` method (no reimplementation).

- Per-component `L = min(k, 64)`.
- Combined `L = min(L_1, L_2)` (= `comb.L` after `Combination._update_stats`).
- A generator is **maximally equidistributed (ME)** iff `Λ_t = 0` for every `t ∈ Φ_1 ∪ Φ_2` per the paper's Section 2.3.2 definition.
- `paper rank` (where shown) is the rank value verbatim from the paper's Table 1.  This is paper data, not our computation; we have not been able to reproduce these specific rank values under the matrix construction the paper itself describes (`t·ℓ × k` matrix from `t` advances).  The discrepancy is unexplained — possibly paper-side typos, possibly an undocumented construction detail.

**Total generators tested: 115**

## Summary

| Paper section | Generators | Paper ME | Kernel ME | Match |
|---|---:|---:|---:|---:|
| Section 3.1 / Table 1 | 1 | 0 | 0 | 1/1 |
| Section 3.1 (weak) | 2 | 0 | 0 | 2/2 |
| Table 3 (s=1) | 1 | 0 | 0 | 1/1 |
| Table 7 (s=2) | 1 | 0 | 0 | 1/1 |
| Table 8 (s=4) | 1 | 0 | 0 | 1/1 |
| Table 9 (s=7) | 1 | 1 | 0 | 0/1 |
| Table 10 (s=8) | 1 | 1 | 0 | 0/1 |
| Table 12 | 49 | 49 | 0 | 0/49 |
| Table 11 | 58 | 46 | 0 | 12/58 |

## Section 3.1 — Weak generators (R1, R3, R4)

### R1 (k=32 single CA)

- **id**: `bhuv2026-r1`
- **family**: single
- **k** = 32; rule-150 positions: `[1, 5, 6, …(16 positions)…, 27, 29, 31]`
- **s** = 1; **L** (combined) = 32
- **Paper claim**: NOT ME
- **Kernel result**: NOT ME (ecart_sum = 56)

| t | ℓ*_t | Λ_t (kernel) | ℓ_t (kernel) | paper rank | t·ℓ*_t | paper equi |
|---:|---:|---:|---:|---:|---:|:---|
| 2 | 16 | 15 | 1 | 18 | 32 | NOT |
| 3 | 10 | 9 | 1 | 13 | 30 | NOT |
| 4 | 8 | 7 | 1 | 12 | 32 | NOT |
| 5 | 6 | 5 | 1 | 11 | 30 | NOT |
| 6 | 5 | 4 | 1 | 11 | 30 | NOT |
| 8 | 4 | 3 | 1 | 12 | 32 | NOT |
| 10 | 3 | 2 | 1 | 13 | 30 | NOT |
| 16 | 2 | 1 | 1 | 17 | 32 | NOT |
| 32 | 1 | 0 | 1 | 32 | 32 | equi |

### R3 (k=35, CA(150′))

- **id**: `bhuv2026-r3`
- **family**: single
- **k** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **s** = 1; **L** (combined) = 35
- **Paper claim**: NOT ME
- **Kernel result**: NOT ME (ecart_sum = 62)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 17 | 16 | 1 |
| 3 | 11 | 10 | 1 |
| 4 | 8 | 7 | 1 |
| 5 | 7 | 6 | 1 |
| 7 | 5 | 4 | 1 |
| 8 | 4 | 3 | 1 |
| 11 | 3 | 2 | 1 |
| 17 | 2 | 1 | 1 |

### R4 (k=64, rule-150 at cells 3,5)

- **id**: `bhuv2026-r4`
- **family**: single
- **k** = 64; rule-150 positions: `[2, 4]`
- **s** = 1; **L** (combined) = 64
- **Paper claim**: NOT ME
- **Kernel result**: NOT ME (ecart_sum = 153)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 32 | 31 | 1 |
| 3 | 21 | 20 | 1 |
| 4 | 16 | 15 | 1 |
| 5 | 12 | 11 | 1 |
| 6 | 10 | 9 | 1 |
| 7 | 9 | 8 | 1 |
| 8 | 8 | 7 | 1 |
| 9 | 7 | 6 | 1 |
| 10 | 6 | 5 | 1 |
| 12 | 5 | 4 | 1 |
| 16 | 4 | 3 | 1 |
| 21 | 3 | 2 | 1 |
| 32 | 2 | 1 | 1 |

## Section 5 — Tables 3, 7–10 (combined (31,32))

### (31,32) s=1

- **id**: `bhuv2026-t3-7-10-s1`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 1; **L** (combined) = 31
- **Paper claim**: NOT ME
- **Kernel result**: NOT ME (ecart_sum = 119)

| t | ℓ*_t | Λ_t (kernel) | kernel equi? | paper equi? | match |
|---:|---:|---:|:---|:---|:---:|
| 2 | 31 | 29 | NOT | NOT | ✓ |
| 3 | 21 | 19 | NOT | NOT | ✓ |
| 4 | 15 | 13 | NOT | NOT | ✓ |
| 5 | 12 | 10 | NOT | NOT | ✓ |
| 6 | 10 | 8 | NOT | NOT | ✓ |
| 7 | 9 | 7 | NOT | NOT | ✓ |
| 9 | 7 | 5 | NOT | NOT | ✓ |
| 10 | 6 | 4 | NOT | NOT | ✓ |
| 12 | 5 | 3 | NOT | NOT | ✓ |
| 15 | 4 | 2 | NOT | NOT | ✓ |
| 21 | 3 | 1 | NOT | NOT | ✓ |
| 31 | 2 | 1 | NOT | equi | ✗ |
| 63 | 1 | 0 | equi | equi | ✓ |

### (31,32) s=2

- **id**: `bhuv2026-t3-7-10-s2`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 2; **L** (combined) = 31
- **Paper claim**: NOT ME
- **Kernel result**: NOT ME (ecart_sum = 86)

| t | ℓ*_t | Λ_t (kernel) | kernel equi? | paper equi? | match |
|---:|---:|---:|:---|:---|:---:|
| 2 | 31 | 27 | NOT | NOT | ✓ |
| 3 | 21 | 17 | NOT | NOT | ✓ |
| 4 | 15 | 11 | NOT | NOT | ✓ |
| 5 | 12 | 8 | NOT | NOT | ✓ |
| 6 | 10 | 6 | NOT | NOT | ✓ |
| 7 | 9 | 5 | NOT | NOT | ✓ |
| 9 | 7 | 3 | NOT | NOT | ✓ |
| 10 | 6 | 2 | NOT | NOT | ✓ |
| 12 | 5 | 1 | NOT | NOT | ✓ |
| 15 | 4 | 1 | NOT | equi | ✗ |
| 21 | 3 | 1 | NOT | equi | ✗ |
| 31 | 2 | 0 | equi | equi | ✓ |
| 63 | 1 | 0 | equi | equi | ✓ |

### (31,32) s=4

- **id**: `bhuv2026-t3-7-10-s4`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 4; **L** (combined) = 31
- **Paper claim**: NOT ME
- **Kernel result**: NOT ME (ecart_sum = 53)

| t | ℓ*_t | Λ_t (kernel) | kernel equi? | paper equi? | match |
|---:|---:|---:|:---|:---|:---:|
| 2 | 31 | 23 | NOT | NOT | ✓ |
| 3 | 21 | 13 | NOT | NOT | ✓ |
| 4 | 15 | 7 | NOT | NOT | ✓ |
| 5 | 12 | 4 | NOT | NOT | ✓ |
| 6 | 10 | 2 | NOT | NOT | ✓ |
| 7 | 9 | 2 | NOT | equi | ✗ |
| 9 | 7 | 1 | NOT | NOT | ✓ |
| 10 | 6 | 0 | equi | NOT | ✗ |
| 12 | 5 | 0 | equi | equi | ✓ |
| 15 | 4 | 0 | equi | equi | ✓ |
| 21 | 3 | 1 | NOT | equi | ✗ |
| 31 | 2 | 0 | equi | equi | ✓ |
| 63 | 1 | 0 | equi | equi | ✓ |

### (31,32) s=7

- **id**: `bhuv2026-t3-7-10-s7`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 7; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 19)

| t | ℓ*_t | Λ_t (kernel) | kernel equi? | paper equi? | match |
|---:|---:|---:|:---|:---|:---:|
| 2 | 31 | 5 | NOT | equi | ✗ |
| 3 | 21 | 6 | NOT | equi | ✗ |
| 4 | 15 | 1 | NOT | equi | ✗ |
| 5 | 12 | 0 | equi | equi | ✓ |
| 6 | 10 | 1 | NOT | equi | ✗ |
| 7 | 9 | 1 | NOT | equi | ✗ |
| 9 | 7 | 1 | NOT | equi | ✗ |
| 10 | 6 | 1 | NOT | equi | ✗ |
| 12 | 5 | 0 | equi | equi | ✓ |
| 15 | 4 | 0 | equi | equi | ✓ |
| 21 | 3 | 1 | NOT | equi | ✗ |
| 31 | 2 | 1 | NOT | equi | ✗ |
| 63 | 1 | 0 | equi | equi | ✓ |

### (31,32) s=8

- **id**: `bhuv2026-t3-7-10-s8`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 8; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

| t | ℓ*_t | Λ_t (kernel) | kernel equi? | paper equi? | match |
|---:|---:|---:|:---|:---|:---:|
| 2 | 31 | 3 | NOT | equi | ✗ |
| 3 | 21 | 2 | NOT | equi | ✗ |
| 4 | 15 | 0 | equi | equi | ✓ |
| 5 | 12 | 0 | equi | equi | ✓ |
| 6 | 10 | 0 | equi | equi | ✓ |
| 7 | 9 | 2 | NOT | equi | ✗ |
| 9 | 7 | 1 | NOT | equi | ✗ |
| 10 | 6 | 1 | NOT | equi | ✗ |
| 12 | 5 | 0 | equi | equi | ✓ |
| 15 | 4 | 0 | equi | equi | ✓ |
| 21 | 3 | 1 | NOT | equi | ✗ |
| 31 | 2 | 1 | NOT | equi | ✗ |
| 63 | 1 | 0 | equi | equi | ✓ |

## Section 6.1 — Table 12 (49 Cattell-Zhang combined CAs)

### (31,32,s=7) — Table 12 [max]

- **id**: `bhuv2026-t12-k31-k32-s7`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 7; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 19)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 5 | 26 |
| 3 | 21 | 6 | 15 |
| 4 | 15 | 1 | 14 |
| 6 | 10 | 1 | 9 |
| 7 | 9 | 1 | 8 |
| 9 | 7 | 1 | 6 |
| 10 | 6 | 1 | 5 |
| 21 | 3 | 1 | 2 |
| 31 | 2 | 1 | 1 |

### (31,32,s=8) — Table 12 [max]

- **id**: `bhuv2026-t12-k31-k32-s8`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 8; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 3 | 28 |
| 3 | 21 | 2 | 19 |
| 7 | 9 | 2 | 7 |
| 9 | 7 | 1 | 6 |
| 10 | 6 | 1 | 5 |
| 21 | 3 | 1 | 2 |
| 31 | 2 | 1 | 1 |

### (31,40,s=8) — Table 12 [max]

- **id**: `bhuv2026-t12-k31-k40-s8`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 40; rule-150 positions: `[7]`
- **k_g** = 71
- **s** = 8; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 9)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 5 | 26 |
| 5 | 14 | 1 | 13 |
| 7 | 10 | 1 | 9 |
| 10 | 7 | 1 | 6 |
| 23 | 3 | 1 | 2 |

### (35,48,s=8) — Table 12 [max]

- **id**: `bhuv2026-t12-k35-k48-s8`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 83
- **s** = 8; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 22)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 35 | 12 | 23 |
| 3 | 27 | 4 | 23 |
| 4 | 20 | 1 | 19 |
| 9 | 9 | 1 | 8 |
| 41 | 2 | 1 | 1 |

### (41,48,s=8) — Table 12 [max]

- **id**: `bhuv2026-t12-k41-k48-s8`
- **family**: combined
- **k1** = 41; rule-150 positions: `[0]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 89
- **s** = 8; **L** (combined) = 41
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 32)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 41 | 18 | 23 |
| 3 | 29 | 6 | 23 |
| 4 | 22 | 4 | 18 |
| 7 | 12 | 1 | 11 |
| 11 | 8 | 1 | 7 |
| 22 | 4 | 1 | 3 |
| 44 | 2 | 1 | 1 |

### (43,48,s=8) — Table 12 [max]

- **id**: `bhuv2026-t12-k43-k48-s8`
- **family**: combined
- **k1** = 43; rule-150 positions: `[2]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 91
- **s** = 8; **L** (combined) = 43
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 45)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 43 | 20 | 23 |
| 3 | 30 | 7 | 23 |
| 4 | 22 | 8 | 14 |
| 5 | 18 | 4 | 14 |
| 6 | 15 | 1 | 14 |
| 7 | 13 | 1 | 12 |
| 9 | 10 | 1 | 9 |
| 18 | 5 | 1 | 4 |

### (47,56,s=8) — Table 12 [max]

- **id**: `bhuv2026-t12-k47-k56-s8`
- **family**: combined
- **k1** = 47; rule-150 positions: `[12]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 103
- **s** = 8; **L** (combined) = 47
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 41)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 47 | 18 | 29 |
| 3 | 34 | 11 | 23 |
| 4 | 25 | 7 | 18 |
| 5 | 20 | 2 | 18 |
| 6 | 17 | 1 | 16 |

### (31,32,s=5) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k31-k32-s5`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 5; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 43)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 21 | 10 |
| 3 | 21 | 11 | 10 |
| 4 | 15 | 5 | 10 |
| 5 | 12 | 2 | 10 |
| 31 | 2 | 1 | 1 |

### (31,32,s=6) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k31-k32-s6`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 6; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 27)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 9 | 22 |
| 3 | 21 | 6 | 15 |
| 4 | 15 | 3 | 12 |
| 7 | 9 | 1 | 8 |
| 9 | 7 | 1 | 6 |
| 21 | 3 | 1 | 2 |
| 31 | 2 | 1 | 1 |

### (31,32,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k31-k32-s9`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 9; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 6)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 1 | 30 |
| 3 | 21 | 1 | 20 |
| 9 | 7 | 1 | 6 |
| 21 | 3 | 1 | 2 |

### (31,32,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k31-k32-s10`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 32; rule-150 positions: `[0, 14]`
- **k_g** = 63
- **s** = 10; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 2 | 29 |
| 3 | 21 | 1 | 20 |
| 6 | 10 | 1 | 9 |
| 7 | 9 | 1 | 8 |
| 9 | 7 | 1 | 6 |
| 10 | 6 | 1 | 5 |
| 15 | 4 | 1 | 3 |

### (31,40,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k31-k40-s9`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 40; rule-150 positions: `[7]`
- **k_g** = 71
- **s** = 9; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 9)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 3 | 28 |
| 3 | 23 | 1 | 22 |
| 5 | 14 | 1 | 13 |
| 10 | 7 | 1 | 6 |
| 23 | 3 | 1 | 2 |
| 35 | 2 | 1 | 1 |

### (31,40,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k31-k40-s10`
- **family**: combined
- **k1** = 31; rule-150 positions: `[10]`
- **k2** = 40; rule-150 positions: `[7]`
- **k_g** = 71
- **s** = 10; **L** (combined) = 31
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 4)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 31 | 1 | 30 |
| 5 | 14 | 1 | 13 |
| 10 | 7 | 1 | 6 |
| 35 | 2 | 1 | 1 |

### (33,40,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k33-k40-s9`
- **family**: combined
- **k1** = 33; rule-150 positions: `[0]`
- **k2** = 40; rule-150 positions: `[7]`
- **k_g** = 73
- **s** = 9; **L** (combined) = 33
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 14)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 33 | 8 | 25 |
| 4 | 18 | 1 | 17 |
| 6 | 12 | 1 | 11 |
| 8 | 9 | 1 | 8 |
| 9 | 8 | 1 | 7 |
| 12 | 6 | 1 | 5 |
| 24 | 3 | 1 | 2 |

### (33,40,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k33-k40-s10`
- **family**: combined
- **k1** = 33; rule-150 positions: `[0]`
- **k2** = 40; rule-150 positions: `[7]`
- **k_g** = 73
- **s** = 10; **L** (combined) = 33
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 13)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 33 | 6 | 27 |
| 3 | 24 | 1 | 23 |
| 4 | 18 | 1 | 17 |
| 6 | 12 | 1 | 11 |
| 8 | 9 | 1 | 8 |
| 18 | 4 | 1 | 3 |
| 36 | 2 | 1 | 1 |

### (33,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k33-k56-s10`
- **family**: combined
- **k1** = 33; rule-150 positions: `[0]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 89
- **s** = 10; **L** (combined) = 33
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 5)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 22 | 1 | 21 |
| 8 | 11 | 1 | 10 |
| 11 | 8 | 1 | 7 |
| 12 | 7 | 1 | 6 |
| 44 | 2 | 1 | 1 |

### (35,48,s=7) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k35-k48-s7`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 83
- **s** = 7; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 28)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 35 | 13 | 22 |
| 3 | 27 | 5 | 22 |
| 4 | 20 | 6 | 14 |
| 5 | 16 | 2 | 14 |
| 10 | 8 | 1 | 7 |
| 41 | 2 | 1 | 1 |

### (35,48,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k35-k48-s9`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 83
- **s** = 9; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 17)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 35 | 11 | 24 |
| 3 | 27 | 3 | 24 |
| 8 | 10 | 1 | 9 |
| 27 | 3 | 1 | 2 |

### (35,48,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k35-k48-s10`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 83
- **s** = 10; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 13)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 35 | 10 | 25 |
| 3 | 27 | 2 | 25 |
| 9 | 9 | 1 | 8 |

### (35,64,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k35-k64-s10`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 99
- **s** = 10; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 13)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 24 | 7 | 17 |
| 5 | 19 | 2 | 17 |
| 7 | 14 | 1 | 13 |
| 11 | 9 | 1 | 8 |
| 12 | 8 | 1 | 7 |
| 33 | 3 | 1 | 2 |

### (39,40,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k39-k40-s9`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 40; rule-150 positions: `[7]`
- **k_g** = 79
- **s** = 9; **L** (combined) = 39
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 21)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 39 | 14 | 25 |
| 3 | 26 | 5 | 21 |
| 4 | 19 | 1 | 18 |
| 26 | 3 | 1 | 2 |

### (39,40,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k39-k40-s10`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 40; rule-150 positions: `[7]`
- **k_g** = 79
- **s** = 10; **L** (combined) = 39
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 22)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 39 | 12 | 27 |
| 3 | 26 | 4 | 22 |
| 6 | 13 | 1 | 12 |
| 7 | 11 | 1 | 10 |
| 26 | 3 | 1 | 2 |

### (39,56,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k39-k56-s9`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 95
- **s** = 9; **L** (combined) = 39
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 14)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 39 | 8 | 31 |
| 5 | 19 | 1 | 18 |
| 19 | 5 | 1 | 4 |
| 47 | 2 | 1 | 1 |

### (39,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k39-k56-s10`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 95
- **s** = 10; **L** (combined) = 39
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 9)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 39 | 6 | 33 |
| 19 | 5 | 1 | 4 |
| 23 | 4 | 1 | 3 |
| 31 | 3 | 1 | 2 |

### (41,48,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k41-k48-s10`
- **family**: combined
- **k1** = 41; rule-150 positions: `[0]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 89
- **s** = 10; **L** (combined) = 41
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 24)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 41 | 16 | 25 |
| 3 | 29 | 4 | 25 |
| 44 | 2 | 1 | 1 |

### (41,64,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k41-k64-s10`
- **family**: combined
- **k1** = 41; rule-150 positions: `[0]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 105
- **s** = 10; **L** (combined) = 41
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 43)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 41 | 8 | 33 |
| 3 | 35 | 8 | 27 |
| 4 | 26 | 9 | 17 |
| 5 | 21 | 4 | 17 |
| 8 | 13 | 1 | 12 |
| 15 | 7 | 1 | 6 |
| 21 | 5 | 1 | 4 |
| 26 | 4 | 1 | 3 |
| 35 | 3 | 1 | 2 |
| 52 | 2 | 1 | 1 |

### (43,48,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k43-k48-s9`
- **family**: combined
- **k1** = 43; rule-150 positions: `[2]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 91
- **s** = 9; **L** (combined) = 43
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 32)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 43 | 19 | 24 |
| 3 | 30 | 6 | 24 |
| 7 | 13 | 1 | 12 |
| 9 | 10 | 1 | 9 |
| 13 | 7 | 1 | 6 |
| 15 | 6 | 1 | 5 |
| 18 | 5 | 1 | 4 |

### (43,48,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k43-k48-s10`
- **family**: combined
- **k1** = 43; rule-150 positions: `[2]`
- **k2** = 48; rule-150 positions: `[14]`
- **k_g** = 91
- **s** = 10; **L** (combined) = 43
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 31)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 43 | 18 | 25 |
| 3 | 30 | 5 | 25 |
| 7 | 13 | 1 | 12 |
| 8 | 11 | 1 | 10 |
| 10 | 9 | 1 | 8 |
| 13 | 7 | 1 | 6 |
| 15 | 6 | 1 | 5 |
| 18 | 5 | 1 | 4 |
| 30 | 3 | 1 | 2 |
| 45 | 2 | 1 | 1 |

### (43,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k43-k56-s10`
- **family**: combined
- **k1** = 43; rule-150 positions: `[2]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 99
- **s** = 10; **L** (combined) = 43
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 15)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 43 | 10 | 33 |
| 3 | 33 | 1 | 32 |
| 9 | 11 | 1 | 10 |
| 11 | 9 | 1 | 8 |
| 14 | 7 | 1 | 6 |
| 49 | 2 | 1 | 1 |

### (45,56,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k45-k56-s9`
- **family**: combined
- **k1** = 45; rule-150 positions: `[8]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 101
- **s** = 9; **L** (combined) = 45
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 23)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 45 | 14 | 31 |
| 3 | 33 | 5 | 28 |
| 4 | 25 | 1 | 24 |
| 5 | 20 | 1 | 19 |
| 50 | 2 | 1 | 1 |

### (45,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k45-k56-s10`
- **family**: combined
- **k1** = 45; rule-150 positions: `[8]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 101
- **s** = 10; **L** (combined) = 45
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 22)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 45 | 12 | 33 |
| 3 | 33 | 4 | 29 |
| 4 | 25 | 1 | 24 |
| 10 | 10 | 1 | 9 |
| 11 | 9 | 1 | 8 |
| 20 | 5 | 1 | 4 |
| 33 | 3 | 1 | 2 |

### (45,64,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k45-k64-s9`
- **family**: combined
- **k1** = 45; rule-150 positions: `[8]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 109
- **s** = 9; **L** (combined) = 45
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 33)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 45 | 17 | 28 |
| 3 | 36 | 8 | 28 |
| 4 | 27 | 1 | 26 |
| 5 | 21 | 1 | 20 |
| 6 | 18 | 1 | 17 |
| 18 | 6 | 1 | 5 |
| 27 | 4 | 1 | 3 |
| 54 | 2 | 1 | 1 |

### (45,64,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k45-k64-s10`
- **family**: combined
- **k1** = 45; rule-150 positions: `[8]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 109
- **s** = 10; **L** (combined) = 45
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 28)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 45 | 16 | 29 |
| 3 | 36 | 7 | 29 |
| 4 | 27 | 1 | 26 |
| 9 | 12 | 1 | 11 |
| 15 | 7 | 1 | 6 |
| 27 | 4 | 1 | 3 |
| 54 | 2 | 1 | 1 |

### (47,56,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k47-k56-s9`
- **family**: combined
- **k1** = 47; rule-150 positions: `[12]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 103
- **s** = 9; **L** (combined) = 47
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 37)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 47 | 16 | 31 |
| 3 | 34 | 7 | 27 |
| 4 | 25 | 7 | 18 |
| 5 | 20 | 2 | 18 |
| 10 | 10 | 1 | 9 |
| 34 | 3 | 1 | 2 |
| 51 | 2 | 1 | 1 |

### (47,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k47-k56-s10`
- **family**: combined
- **k1** = 47; rule-150 positions: `[12]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 103
- **s** = 10; **L** (combined) = 47
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 25)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 47 | 14 | 33 |
| 3 | 34 | 5 | 29 |
| 4 | 25 | 2 | 23 |
| 6 | 17 | 1 | 16 |
| 11 | 9 | 1 | 8 |
| 17 | 6 | 1 | 5 |
| 34 | 3 | 1 | 2 |

### (47,64,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k47-k64-s9`
- **family**: combined
- **k1** = 47; rule-150 positions: `[12]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 111
- **s** = 9; **L** (combined) = 47
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 55)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 47 | 17 | 30 |
| 3 | 37 | 7 | 30 |
| 4 | 27 | 14 | 13 |
| 5 | 22 | 9 | 13 |
| 6 | 18 | 5 | 13 |
| 7 | 15 | 2 | 13 |
| 37 | 3 | 1 | 2 |

### (47,64,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k47-k64-s10`
- **family**: combined
- **k1** = 47; rule-150 positions: `[12]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 111
- **s** = 10; **L** (combined) = 47
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 26)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 47 | 15 | 32 |
| 3 | 37 | 5 | 32 |
| 11 | 10 | 1 | 9 |
| 12 | 9 | 1 | 8 |
| 37 | 3 | 1 | 2 |

### (47,72,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k47-k72-s10`
- **family**: combined
- **k1** = 47; rule-150 positions: `[12]`
- **k2** = 72; rule-150 positions: `[5, 54]`
- **k_g** = 119
- **s** = 10; **L** (combined) = 47
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 28)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 47 | 15 | 32 |
| 3 | 39 | 7 | 32 |
| 4 | 29 | 2 | 27 |
| 5 | 23 | 2 | 21 |
| 17 | 7 | 1 | 6 |

### (51,56,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k51-k56-s9`
- **family**: combined
- **k1** = 51; rule-150 positions: `[0]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 107
- **s** = 9; **L** (combined) = 51
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 49)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 51 | 28 | 23 |
| 3 | 35 | 12 | 23 |
| 4 | 26 | 4 | 22 |
| 5 | 21 | 1 | 20 |

### (51,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k51-k56-s10`
- **family**: combined
- **k1** = 51; rule-150 positions: `[0]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 107
- **s** = 10; **L** (combined) = 51
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 28)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 51 | 18 | 33 |
| 3 | 35 | 7 | 28 |
| 4 | 26 | 1 | 25 |
| 17 | 6 | 1 | 5 |
| 53 | 2 | 1 | 1 |

### (53,56,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k53-k56-s9`
- **family**: combined
- **k1** = 53; rule-150 positions: `[0]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 109
- **s** = 9; **L** (combined) = 53
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 61)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 53 | 30 | 23 |
| 3 | 36 | 13 | 23 |
| 4 | 27 | 5 | 22 |
| 5 | 21 | 1 | 20 |
| 6 | 18 | 2 | 16 |
| 9 | 12 | 1 | 11 |
| 12 | 9 | 1 | 8 |
| 18 | 6 | 1 | 5 |
| 27 | 4 | 1 | 3 |
| 36 | 3 | 1 | 2 |

### (53,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k53-k56-s10`
- **family**: combined
- **k1** = 53; rule-150 positions: `[0]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 109
- **s** = 10; **L** (combined) = 53
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 36)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 53 | 20 | 33 |
| 3 | 36 | 9 | 27 |
| 4 | 27 | 3 | 24 |
| 9 | 12 | 1 | 11 |
| 27 | 4 | 1 | 3 |
| 36 | 3 | 1 | 2 |

### (55,56,s=9) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k55-k56-s9`
- **family**: combined
- **k1** = 55; rule-150 positions: `[16]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 111
- **s** = 9; **L** (combined) = 55
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 55)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 55 | 23 | 32 |
| 3 | 37 | 17 | 20 |
| 4 | 27 | 8 | 19 |
| 5 | 22 | 3 | 19 |
| 6 | 18 | 1 | 17 |
| 11 | 10 | 1 | 9 |
| 18 | 6 | 1 | 5 |
| 37 | 3 | 1 | 2 |

### (55,56,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k55-k56-s10`
- **family**: combined
- **k1** = 55; rule-150 positions: `[16]`
- **k2** = 56; rule-150 positions: `[3, 13]`
- **k_g** = 111
- **s** = 10; **L** (combined) = 55
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 46)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 55 | 22 | 33 |
| 3 | 37 | 15 | 22 |
| 4 | 27 | 6 | 21 |
| 5 | 22 | 1 | 21 |
| 22 | 5 | 1 | 4 |

### (59,64,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k59-k64-s10`
- **family**: combined
- **k1** = 59; rule-150 positions: `[3, 14]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 123
- **s** = 10; **L** (combined) = 59
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 46)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 59 | 25 | 34 |
| 3 | 41 | 12 | 29 |
| 4 | 30 | 5 | 25 |
| 5 | 24 | 1 | 23 |
| 24 | 5 | 1 | 4 |
| 41 | 3 | 1 | 2 |
| 61 | 2 | 1 | 1 |

### (63,64,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k63-k64-s10`
- **family**: combined
- **k1** = 63; rule-150 positions: `[30]`
- **k2** = 64; rule-150 positions: `[2, 4]`
- **k_g** = 127
- **s** = 10; **L** (combined) = 63
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 79)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 63 | 39 | 24 |
| 3 | 42 | 21 | 21 |
| 4 | 31 | 11 | 20 |
| 5 | 25 | 5 | 20 |
| 6 | 21 | 1 | 20 |
| 7 | 18 | 1 | 17 |
| 21 | 6 | 1 | 5 |

### (63,80,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k63-k80-s10`
- **family**: combined
- **k1** = 63; rule-150 positions: `[30]`
- **k2** = 80; rule-150 positions: `[0, 70]`
- **k_g** = 143
- **s** = 10; **L** (combined) = 63
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 101)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 63 | 43 | 20 |
| 3 | 47 | 27 | 20 |
| 4 | 35 | 15 | 20 |
| 5 | 28 | 8 | 20 |
| 6 | 23 | 3 | 20 |
| 8 | 17 | 2 | 15 |
| 11 | 13 | 1 | 12 |
| 13 | 11 | 1 | 10 |
| 28 | 5 | 1 | 4 |

### (67,72,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k67-k72-s10`
- **family**: combined
- **k1** = 67; rule-150 positions: `[14]`
- **k2** = 72; rule-150 positions: `[5, 54]`
- **k_g** = 139
- **s** = 10; **L** (combined) = 64
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 80)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 64 | 36 | 28 |
| 3 | 46 | 21 | 25 |
| 4 | 34 | 13 | 21 |
| 5 | 27 | 6 | 21 |
| 6 | 23 | 2 | 21 |
| 23 | 6 | 1 | 5 |
| 46 | 3 | 1 | 2 |

### (71,72,s=10) — Table 12 [reduced]

- **id**: `bhuv2026-t12-k71-k72-s10`
- **family**: combined
- **k1** = 71; rule-150 positions: `[16]`
- **k2** = 72; rule-150 positions: `[5, 54]`
- **k_g** = 143
- **s** = 10; **L** (combined) = 64
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 125)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 64 | 48 | 16 |
| 3 | 47 | 31 | 16 |
| 4 | 35 | 19 | 16 |
| 5 | 28 | 12 | 16 |
| 6 | 23 | 7 | 16 |
| 7 | 20 | 4 | 16 |
| 8 | 17 | 1 | 16 |
| 11 | 13 | 1 | 12 |
| 13 | 11 | 1 | 10 |
| 47 | 3 | 1 | 2 |

## Section 6.2 — Table 11 (Adak-Das CA(90′)/CA(150′) combos)

### CA(90′) k=26 × CA(90′) k=29, s=10

- **id**: `bhuv2026-t11-90prime26-90prime29-s10`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[0]`
- **k_g** = 55
- **s** = 10; **L** (combined) = 26
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 35)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 26 | 13 | 13 |
| 3 | 18 | 10 | 8 |
| 4 | 13 | 5 | 8 |
| 5 | 11 | 3 | 8 |
| 6 | 9 | 1 | 8 |
| 7 | 7 | 1 | 6 |

### CA(90′) k=26 × CA(90′) k=35, s=10

- **id**: `bhuv2026-t11-90prime26-90prime35-s10`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[0]`
- **k_g** = 61
- **s** = 10; **L** (combined) = 26
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 23)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 26 | 8 | 18 |
| 3 | 20 | 6 | 14 |
| 4 | 15 | 5 | 10 |
| 5 | 12 | 2 | 10 |
| 12 | 5 | 1 | 4 |
| 20 | 3 | 1 | 2 |

### CA(90′) k=29 × CA(90′) k=35, s=10

- **id**: `bhuv2026-t11-90prime29-90prime35-s10`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[0]`
- **k_g** = 64
- **s** = 10; **L** (combined) = 29
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 35)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 29 | 13 | 16 |
| 3 | 21 | 9 | 12 |
| 4 | 16 | 6 | 10 |
| 5 | 12 | 2 | 10 |
| 7 | 9 | 1 | 8 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 16 | 4 | 1 | 3 |

### CA(150′) k=26 × CA(150′) k=29, s=10

- **id**: `bhuv2026-t11-150prime26-150prime29-s10`
- **family**: combined
- **k1** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 55
- **s** = 10; **L** (combined) = 26
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 32)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 26 | 13 | 13 |
| 3 | 18 | 8 | 10 |
| 4 | 13 | 3 | 10 |
| 5 | 11 | 1 | 10 |
| 6 | 9 | 1 | 8 |
| 9 | 6 | 1 | 5 |
| 11 | 5 | 1 | 4 |

### CA(150′) k=26 × CA(150′) k=35, s=10

- **id**: `bhuv2026-t11-150prime26-150prime35-s10`
- **family**: combined
- **k1** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 61
- **s** = 10; **L** (combined) = 26
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 26)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 26 | 10 | 16 |
| 3 | 20 | 6 | 14 |
| 4 | 15 | 5 | 10 |
| 5 | 12 | 2 | 10 |
| 10 | 6 | 1 | 5 |
| 12 | 5 | 1 | 4 |

### CA(150′) k=29 × CA(150′) k=35, s=10

- **id**: `bhuv2026-t11-150prime29-150prime35-s10`
- **family**: combined
- **k1** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 64
- **s** = 10; **L** (combined) = 29
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 32)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 29 | 13 | 16 |
| 3 | 21 | 10 | 11 |
| 4 | 16 | 6 | 10 |
| 5 | 12 | 2 | 10 |
| 16 | 4 | 1 | 3 |

### CA(90′) k=26 × CA(150′) k=29, s=5

- **id**: `bhuv2026-t11-90prime26-150prime29-s5`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 55
- **s** = 5; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 16)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 18 | 5 | 13 |
| 4 | 13 | 3 | 10 |
| 5 | 11 | 2 | 9 |
| 6 | 9 | 1 | 8 |
| 9 | 6 | 1 | 5 |
| 11 | 5 | 1 | 4 |
| 18 | 3 | 1 | 2 |

### CA(90′) k=29 × CA(150′) k=26, s=5

- **id**: `bhuv2026-t11-90prime29-150prime26-s5`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 55
- **s** = 5; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 18 | 5 | 13 |
| 4 | 13 | 3 | 10 |
| 5 | 11 | 2 | 9 |
| 9 | 6 | 1 | 5 |

### CA(90′) k=26 × CA(150′) k=29, s=7

- **id**: `bhuv2026-t11-90prime26-150prime29-s7`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 55
- **s** = 7; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 5)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 18 | 1 | 17 |
| 5 | 11 | 1 | 10 |
| 9 | 6 | 1 | 5 |
| 11 | 5 | 1 | 4 |
| 18 | 3 | 1 | 2 |

### CA(90′) k=29 × CA(150′) k=26, s=7

- **id**: `bhuv2026-t11-90prime29-150prime26-s7`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 55
- **s** = 7; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 6)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 26 | 1 | 25 |
| 3 | 18 | 1 | 17 |
| 5 | 11 | 1 | 10 |
| 9 | 6 | 1 | 5 |
| 18 | 3 | 1 | 2 |
| 27 | 2 | 1 | 1 |

### CA(90′) k=26 × CA(150′) k=29, s=8

- **id**: `bhuv2026-t11-90prime26-150prime29-s8`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 55
- **s** = 8; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 3)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 5 | 11 | 1 | 10 |
| 6 | 9 | 1 | 8 |
| 18 | 3 | 1 | 2 |

### CA(90′) k=29 × CA(150′) k=26, s=8

- **id**: `bhuv2026-t11-90prime29-150prime26-s8`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 55
- **s** = 8; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 5)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 5 | 11 | 1 | 10 |
| 6 | 9 | 1 | 8 |
| 13 | 4 | 1 | 3 |
| 18 | 3 | 1 | 2 |

### CA(90′) k=26 × CA(150′) k=29, s=10

- **id**: `bhuv2026-t11-90prime26-150prime29-s10`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 55
- **s** = 10; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 3)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 18 | 1 | 17 |
| 6 | 9 | 1 | 8 |
| 9 | 6 | 1 | 5 |

### CA(90′) k=29 × CA(150′) k=26, s=10

- **id**: `bhuv2026-t11-90prime29-150prime26-s10`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 55
- **s** = 10; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 2)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 9 | 6 | 1 | 5 |
| 11 | 5 | 1 | 4 |

### CA(90′) k=26 × CA(150′) k=35, s=5

- **id**: `bhuv2026-t11-90prime26-150prime35-s5`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 61
- **s** = 5; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 10)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 20 | 2 | 18 |
| 4 | 15 | 2 | 13 |
| 5 | 12 | 2 | 10 |
| 10 | 6 | 1 | 5 |
| 15 | 4 | 1 | 3 |
| 20 | 3 | 1 | 2 |

### CA(90′) k=35 × CA(150′) k=26, s=5

- **id**: `bhuv2026-t11-90prime35-150prime26-s5`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 61
- **s** = 5; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 14)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 20 | 3 | 17 |
| 4 | 15 | 2 | 13 |
| 5 | 12 | 2 | 10 |
| 6 | 10 | 1 | 9 |
| 10 | 6 | 1 | 5 |
| 15 | 4 | 1 | 3 |
| 30 | 2 | 1 | 1 |

### CA(90′) k=26 × CA(150′) k=35, s=7

- **id**: `bhuv2026-t11-90prime26-150prime35-s7`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 61
- **s** = 7; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 1)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 30 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=26, s=7

- **id**: `bhuv2026-t11-90prime35-150prime26-s7`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 61
- **s** = 7; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 6)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 20 | 1 | 19 |
| 5 | 12 | 1 | 11 |
| 6 | 10 | 1 | 9 |
| 10 | 6 | 1 | 5 |
| 12 | 5 | 1 | 4 |
| 15 | 4 | 1 | 3 |

### CA(90′) k=26 × CA(150′) k=35, s=8

- **id**: `bhuv2026-t11-90prime26-150prime35-s8`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 61
- **s** = 8; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 8)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 20 | 2 | 18 |
| 4 | 15 | 1 | 14 |
| 5 | 12 | 1 | 11 |
| 6 | 10 | 1 | 9 |
| 10 | 6 | 1 | 5 |
| 15 | 4 | 1 | 3 |

### CA(90′) k=35 × CA(150′) k=26, s=8

- **id**: `bhuv2026-t11-90prime35-150prime26-s8`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 61
- **s** = 8; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 4)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 20 | 2 | 18 |
| 10 | 6 | 1 | 5 |

### CA(90′) k=26 × CA(150′) k=35, s=10

- **id**: `bhuv2026-t11-90prime26-150prime35-s10`
- **family**: combined
- **k1** = 26; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 61
- **s** = 10; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 5)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 20 | 1 | 19 |
| 4 | 15 | 1 | 14 |
| 6 | 10 | 1 | 9 |
| 15 | 4 | 1 | 3 |
| 20 | 3 | 1 | 2 |

### CA(90′) k=35 × CA(150′) k=26, s=10

- **id**: `bhuv2026-t11-90prime35-150prime26-s10`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 26; rule-150 positions: `[1, 2, 3, …(25 positions)…, 23, 24, 25]`
- **k_g** = 61
- **s** = 10; **L** (combined) = 26
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 4)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 6 | 10 | 1 | 9 |
| 10 | 6 | 1 | 5 |
| 12 | 5 | 1 | 4 |
| 20 | 3 | 1 | 2 |

### CA(90′) k=29 × CA(150′) k=35, s=5

- **id**: `bhuv2026-t11-90prime29-150prime35-s5`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 64
- **s** = 5; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 17)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 21 | 7 | 14 |
| 4 | 16 | 3 | 13 |
| 5 | 12 | 2 | 10 |
| 8 | 8 | 1 | 7 |
| 16 | 4 | 1 | 3 |
| 21 | 3 | 1 | 2 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=29, s=5

- **id**: `bhuv2026-t11-90prime35-150prime29-s5`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 64
- **s** = 5; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 16)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 21 | 7 | 14 |
| 4 | 16 | 3 | 13 |
| 5 | 12 | 2 | 10 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 16 | 4 | 1 | 3 |

### CA(90′) k=29 × CA(150′) k=35, s=6

- **id**: `bhuv2026-t11-90prime29-150prime35-s6`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 64
- **s** = 6; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 21 | 3 | 18 |
| 4 | 16 | 2 | 14 |
| 6 | 10 | 1 | 9 |
| 7 | 9 | 1 | 8 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 16 | 4 | 1 | 3 |
| 21 | 3 | 1 | 2 |

### CA(90′) k=35 × CA(150′) k=29, s=6

- **id**: `bhuv2026-t11-90prime35-150prime29-s6`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 64
- **s** = 6; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 21 | 3 | 18 |
| 4 | 16 | 2 | 14 |
| 7 | 9 | 1 | 8 |
| 8 | 8 | 1 | 7 |
| 12 | 5 | 1 | 4 |
| 16 | 4 | 1 | 3 |
| 21 | 3 | 1 | 2 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=29 × CA(150′) k=35, s=7

- **id**: `bhuv2026-t11-90prime29-150prime35-s7`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 64
- **s** = 7; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 9)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 29 | 1 | 28 |
| 3 | 21 | 1 | 20 |
| 4 | 16 | 1 | 15 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=29, s=7

- **id**: `bhuv2026-t11-90prime35-150prime29-s7`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 64
- **s** = 7; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 7)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 21 | 1 | 20 |
| 4 | 16 | 1 | 15 |
| 8 | 8 | 1 | 7 |
| 16 | 4 | 1 | 3 |
| 21 | 3 | 1 | 2 |

### CA(90′) k=29 × CA(150′) k=35, s=8

- **id**: `bhuv2026-t11-90prime29-150prime35-s8`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 64
- **s** = 8; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 16 | 2 | 14 |
| 5 | 12 | 1 | 11 |
| 7 | 9 | 1 | 8 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 16 | 4 | 1 | 3 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=29, s=8

- **id**: `bhuv2026-t11-90prime35-150prime29-s8`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 64
- **s** = 8; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 9)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 16 | 1 | 15 |
| 7 | 9 | 1 | 8 |
| 8 | 8 | 1 | 7 |
| 16 | 4 | 1 | 3 |
| 21 | 3 | 1 | 2 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=29 × CA(150′) k=35, s=9

- **id**: `bhuv2026-t11-90prime29-150prime35-s9`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 64
- **s** = 9; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 6)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 16 | 1 | 15 |
| 6 | 10 | 1 | 9 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 21 | 3 | 1 | 2 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=29, s=9

- **id**: `bhuv2026-t11-90prime35-150prime29-s9`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 64
- **s** = 9; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 3)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 16 | 1 | 15 |
| 8 | 8 | 1 | 7 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=29 × CA(150′) k=35, s=10

- **id**: `bhuv2026-t11-90prime29-150prime35-s10`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 64
- **s** = 10; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 7)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 7 | 9 | 1 | 8 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 16 | 4 | 1 | 3 |
| 21 | 3 | 1 | 2 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=29, s=10

- **id**: `bhuv2026-t11-90prime35-150prime29-s10`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 64
- **s** = 10; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 6)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 16 | 1 | 15 |
| 8 | 8 | 1 | 7 |
| 16 | 4 | 1 | 3 |
| 21 | 3 | 1 | 2 |
| 32 | 2 | 1 | 1 |

### CA(90′) k=29 × CA(150′) k=39, s=6

- **id**: `bhuv2026-t11-90prime29-150prime39-s6`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 68
- **s** = 6; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 5)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 22 | 1 | 21 |
| 5 | 13 | 1 | 12 |
| 22 | 3 | 1 | 2 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=29, s=6

- **id**: `bhuv2026-t11-90prime39-150prime29-s6`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 68
- **s** = 6; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 8)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 22 | 1 | 21 |
| 5 | 13 | 1 | 12 |
| 6 | 11 | 1 | 10 |
| 13 | 5 | 1 | 4 |
| 17 | 4 | 1 | 3 |
| 22 | 3 | 1 | 2 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=29 × CA(150′) k=39, s=8

- **id**: `bhuv2026-t11-90prime29-150prime39-s8`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 68
- **s** = 8; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 11)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 22 | 1 | 21 |
| 4 | 17 | 1 | 16 |
| 5 | 13 | 1 | 12 |
| 17 | 4 | 1 | 3 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=29, s=8

- **id**: `bhuv2026-t11-90prime39-150prime29-s8`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 68
- **s** = 8; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 10)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 22 | 1 | 21 |
| 4 | 17 | 1 | 16 |
| 17 | 4 | 1 | 3 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=29 × CA(150′) k=39, s=9

- **id**: `bhuv2026-t11-90prime29-150prime39-s9`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 68
- **s** = 9; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 8)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 17 | 1 | 16 |
| 6 | 11 | 1 | 10 |
| 8 | 8 | 1 | 7 |
| 9 | 7 | 1 | 6 |
| 17 | 4 | 1 | 3 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=29, s=9

- **id**: `bhuv2026-t11-90prime39-150prime29-s9`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 68
- **s** = 9; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 5)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 17 | 1 | 16 |
| 6 | 11 | 1 | 10 |
| 13 | 5 | 1 | 4 |
| 17 | 4 | 1 | 3 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=29 × CA(150′) k=39, s=10

- **id**: `bhuv2026-t11-90prime29-150prime39-s10`
- **family**: combined
- **k1** = 29; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 68
- **s** = 10; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 5)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 17 | 1 | 16 |
| 5 | 13 | 1 | 12 |
| 17 | 4 | 1 | 3 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=29, s=10

- **id**: `bhuv2026-t11-90prime39-150prime29-s10`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 29; rule-150 positions: `[1, 2, 3, …(28 positions)…, 26, 27, 28]`
- **k_g** = 68
- **s** = 10; **L** (combined) = 29
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 7)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 4 | 17 | 2 | 15 |
| 17 | 4 | 1 | 3 |
| 34 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=39, s=6

- **id**: `bhuv2026-t11-90prime35-150prime39-s6`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 74
- **s** = 6; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 20)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 24 | 11 | 13 |
| 4 | 18 | 5 | 13 |
| 5 | 14 | 2 | 12 |
| 8 | 9 | 1 | 8 |
| 37 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=35, s=6

- **id**: `bhuv2026-t11-90prime39-150prime35-s6`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 74
- **s** = 6; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 18)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 24 | 11 | 13 |
| 4 | 18 | 5 | 13 |
| 5 | 14 | 2 | 12 |

### CA(90′) k=35 × CA(150′) k=39, s=8

- **id**: `bhuv2026-t11-90prime35-150prime39-s8`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 74
- **s** = 8; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 26)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 35 | 4 | 31 |
| 3 | 24 | 5 | 19 |
| 4 | 18 | 3 | 15 |
| 5 | 14 | 2 | 12 |
| 6 | 12 | 1 | 11 |
| 8 | 9 | 1 | 8 |
| 9 | 8 | 1 | 7 |
| 12 | 6 | 1 | 5 |
| 14 | 5 | 1 | 4 |
| 24 | 3 | 1 | 2 |
| 37 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=35, s=8

- **id**: `bhuv2026-t11-90prime39-150prime35-s8`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 74
- **s** = 8; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 24)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 2 | 35 | 4 | 31 |
| 3 | 24 | 5 | 19 |
| 4 | 18 | 3 | 15 |
| 5 | 14 | 1 | 13 |
| 6 | 12 | 1 | 11 |
| 7 | 10 | 1 | 9 |
| 8 | 9 | 1 | 8 |
| 12 | 6 | 1 | 5 |
| 14 | 5 | 1 | 4 |
| 24 | 3 | 1 | 2 |
| 37 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=39, s=9

- **id**: `bhuv2026-t11-90prime35-150prime39-s9`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 74
- **s** = 9; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 10)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 24 | 3 | 21 |
| 5 | 14 | 1 | 13 |
| 9 | 8 | 1 | 7 |
| 18 | 4 | 1 | 3 |
| 37 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=35, s=9

- **id**: `bhuv2026-t11-90prime39-150prime35-s9`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 74
- **s** = 9; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 13)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 24 | 3 | 21 |
| 4 | 18 | 1 | 17 |
| 8 | 9 | 1 | 8 |
| 9 | 8 | 1 | 7 |
| 12 | 6 | 1 | 5 |
| 18 | 4 | 1 | 3 |
| 37 | 2 | 1 | 1 |

### CA(90′) k=35 × CA(150′) k=39, s=10

- **id**: `bhuv2026-t11-90prime35-150prime39-s10`
- **family**: combined
- **k1** = 35; rule-150 positions: `[0]`
- **k2** = 39; rule-150 positions: `[1, 2, 3, …(38 positions)…, 36, 37, 38]`
- **k_g** = 74
- **s** = 10; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 3)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 24 | 1 | 23 |
| 8 | 9 | 1 | 8 |
| 37 | 2 | 1 | 1 |

### CA(90′) k=39 × CA(150′) k=35, s=10

- **id**: `bhuv2026-t11-90prime39-150prime35-s10`
- **family**: combined
- **k1** = 39; rule-150 positions: `[0]`
- **k2** = 35; rule-150 positions: `[1, 2, 3, …(34 positions)…, 32, 33, 34]`
- **k_g** = 74
- **s** = 10; **L** (combined) = 35
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 4)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 24 | 1 | 23 |
| 5 | 14 | 1 | 13 |
| 6 | 12 | 1 | 11 |
| 12 | 6 | 1 | 5 |

### CA(90′) k=65 × CA(150′) k=69, s=10

- **id**: `bhuv2026-t11-90prime65-150prime69-s10`
- **family**: combined
- **k1** = 65; rule-150 positions: `[0]`
- **k2** = 69; rule-150 positions: `[1, 2, 3, …(68 positions)…, 66, 67, 68]`
- **k_g** = 134
- **s** = 10; **L** (combined) = 64
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 48)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 44 | 24 | 20 |
| 4 | 33 | 13 | 20 |
| 5 | 26 | 6 | 20 |
| 6 | 22 | 2 | 20 |
| 7 | 19 | 1 | 18 |
| 12 | 11 | 1 | 10 |
| 67 | 2 | 1 | 1 |

### CA(90′) k=69 × CA(150′) k=65, s=10

- **id**: `bhuv2026-t11-90prime69-150prime65-s10`
- **family**: combined
- **k1** = 69; rule-150 positions: `[0]`
- **k2** = 65; rule-150 positions: `[1, 2, 3, …(64 positions)…, 62, 63, 64]`
- **k_g** = 134
- **s** = 10; **L** (combined) = 64
- **Paper claim**: ME
- **Kernel result**: NOT ME (ecart_sum = 46)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 44 | 24 | 20 |
| 4 | 33 | 13 | 20 |
| 5 | 26 | 6 | 20 |
| 6 | 22 | 2 | 20 |
| 19 | 7 | 1 | 6 |

### CA(90′) k=105 × CA(150′) k=113, s=9

- **id**: `bhuv2026-t11-90prime105-150prime113-s9`
- **family**: combined
- **k1** = 105; rule-150 positions: `[0]`
- **k2** = 113; rule-150 positions: `[1, 2, 3, …(112 positions)…, 110, 111, 112]`
- **k_g** = 218
- **s** = 9; **L** (combined) = 64
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 160)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 64 | 46 | 18 |
| 4 | 54 | 36 | 18 |
| 5 | 43 | 25 | 18 |
| 6 | 36 | 18 | 18 |
| 7 | 31 | 13 | 18 |
| 8 | 27 | 9 | 18 |
| 9 | 24 | 6 | 18 |
| 10 | 21 | 3 | 18 |
| 11 | 19 | 1 | 18 |
| 12 | 18 | 1 | 17 |
| 31 | 7 | 1 | 6 |
| 43 | 5 | 1 | 4 |

### CA(90′) k=113 × CA(150′) k=105, s=9

- **id**: `bhuv2026-t11-90prime113-150prime105-s9`
- **family**: combined
- **k1** = 113; rule-150 positions: `[0]`
- **k2** = 105; rule-150 positions: `[1, 2, 3, …(104 positions)…, 102, 103, 104]`
- **k_g** = 218
- **s** = 9; **L** (combined) = 64
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 161)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 64 | 46 | 18 |
| 4 | 54 | 36 | 18 |
| 5 | 43 | 25 | 18 |
| 6 | 36 | 18 | 18 |
| 7 | 31 | 13 | 18 |
| 8 | 27 | 9 | 18 |
| 9 | 24 | 6 | 18 |
| 10 | 21 | 3 | 18 |
| 11 | 19 | 1 | 18 |
| 18 | 12 | 1 | 11 |
| 31 | 7 | 1 | 6 |
| 36 | 6 | 1 | 5 |
| 109 | 2 | 1 | 1 |

### CA(90′) k=105 × CA(150′) k=113, s=10

- **id**: `bhuv2026-t11-90prime105-150prime113-s10`
- **family**: combined
- **k1** = 105; rule-150 positions: `[0]`
- **k2** = 113; rule-150 positions: `[1, 2, 3, …(112 positions)…, 110, 111, 112]`
- **k_g** = 218
- **s** = 10; **L** (combined) = 64
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 143)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 64 | 44 | 20 |
| 4 | 54 | 34 | 20 |
| 5 | 43 | 23 | 20 |
| 6 | 36 | 16 | 20 |
| 7 | 31 | 11 | 20 |
| 8 | 27 | 7 | 20 |
| 9 | 24 | 4 | 20 |
| 10 | 21 | 1 | 20 |
| 12 | 18 | 1 | 17 |
| 18 | 12 | 1 | 11 |
| 43 | 5 | 1 | 4 |

### CA(90′) k=113 × CA(150′) k=105, s=10

- **id**: `bhuv2026-t11-90prime113-150prime105-s10`
- **family**: combined
- **k1** = 113; rule-150 positions: `[0]`
- **k2** = 105; rule-150 positions: `[1, 2, 3, …(104 positions)…, 102, 103, 104]`
- **k_g** = 218
- **s** = 10; **L** (combined) = 64
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 142)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 64 | 44 | 20 |
| 4 | 54 | 34 | 20 |
| 5 | 43 | 23 | 20 |
| 6 | 36 | 16 | 20 |
| 7 | 31 | 11 | 20 |
| 8 | 27 | 7 | 20 |
| 9 | 24 | 4 | 20 |
| 10 | 21 | 1 | 20 |
| 109 | 2 | 1 | 1 |

### CA(90′) k=113 × CA(150′) k=119, s=10

- **id**: `bhuv2026-t11-90prime113-150prime119-s10`
- **family**: combined
- **k1** = 113; rule-150 positions: `[0]`
- **k2** = 119; rule-150 positions: `[1, 2, 3, …(118 positions)…, 116, 117, 118]`
- **k_g** = 232
- **s** = 10; **L** (combined) = 64
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 160)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 64 | 44 | 20 |
| 4 | 58 | 38 | 20 |
| 5 | 46 | 26 | 20 |
| 6 | 38 | 18 | 20 |
| 7 | 33 | 13 | 20 |
| 8 | 29 | 9 | 20 |
| 9 | 25 | 5 | 20 |
| 10 | 23 | 3 | 20 |
| 11 | 21 | 1 | 20 |
| 29 | 8 | 1 | 7 |
| 46 | 5 | 1 | 4 |
| 58 | 4 | 1 | 3 |

### CA(90′) k=119 × CA(150′) k=113, s=10

- **id**: `bhuv2026-t11-90prime119-150prime113-s10`
- **family**: combined
- **k1** = 119; rule-150 positions: `[0]`
- **k2** = 113; rule-150 positions: `[1, 2, 3, …(112 positions)…, 110, 111, 112]`
- **k_g** = 232
- **s** = 10; **L** (combined) = 64
- **Paper claim**: almost_ME
- **Kernel result**: NOT ME (ecart_sum = 159)

**Deficiencies** (Λ_t > 0):

| t | ℓ*_t | Λ_t | kernel ℓ_t |
|---:|---:|---:|---:|
| 3 | 64 | 44 | 20 |
| 4 | 58 | 38 | 20 |
| 5 | 46 | 26 | 20 |
| 6 | 38 | 18 | 20 |
| 7 | 33 | 13 | 20 |
| 8 | 29 | 9 | 20 |
| 9 | 25 | 5 | 20 |
| 10 | 23 | 3 | 20 |
| 11 | 21 | 1 | 20 |
| 29 | 8 | 1 | 7 |
| 58 | 4 | 1 | 3 |
