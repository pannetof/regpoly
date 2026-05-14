# Generators

REGPOLY implements 18 F₂-linear generator families plus a `CombinedGenerator` that XORs J of them. Each family page below documents the recurrence, the parameter space, common parameter sets from the literature, and source-file pointers.

> **Status**: Phase 0 ports the existing per-family pages as-is. Phase 6 rewrites them onto a common template (`_template.md`, also pending) with paper references, equidistribution properties, and notebook links.

## Mersenne Twister family

- [MT (Mersenne Twister)](MTGen.md)
- [SFMT (SIMD Fast MT)](SFMTGen.md)
- [dSFMT (Double-precision SFMT)](DSFMTGen.md)
- [MTGP (MT for GPU)](MTGPGen.md)
- [TinyMT32](TinyMT32Gen.md)
- [RMT64 (Rotated MT 64-bit)](RMT64Gen.md)
- [Matsumoto](MatsumotoGen.md)

## WELL and Tausworthe-style

- [WELL](WELLGen.md)
- [Tausworthe](TauswortheGen.md)
- [TGFSR (Twisted GFSR)](TGFSRGen.md)
- [PolyLCG](PolyLCGGen.md)

## GF(2^w) extension-field generators

- [F2w-LFSR](F2wLFSRGen.md)
- [F2w-PolyLCG](F2wPolyLCGGen.md)

## Xorshift family

- [Marsaglia XORShift](MarsaXorshiftGen.md)

## Other

- [AC1D (Additive Cellular 1D)](AC1DGen.md)
- [MELG (Multiplicative-Equidistribution LCG)](MELGGen.md)

## Combined generator

- *Page to be authored in Phase 1.*
