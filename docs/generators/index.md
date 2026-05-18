# Generators

REGPOLY implements 17 F₂-linear generator families plus a `CombinedGenerator` that XORs J of them. Each family page below documents the recurrence, the parameter space, common parameter sets from the literature, and source-file pointers.

## Mersenne Twister family

- [MT (Mersenne Twister)](MTGen.md)
- [SFMT (SIMD Fast MT)](SFMTGen.md)
- [dSFMT (Double-precision SFMT)](DSFMTGen.md)
- [MTGP (MT for GPU)](MTGPGen.md)
- [TinyMT32](TinyMT32Gen.md)
- [RMT64 (Reducible MT 64-bit)](RMT64Gen.md)

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
- [Xoroshiro (Blackman–Vigna)](XoroshiroGen.md)
- [Xoshiro (Blackman–Vigna)](XoshiroGen.md)

## Cellular automata

- [Two-rule (90/150) cellular automata](CellularAutomataGen.md)

## MELG

- [MELG (Maximally Equidistributed F₂-Linear)](MELGGen.md)
