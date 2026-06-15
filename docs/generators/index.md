# Generators

REGPOLY implements 17 F₂-linear PRNG families, 2 F₂ digital-net
families, and a `CombinedF2LinearSource` composition wrapper that XORs
J of them. Each family page below documents the recurrence (or matrix
construction, for digital nets), the parameter space, common parameter
sets from the literature, and source-file pointers.

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

## F₂ digital nets

These inherit from the `DigitalNet` base class (a direct
`F2LinearSource` subclass — sibling of `Recurrence` under the post-
refactor hierarchy), so the existing equidistribution machinery runs
on them unchanged. They pair naturally with the
[t-value test](../theory/tvalue-spec.md).

- [SobolNet (Sobol, Joe–Kuo direction numbers)](SobolNet.md)
- [NiederreiterF2Net (Niederreiter 1988 digital sequence in F₂)](NiederreiterF2Net.md)
