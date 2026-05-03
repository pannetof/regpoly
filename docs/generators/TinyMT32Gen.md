# TinyMT32Gen — Tiny Mersenne Twister (32-bit)

`TinyMT32Gen` is the 127-bit-state Tiny MT variant of Saito and
Matsumoto (2011): a compact MT-style F₂-linear generator with state
sized for embedded contexts and parallel-stream applications, yet
retaining an $2^{127}-1$ period.

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `mat1` | uint32 | First per-stream parameter (matrix entry). |
| `mat2` | uint32 | Second per-stream parameter. |
| `tmat` | uint32 | Tempering matrix entry. |

The structural parameters (state size, shift constants) are fixed
constants of the family; `(mat1, mat2, tmat)` is the per-stream tuple
that distinguishes parallel TinyMT instances.

## Provenance

- **Reference**: Saito, M., Matsumoto, M., "TinyMT: a small-sized
  variant of Mersenne Twister", 2011 (released alongside the MTGP paper).
- **Cross-check fixtures**: `docs/library/tinymt_params.yaml`.
- **Library entry** (planned): `docs/library/saito-matsumoto-2011-tinymt.yaml`.
