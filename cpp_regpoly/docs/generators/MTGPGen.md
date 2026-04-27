# MTGPGen — Mersenne Twister for Graphic Processors

`MTGPGen` is the GPU-oriented MT variant of Saito and Matsumoto (2013).
The state is partitioned across small "blocks" so the recurrence can be
evaluated cooperatively by GPU threads, while preserving an MT-style
Mersenne-prime period.

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `mexp` | int | Mersenne exponent. |
| `pos1`, `sh1`, `sh2` | int | Recurrence offset and shift counts. |
| `mask` | uint | Element mask (per-block). |
| `tbl` | uint[16] | Tempering / output table. |
| `tmp` | uint[16] | Temporary tempering table. |
| `flt` | uint[16] | Floating-point output table. |

## Provenance

- **Reference**: Saito, M., Matsumoto, M., "Variants of Mersenne
  Twister suitable for graphic processors", ACM Transactions on
  Mathematical Software (TOMS), 2013.
- **Cross-check fixtures**: `docs/library/mtgp_params.yaml` (deferred —
  upstream MTToolBox does not ship a per-mexp `calc_equidist` for MTGP).
