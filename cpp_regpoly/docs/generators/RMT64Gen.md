# RMT64Gen — Reducible Mersenne Twister (64-bit)

`RMT64Gen` is a 64-bit Mersenne-Twister variant with a *reducible*
characteristic polynomial (composite period rather than prime). It is
shipped by MTToolBox as `samples/RMT/rmt64dc` and uses the
non-primitive equidistribution analysis path (`notprimitive_de`).

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `mexp` | int | Total state degree. |
| `pos` | int | Recurrence offset in 64-bit words. |
| `mata` | uint64 | Matrix-A entry of the recurrence. |
| `parity` | uint64 | Parity vector (period certification). |
| `maskb`, `maskc` | uint64 | Tempering masks. |

## Provenance

- **Reference**: distributed with MTToolBox (`samples/RMT/`). A
  canonical academic reference is being tracked.
- **Cross-check fixtures**: `docs/library/rmt_params.yaml` (per-mexp
  parameter sets generated with `rmt64dc`).
