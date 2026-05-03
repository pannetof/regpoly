# XorShift128Gen — Marsaglia 128-bit Xor-shift

`XorShift128Gen` is the 128-bit-state Xor-shift generator of Marsaglia
(2003): four 32-bit words rotated through a triple-shift recurrence,
producing 32-bit output. Two recurrence shapes are supported via the
`pattern` field, matching MTToolBox's `samples/XORSHIFT/xorshift-2.cpp`
(canonical Marsaglia 2003 form) and `xorshift-5.cpp` (search-result
form).

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `a`, `b`, `c` | int | Shift triple (the search axis). |
| `pattern` | int | 0 = xorshift-2 form, 1 = xorshift-5 form. |

State is initialized as `w = ~v; x = y = z = v` for seed `v`.

## Provenance

- **Reference**: Marsaglia, G., "Xorshift RNGs", Journal of Statistical
  Software 8(14), 2003.
- **Cross-check fixtures**: `docs/library/xorshift_params.yaml`.
- **Library entry** (planned): `docs/library/marsaglia-2003-xorshift.yaml`.
