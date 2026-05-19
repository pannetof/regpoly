# {Family display name}

<!--
  Per-family page template.
  Copy this file to <Family>Gen.md, fill in each section against the
  C++ class in packages/regpoly-cpp/src/generators/, and add it to the
  toctree in docs/generators/index.md (Sphinx + MyST picks it up
  automatically once the file exists and is listed).
  Already-authored examples:
    - generators/MTGen.md           (full, canonical)
    - generators/WELLGen.md         (full)
    - generators/MELGGen.md         (compact)
  Sections marked OPTIONAL can be omitted on simple families.
-->

**C++ class:** `{ClassName}`
**Name in code:** `"{string accepted by Generator.create()}"`
**Legacy aliases:** `{old}, {names}` *(omit if none)*

One-paragraph orientation: what the family is, the historical paper(s)
it traces back to, and the main reason a user would choose it
(period, equidistribution, speed, SIMD layout, etc).

## Mathematical recurrence

The state is `{describe shape: r words of w bits, or k-bit packed
vector, etc}`. At each step:

```
{the family's word-sized recurrence in pseudo-code,
 mirroring the description in the source comments}
```

If the family's $A$ matrix has a closed-form description (companion
matrix, twist matrix, polyomial recurrence), give it. Otherwise
point at the C++ source file.

### Characteristic polynomial

Either a closed-form expression (for TGFSR-like families) or a
short statement of how it is computed at runtime (Berlekamp–Massey
on the bit-stream, factor of a known polynomial, etc).

## Parameters

| Name | Type  | Role       | Rand type | Rand args | Description |
|------|-------|------------|-----------|-----------|-------------|
| `{p}`| `int` | structural | --        | --        | {what it controls} |
| `{q}`| `int` | search     | `range`   | `1, k-1`  | {what gets randomized over} |

Roles:

- **structural** — fixed by the user; defines the shape of $A$.
- **search** — randomized by the search loop; the values that get
  reported when a candidate passes equidistribution.

The `Rand type` and `Rand args` columns name the C++ sampler in
`regpoly_random::sample_param_into`.

## State size (period exponent)

```
k = {formula in terms of structural params}
```

Classical examples (one bullet per published parameter set):

- **{Name}:** {structural params} ⇒ $k = {value}$.

## Search examples

### {Short title}

```yaml
search:
  family: {ClassName}
  L: 32
  limit:
    max_tries: 100000

structural_params:
  {p}: {value}

fixed_params:
  {q}:    # randomized in {Rand args}
```

OPTIONAL: one or two more YAML snippets covering edge cases the family
exposes (smaller period, fixed search params, alternate word width).

## Tempering

OPTIONAL section: list the canonical tempering chain associated with
the family (e.g. MT's two shift-mask steps, MELG's lagged tempering,
SFMT's vector tempering). Link to
[`theory/tempering_optimization.md`](../theory/tempering_optimization.md)
for the search-side details.

## Implementation notes

OPTIONAL: any non-obvious thing about the implementation that a user
porting parameters between REGPOLY and an upstream reference
implementation needs to know — bit-order conventions, packed-word
layout, SIMD lane assumptions, dSFMT's double-precision splice, etc.

---

## References

- {First author}, {Year}. *{Paper title}.*
  {Journal} {volume} ({year}), {pages}.
  [DOI](https://doi.org/{...}).
  Published as [`{lib_id}`](../papers/{paper_id}.md) in the library.
- {Other paper if relevant}.
