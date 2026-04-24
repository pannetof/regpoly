# Primitive Search YAML Format

Search for generators with full period (primitive characteristic polynomial).

## Format

```yaml
search:
  family: <C++ class name or legacy alias>
  L: <output resolution in bits>
  output: <optional output file path, default: <input>_results.yaml>
  limit:
    max_tries: <stop after N attempts>
    max_seconds: <stop after N seconds>

structural_params:
  <param>: <value>
  ...

fixed_params:
  <param>: <value>       # fixed to this value
  <param>:               # omitted value → randomized
```

## Sections

### `search`

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `family` | string | yes | Generator family name |
| `L` | int | yes | Output resolution in bits |
| `output` | string | no | Path for results YAML file |
| `limit.max_tries` | int | no | Maximum number of attempts |
| `limit.max_seconds` | float | no | Maximum wall-clock time |

Both limits are optional. If both are given, whichever is reached first
stops the search. If neither is given, the search runs forever (Ctrl-C
to stop).

### `structural_params`

Parameters that define the period k. These must always be provided and
are never randomized.

Call `Generateur.parameters(family)` to see which parameters are
structural for a given family.

### `fixed_params`

Parameters that are varied during the search. Two cases:

- **Key with a value** (`m: 1`): fixed to that value for every try.
- **Key with no value** (`a:`): randomized on each try using the
  generator's built-in randomization hints.

Parameters not listed in `fixed_params` that are non-structural and
have randomization hints will also be randomized automatically.

## Output Format

The search writes a YAML file with all found generators:

```yaml
family: TGFSR
L: 32
structural_params:
  w: 32
  r: 3
summary:
  tries: 5000
  found: 12
  elapsed_seconds: 3.45
results:
  - m: 1
    a: 2571403067
    k: 96
  - m: 2
    a: 3522167137
    k: 96
```

## Running

```bash
regpoly fullperiodsearch.TGFSR.w32r3.yaml
```

Or from Python:

```python
from regpoly.search_primitive import PrimitiveSearch

search = PrimitiveSearch.from_yaml("fullperiodsearch.TGFSR.w32r3.yaml")
results = search.run()
```

## Generator Family Names

| Family | C++ class name | Legacy aliases |
|--------|---------------|----------------|
| Polynomial LCG | `PolyLCG` | `polylcg` |
| Tausworthe | `Tausworthe` | `taus`, `taus2` |
| Twisted GFSR | `TGFSR` | `tgfsr`, `TGFSRGen` |
| Mersenne Twister | `MersenneTwister` | `MT` |
| GF(2^w) Poly LCG | `GenF2wPolyLCG` | `genf2w` |
| GF(2^w) LFSR | `GenF2wLFSR` | `genf2w` (with type=lfsr) |
| Matsumoto | `MatsumotoGen` | `matsumoto` |
| Marsaglia Xor-shift | `MarsaXorshiftGen` | `marsaxorshift` |
| AC-1D | `AC1DGen` | `AC1D` |
| WELL RNG | `WELLRNG` | `carry`, `Carry2Gen` |

See the per-generator documentation in `docs/generators/` for parameter
details, mathematical formulas, and search examples.
