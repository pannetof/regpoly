# Generator parameter catalogs

This directory holds two kinds of YAML files; the library loader
distinguishes them by filename (`*_params.yaml` is treated as a
cross-check fixture, never as a paper).

1. **Paper-organised catalogs** (`saito-matsumoto-2008.yaml`,
   `matsumoto-nishimura-1998.yaml`, …): one file per published paper,
   listing the generator instances introduced or tabulated in that
   paper. These are loaded by `regpoly.library.Catalog` and surfaced
   via the web `/library` UI for academic reference and reproductions.

   - Paper YAMLs may set `deferred: true` to indicate that no
     ready-to-run generator entries exist yet (e.g. when MTToolBox
     does not ship a per-MEXP `calc_equidist` binary for the
     family). Deferred papers are loaded with an empty `generators:`
     list; they appear in the catalog for citation but are not
     surfaced as "instantiable" entries on the web UI.

2. **Family-organised cross-check catalogs** (`*_params.yaml`,
   e.g. `sfmt_params.yaml`, `dsfmt_params.yaml`): one file per
   F₂-linear generator family, containing every parameter set we
   cross-check against MTToolBox `calc_equidist`-style binaries.
   These are the source of truth for
   `tests/_gen_mttoolbox_reference.py` (which materialises the
   per-family `tests/*_mttoolbox_reference.py` modules) and for
   `tests/test_mttoolbox_crosscheck.py`. The library loader skips
   files whose name ends in `_params.yaml`.

## Cross-check catalog schema (`*_params.yaml`)

```yaml
family: SFMT                # regpoly Generateur family name
mttoolbox_binary: sfmtdc/calc_equidist   # path under MTToolBox/samples/
slow_threshold_mexp: 20000  # entries with mexp > this use pytest -m slow
generators:
  - id: sfmt607             # arbitrary stable identifier
    mexp: 607
    # --- params consumed by regpoly's Generateur.create() ---
    regpoly:
      pos1: 2
      sl1: 15
      sl2: 3
      sr1: 13
      sr2: 3
      msk1: 0xfdff37ff
      msk2: 0xef7f3f7d
      msk3: 0xff777b7d
      msk4: 0x7ff7fb2f
    # --- argv format consumed by calc_equidist (one positional arg, comma-sep) ---
    mttoolbox_argv:
      - "607,2,15,3,13,3,0xfdff37ff,0xef7f3f7d,0xff777b7d,0x7ff7fb2f,0x00000001,0x00000000,0x00000000,0x5986f054"
```

### Optional fields

- `cross_check_xfail: <reason>` — translates into
  `pytest.mark.xfail(strict=False)`.  Use for entries where regpoly's
  d(v) differs from MTToolBox by a known (and accepted) margin.
- `tempering: [...]` — list of `{type, params}` transformations applied
  to the generator output before computing d(v).  Each entry is fed
  through `Transformation.create(type, **params)` and wired through
  `Combinaison`.  Used e.g. by MT19937 to mirror MTToolBox/mt1's
  built-in MT19937 tempering pipeline.

```yaml
    tempering:
      - type: tempMK2
        params:
          w: 32
          u: 11
          eta: 7
          mu: 15
          l: 18
          b: 0x9D2C5680
          c: 0xEFC60000
```

Source of truth for parameter values is upstream headers
(`SFMT-paramsNNN.h`, `dSFMT-paramsNNN.h`, etc.). When updating,
re-run `python tests/_gen_mttoolbox_reference.py` to refresh the
per-family `*_mttoolbox_reference.py` modules.

## Adding a new family

1. Drop the upstream parameter header into a comment block at the
   top of the new YAML, citing the source.
2. List one entry per `(mexp, params)` pair you want cross-checked.
3. Pick a `mttoolbox_binary` (build it via the recipe in
   `tests/_mttoolbox_build.md` if not already built).
4. Run `python tests/_gen_mttoolbox_reference.py` to populate the
   reference module.
5. Add an arm to `tests/test_mttoolbox_crosscheck.py` that loads
   the new YAML and dispatches to the right regpoly equidistribution
   method.
