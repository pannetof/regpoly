# Tests

## Default lane

```bash
pytest                                       # all default-lane tests
pytest tests/test_mttoolbox_crosscheck.py    # just the cross-check
```

The default lane skips heavy cross-checks (entries with mexp > 20000),
configured via `pyproject.toml`'s `addopts = "-m 'not slow'"`.

## Slow lane

```bash
pytest -m slow                               # only slow tests
pytest -m "slow or not slow"                 # everything
```

Slow tests run MTToolBox / regpoly equidistribution on large mexps
(SFMT-44497, SFMT-86243, etc.); each may take many minutes per case.
Reference data for slow entries is captured by:

```bash
python tests/_gen_mttoolbox_reference.py --slow sfmt
```

(omit `--slow` to refresh only the fast entries; both passes are
incremental — previously-captured slow entries are preserved across
fast-only re-runs).

## Refreshing MTToolBox reference data

The cross-check test compares regpoly's d(1..5) against values
captured from MTToolBox's `calc_equidist`-style binaries. To
regenerate the reference modules under `tests/*_mttoolbox_d5.py`:

```bash
# Build the binaries once (see _mttoolbox_build.md):
bash <recipes from _mttoolbox_build.md>

# Capture default-lane entries:
python tests/_gen_mttoolbox_reference.py

# Capture slow-lane entries (long-running):
python tests/_gen_mttoolbox_reference.py --slow

# Or per family:
python tests/_gen_mttoolbox_reference.py sfmt mt xorshift

# Force re-capture of entries already on disk:
python tests/_gen_mttoolbox_reference.py --slow --refresh sfmt
```

The capture script writes the reference module incrementally — after
every successful entry — so a long --slow pass that hits a per-binary
timeout still leaves earlier successes on disk.  Re-invoking with the
same arguments resumes from the last unfinished entry by default
(skip-already-captured).  To override per-binary timeouts, set
`MTTOOLBOX_REF_TIMEOUT_S` (default 1800).

## Adding a new generator family

1. Drop the upstream params header into a new
   `docs/library/<family>_params.yaml` (see existing files for schema).
2. If the family needs a new MTToolBox binary, document the build
   recipe in `tests/_mttoolbox_build.md`.
3. Run `python tests/_gen_mttoolbox_reference.py <family>` to
   materialise `tests/<family>_mttoolbox_d5.py`.
4. The cross-check test will pick the new family up automatically
   (it walks every `*_params.yaml`).
