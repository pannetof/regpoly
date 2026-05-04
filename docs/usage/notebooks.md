# Notebooks

The doc site renders Jupyter notebooks via [`mkdocs-jupyter`](https://github.com/danielfrg/mkdocs-jupyter). Notebooks live under `docs/notebooks/` and fall in three categories:

- `docs/notebooks/families/` — one notebook per generator family. Phase 7 authors all 18 from a common skeleton (`_template.ipynb`); each notebook ports a familiar parameter set from the catalog, verifies the characteristic polynomial and full-period claim, computes equidistribution `d(1..L)`, cross-checks against published values, and runs a small search for new parameters of the same shape.
- `docs/notebooks/research/` — preserved, including:
    - `melg_paper_verification.ipynb` — reproduces the MELG-19937 reference numbers.
    - `melg_charpoly_derivation.ipynb` — derives the MELG characteristic polynomial step by step.
- `docs/notebooks/utilities/` — preserved, including:
    - `bitmatrix_display.ipynb` — visualisation helpers for GF(2) matrices.
    - `primitive.ipynb` — primitivity-test utilities.

## Running locally

Open them with Jupyter:

```bash
uv sync --group docs
uv run --group docs jupyter lab
```

Or serve the built doc site (notebooks render to HTML):

```bash
uv run --group docs mkdocs serve
```

## Notebook tests

The `nbmake` plugin re-executes each notebook against tiny CI-budget overrides so that "API changed and broke the demo" regressions surface in the slow lane:

```bash
uv run pytest -m slow --nbmake docs/notebooks/families/
```

The slow lane is gated out of the default `pytest` invocation (`addopts = "-m 'not slow and not e2e'"`).

## See also

- [Python usage](python.md) — programmatic access to the same APIs the notebooks use.
- [Generator catalogue](../generators/index.md) — per-family theory pages each notebook links back to.
