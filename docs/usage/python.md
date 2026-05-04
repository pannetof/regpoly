# Using REGPOLY from Python

The `regpoly` package is a thin Python wrapper around the C++ core (`regpoly_cpp`). Use it when you want full search-loop ergonomics — DataFrames, plots, notebook integration — on top of the same algorithms the standalone CLI runs.

## Install

```bash
git clone https://github.com/pannetof/regpoly
cd regpoly
uv sync
```

`uv sync` installs the workspace in editable mode and builds the C++ extension via `scikit-build-core`. After it completes, the `regpoly` and `regpoly-web` console scripts are on the venv's `PATH`.

## Run a search from a YAML config

```bash
uv run regpoly shared/yaml/equidist/well19937a.yml
```

Test config templates live under [`shared/yaml/`](https://github.com/pannetof/regpoly/tree/master/shared/yaml). The search loop iterates a `Combination` of generator pools, runs the configured analysis (`equidistribution`, `collision_free`, `tuplets`), and writes a tested-generator YAML for any candidate that passes the acceptance predicate.

## Use the library programmatically

```python
from regpoly.core.combination import Combination
from regpoly.core.combination_build import build_combinaison_inputs
from regpoly.library import Catalog

cat = Catalog("docs/library")
cat.load()
_, gen = cat.generator("mt19937")

gen_lists, temperings = build_combinaison_inputs(gen.components, gen.Lmax)
comb = Combination.CreateFromFiles(gen_lists, gen.Lmax, temperings)
comb.reset()
```

`comb` is now a fully instantiated combined generator. From here:

- Run an equidistribution analysis on the first component:

    ```python
    from regpoly.analyses.pis import analyze_single_generator
    res = analyze_single_generator(comb[0])
    print(res["se"], res["gaps"])
    ```

- Iterate a Seek search:

    ```python
    from regpoly.search.seek import Seek
    seek = Seek.from_yaml("shared/yaml/equidist/well19937a.yml")
    seek.run()
    ```

## Layer boundaries

The `regpoly` package may import `regpoly_cpp`, but the web app must go through `regpoly` only — never `_cpp` directly. `regpoly.introspection` re-exports the small set of catalog/parameter helpers the web layer needs. This is enforced by `import-linter` contracts declared in the workspace `pyproject.toml`.

## See also

- [C++ usage](cpp.md) — for users who do not want a Python runtime.
- [Web UI](web.md) — for browser-based search and result browsing.
- [Notebooks](notebooks.md) — for per-family demos and equidistribution exploration.
- [Architecture](../dev/architecture.md) — for a deeper look at the package layout.
