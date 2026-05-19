# Notebooks

Worked Jupyter notebooks for REGPOLY. Rendered by Sphinx via
[MyST-NB](https://myst-nb.readthedocs.io/) with hash-based caching
(`nb_execution_mode = "cache"`) — only re-executed when their inputs
change. A first build on a fresh checkout takes several minutes; cached
builds are near-instant.

## Per-family demos

One notebook per generator family (or close to it). Each demo loads a
representative parameter set, verifies the characteristic polynomial,
computes the equidistribution profile, and cross-checks against the
published reference values. Stamped from a parametric template
(`docs/notebooks/families/_stamp.py`).

```{toctree}
:maxdepth: 1
:glob:

families/*Gen
families/CombinedGenerator
```

## Reference notebooks

Pre-v2 notebooks kept as reference. Searches, equidistribution
analyses, characteristic-polynomial derivations, primitivity tests.

```{toctree}
:maxdepth: 1
:glob:

reference/*
```
