# Using REGPOLY from Python

> **Status**: Phase 0 placeholder. Final content authored in Phase 6 against the v2.0 Python surface.

Quick orientation while migration is underway:

```bash
uv sync
uv run regpoly shared/yaml/equidist/well19937a.yml
```

The Python `regpoly` package is a thin wrapper around the C++ core. After Phase 2 most calls will be one-line delegations into `regpoly_cpp`; today there is some Python orchestration logic that is still in the process of being moved.
