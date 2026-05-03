# Library

Published parameter sets indexed by paper. Each YAML in this directory describes either:

- a **paper-organised catalog** (e.g. `matsumoto-nishimura-1998.yaml`): one file per published paper, listing the generator instances introduced or tabulated in that paper. Loaded by `regpoly.library.Catalog` and surfaced via the web `/library` UI.
- a **family-organised cross-check catalog** (e.g. `sfmt_params.yaml`): one file per generator family, containing every parameter set we cross-check against MTToolBox `calc_equidist`-style binaries.

See [About the catalog](about.md) for the full schema specification.

> **Status**: Phase 0 imports the existing catalog files. Phase 6 wires `mkdocs-gen-files` to auto-stub one rendered page per paper from the YAML metadata.
