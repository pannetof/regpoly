# Library

Published parameter sets indexed by paper. Each YAML in this directory describes either:

- a **paper-organised catalog** (e.g. `matsumoto-nishimura-1998.yaml`): one file per published paper, listing the generator instances introduced or tabulated in that paper. Loaded by `regpoly.library.Catalog` and surfaced via the web `/library` UI.
- a **family-organised cross-check catalog** (e.g. `sfmt_params.yaml`): one file per generator family, containing every parameter set we cross-check against MTToolBox `calc_equidist`-style binaries.

See [About the catalog](about.md) for the full schema specification.

Each paper YAML is auto-rendered as a page under [Papers → Catalog-backed (auto-generated)](../papers/index.md) at site-build time, via `docs/gen_paper_pages.py` and the `mkdocs-gen-files` plugin. Editing a `docs/library/*.yaml` updates the corresponding paper page on the next build — there is no on-disk markdown to keep in sync.
