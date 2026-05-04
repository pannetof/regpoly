"""Phase 7 — shared per-family notebook helpers.

Each family notebook is intentionally small: it sets a `FAMILY` and
optional `EXAMPLE_LIB_ID`, then calls into this module. Per-family
quirks (SIMD lanes for SFMT, dSFMT's double-precision splice,
MELG's lagged tempering) get their own narrative cells appended by
hand on top of the generated skeleton.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

# Anchors the catalog at the in-tree library/, so the notebooks work
# whether they run from the repo root or from docs/notebooks/families/.
DEFAULT_LIBRARY_DIR = Path(__file__).resolve().parents[2] / "library"


def load_catalog(library_dir: Path | None = None):
    """Open the in-tree catalog and return the loaded `Catalog`."""
    from regpoly.library import Catalog

    cat = Catalog(str(library_dir or DEFAULT_LIBRARY_DIR))
    cat.load()
    return cat


def find_example_for_family(cat, family: str) -> tuple[str, Any] | None:
    """Walk the catalog for the first generator whose first component
    is the requested family. Returns (library_id, generator) or None.
    """
    for paper in cat.papers():
        for g in getattr(paper, "generators", []) or []:
            comps = getattr(g, "components", None) or []
            if not comps:
                continue
            comp_family = getattr(comps[0], "family", None)
            if comp_family == family:
                return g.id, g
    return None


def instantiate(g) -> tuple[Any, Any]:
    """Build a Combination from a catalog entry. Returns (comb, info)."""
    from regpoly.core.combination import Combination
    from regpoly.core.combination_build import build_combinaison_inputs

    gen_lists, temperings = build_combinaison_inputs(g.components, g.Lmax)
    comb = Combination.CreateFromFiles(gen_lists, g.Lmax, temperings)
    comb.reset()
    info = {
        "Lmax": g.Lmax,
        "k_g": comb.k_g,
        "J": comb.J,
        "components": [
            {"family": c.family, "L": getattr(c, "L", None)}
            for c in g.components
        ],
    }
    return comb, info


def analyze_first_component(comb) -> dict[str, Any]:
    """Run an equidistribution analysis on the first component
    (uses regpoly.analyses.pis.analyze_single_generator)."""
    from regpoly.analyses.pis import analyze_single_generator

    t0 = time.time()
    res = analyze_single_generator(comb[0])
    res["walltime"] = time.time() - t0
    return res


def summarise(family: str, lib_id: str | None, info: dict, ana: dict) -> str:
    """One-line per-key summary suitable for `print()` in a notebook."""
    lines = [
        f"family            : {family}",
        f"library entry     : {lib_id or '(none — using fallback params)'}",
        f"k_g (combined)    : {info['k_g']}",
        f"J (components)    : {info['J']}",
        f"Lmax              : {info['Lmax']}",
        f"se (Σ gaps)       : {ana['se']}",
        f"k (component)     : {ana['k']}",
        f"L (component)     : {ana['L']}",
        f"hamming weight    : {ana['hamming_weight']}",
        f"PIS walltime (s)  : {ana['walltime']:.4f}",
        f"first 8 gaps      : {ana['gaps'][:8]}",
    ]
    return "\n".join(lines)


def run_for_family(family: str, library_id: str | None = None) -> dict:
    """End-to-end driver. Loads catalog, picks an example for `family`
    (or the named library entry if `library_id` is set), instantiates,
    analyses, returns the assembled result dict."""
    cat = load_catalog()
    if library_id:
        loc = cat.generator(library_id)
        if loc is None:
            raise RuntimeError(f"library_id not in catalog: {library_id}")
        _, g = loc
    else:
        match = find_example_for_family(cat, family)
        if match is None:
            raise RuntimeError(
                f"No catalog entry has {family} as its first component. "
                f"Provide library_id explicitly or supply fallback params."
            )
        library_id, g = match

    comb, info = instantiate(g)
    ana = analyze_first_component(comb)
    return {
        "family": family,
        "library_id": library_id,
        "info": info,
        "analysis": ana,
        "summary": summarise(family, library_id, info, ana),
    }
