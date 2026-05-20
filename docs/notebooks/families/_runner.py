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

    Resolves legacy family-name aliases (e.g. catalog stores
    ``MersenneTwister``; the requested family is ``MTGen``).
    """
    from regpoly.core.generator import resolve_family

    target = resolve_family(family)
    for paper in cat.papers():
        for g in getattr(paper, "generators", []) or []:
            comps = getattr(g, "components", None) or []
            if not comps:
                continue
            # Components surfaced via the C++ catalog are dicts with the
            # keys {"family", "L", "params", "tempering"}.
            comp = comps[0]
            comp_family = comp["family"] if isinstance(comp, dict) else getattr(comp, "family", None)
            if comp_family is None:
                continue
            if resolve_family(comp_family) == target:
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
            {
                "family": c["family"] if isinstance(c, dict) else getattr(c, "family", None),
                "L": (c.get("L") if isinstance(c, dict) else getattr(c, "L", None)),
            }
            for c in g.components
        ],
    }
    return comb, info


def analyze_first_component(comb) -> dict[str, Any]:
    """Run an equidistribution analysis on the first component
    (uses :func:`regpoly.analyses.pis.analyze_single_generator`).

    Returns a result dict on success; on a known algorithm-level
    limitation (PISCache requires a primitive characteristic
    polynomial — SIMD families like SFMT/dSFMT and the
    lagged-feedback families like MELG don't satisfy that today and
    need the notprimitive / simd_notprimitive methods, which the
    `analyze_single_generator` shim does not yet dispatch on
    automatically), returns a stub dict with the failure message so
    the notebook continues to render instead of crashing.
    """
    from regpoly.analyses.pis import analyze_single_generator

    t0 = time.time()
    try:
        res = analyze_single_generator(comb[0])
        res["walltime"] = time.time() - t0
        res["status"] = "ok"
    except RuntimeError as exc:
        res = {
            "status": "skipped",
            "reason": str(exc).split("\n")[0],
            "walltime": time.time() - t0,
            "se": None,
            "k": comb[0].k,
            "L": comb[0].L,
            "gaps": [],
            "hamming_weight": None,
        }
    return res


def summarise(family: str, lib_id: str | None, info: dict, ana: dict) -> str:
    """One-line per-key summary suitable for `print()` in a notebook."""
    base = [
        f"family            : {family}",
        f"library entry     : {lib_id or '(none — using fallback params)'}",
        f"k_g (combined)    : {info['k_g']}",
        f"J (components)    : {info['J']}",
        f"Lmax              : {info['Lmax']}",
    ]
    if ana.get("status") == "skipped":
        base.extend([
            f"k (component)     : {ana['k']}",
            f"L (component)     : {ana['L']}",
            f"analysis          : SKIPPED — {ana['reason']}",
            f"                    (this family needs the notprimitive /",
            f"                    simd_notprimitive equidistribution",
            f"                    method, which the auto-stub does not",
            f"                    dispatch on automatically yet.)",
        ])
    else:
        base.extend([
            f"se (Σ gaps)       : {ana['se']}",
            f"k (component)     : {ana['k']}",
            f"L (component)     : {ana['L']}",
            f"hamming weight    : {ana['hamming_weight']}",
            f"PIS walltime (s)  : {ana['walltime']:.4f}",
            f"first 8 gaps      : {ana['gaps'][:8]}",
        ])
    return "\n".join(base)


def run_for_family(family: str, library_id: str | None = None) -> dict:
    """End-to-end driver. Loads catalog, picks an example for `family`
    (or the named library entry if `library_id` is set), instantiates,
    analyses, returns the assembled result dict.

    If neither the named library entry nor a catalog walk finds an
    example for the family, returns a stub result with status
    ``"no_example"`` so the notebook continues to render.
    """
    cat = load_catalog()
    g = None
    if library_id:
        loc = cat.generator(library_id)
        if loc is not None:
            _, g = loc
    if g is None:
        match = find_example_for_family(cat, family)
        if match is not None:
            library_id, g = match

    if g is None:
        info = {"Lmax": None, "k_g": None, "J": None, "components": []}
        ana = {
            "status": "no_example",
            "reason": (
                f"No catalog entry exercises {family} as its first "
                f"component (yet)."
            ),
            "walltime": 0.0,
            "se": None, "k": None, "L": None,
            "gaps": [], "hamming_weight": None,
        }
        return {
            "family": family,
            "library_id": None,
            "info": info,
            "analysis": ana,
            "summary": (
                f"family            : {family}\n"
                f"library entry     : (none — {ana['reason']})\n"
                f"analysis          : SKIPPED"
            ),
        }

    comb, info = instantiate(g)
    ana = analyze_first_component(comb)
    return {
        "family": family,
        "library_id": library_id,
        "info": info,
        "analysis": ana,
        "summary": summarise(family, library_id, info, ana),
    }
