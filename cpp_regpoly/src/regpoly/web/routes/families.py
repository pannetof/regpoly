"""API endpoints for generator families, transformations, and tests.

These are read-only introspection routes: the web frontend uses them to
build dynamic parameter forms.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, HTTPException

import regpoly._regpoly_cpp as _cpp
from regpoly.generateur import _FAMILY_ALIASES
from regpoly.web.config import find_docs_dir


router = APIRouter()


# All C++ generator class names exposed by the factory.
# (The factory also accepts legacy aliases — see _FAMILY_ALIASES.)
KNOWN_FAMILIES: list[str] = [
    "PolyLCG",
    "Tausworthe",
    "TGFSR",
    "MersenneTwister",
    "MatsumotoGen",
    "MarsaXorshiftGen",
    "AC1DGen",
    "WELLRNG",
    "GenF2wLFSR",
    "GenF2wPolyLCG",
    "MELG",
    "SFMT",
    "MTGP",
]

# Transformation types exposed by the factory.
_KNOWN_TRANSFORMATIONS: list[str] = [
    "tempMK",
    "tempMK2",
    "laggedTempering",
    "permut",
]

# Available test types (read by AbstractTest._from_params dispatchers).
_KNOWN_TESTS: list[str] = [
    "equidistribution",
    "collision_free",
    "tuplets",
]


def _aliases_for(family: str) -> list[str]:
    return [a for a, canon in _FAMILY_ALIASES.items() if canon == family]


@router.get("/families")
async def list_families() -> list[dict]:
    result = []
    for fam in KNOWN_FAMILIES:
        try:
            specs = _cpp.get_gen_param_specs(fam)
        except Exception as e:
            # Skip families whose C++ bindings are missing for this build.
            continue
        result.append({
            "name": fam,
            "aliases": _aliases_for(fam),
            "params": list(specs),
        })
    return result


@router.get("/families/{family}")
async def family_detail(family: str) -> dict:
    try:
        specs = _cpp.get_gen_param_specs(family)
    except Exception:
        raise HTTPException(404, f"Unknown family: {family}")
    try:
        enumerable = bool(_cpp.family_is_enumerable(family))
    except Exception:
        enumerable = False
    return {
        "name": family,
        "aliases": _aliases_for(family),
        "params": list(specs),
        "enumerable": enumerable,
    }


@router.get("/families/{family}/docs")
async def family_docs(family: str) -> dict:
    docs_dir = find_docs_dir()
    if docs_dir is None:
        return {"family": family, "docs_md": "", "html": ""}
    md_path = Path(docs_dir) / f"{family}.md"
    if not md_path.exists():
        return {"family": family, "docs_md": "", "html": ""}
    md = md_path.read_text()
    html = _markdown_to_html(md)
    return {"family": family, "docs_md": md, "html": html}


@router.get("/transformations")
async def list_transformations() -> list[dict]:
    result = []
    for t in _KNOWN_TRANSFORMATIONS:
        try:
            specs = _cpp.get_trans_param_specs(t)
        except Exception:
            continue
        result.append({"name": t, "params": list(specs)})
    return result


@router.get("/transformations/{name}")
async def transformation_detail(name: str) -> dict:
    try:
        specs = _cpp.get_trans_param_specs(name)
    except Exception:
        raise HTTPException(404, f"Unknown transformation: {name}")
    return {"name": name, "params": list(specs)}


@router.get("/tests")
async def list_tests() -> list[dict]:
    """Static metadata about available tests and their YAML parameters."""
    return [
        {
            "name": "equidistribution",
            "description": "Dimension-of-equidistribution test (ME)",
            "params": [
                {"name": "max_gap_sum", "type": "int", "required": False,
                 "description": "Maximum allowed total defect (mse)"},
                {"name": "method", "type": "str", "required": False,
                 "default": "matricial",
                 "choices": ["matricial", "lattice", "harase", "nothing"]},
                {"name": "delta", "type": "list",
                 "description": "Per-resolution gap bounds",
                 "required": False},
            ],
        },
        {
            "name": "collision_free",
            "description": "Collision-free test (CF)",
            "params": [
                {"name": "max_gap_sum", "type": "int", "required": False,
                 "description": "msecf threshold"},
            ],
        },
        {
            "name": "tuplets",
            "description": "Tuplets criterion test",
            "params": [
                {"name": "dimensions", "type": "list", "required": True,
                 "description": "Per-level dimension array (s)"},
                {"name": "threshold", "type": "float", "required": False,
                 "description": "Acceptance threshold (mDD)"},
                {"name": "test_type", "type": "str", "required": False,
                 "default": "max", "choices": ["max", "sum"]},
            ],
        },
    ]


# ── Markdown rendering ────────────────────────────────────────────────────

def _markdown_to_html(md: str) -> str:
    """Render a markdown string, preserving LaTeX math spans verbatim.

    Python-Markdown will happily mangle LaTeX: underscores become
    <em>, backslash-commands get stripped, and display math wrapped in
    `$$…$$` may be split across `<p>` boundaries.  We pull math spans
    out before conversion, replace them with opaque placeholder tokens,
    run markdown, then reinsert the originals so KaTeX auto-render can
    typeset them in the browser.
    """
    import re
    try:
        import markdown
    except ImportError:
        return f"<pre>{md}</pre>"

    spans: list[str] = []

    def _stash(match: "re.Match[str]") -> str:
        spans.append(match.group(0))
        return f"@@MATH{len(spans) - 1}@@"

    # Order matters: match the longest delimiter first.
    patterns = [
        r"\$\$[\s\S]+?\$\$",           # $$ ... $$ (display)
        r"\\\[[\s\S]+?\\\]",           # \[ ... \] (display)
        r"\\\([\s\S]+?\\\)",           # \( ... \) (inline)
        r"(?<!\\)\$(?!\$)[^\n$]+?(?<!\\)\$",  # $ ... $ (inline)
    ]
    combined = re.compile("|".join(patterns))
    stashed = combined.sub(_stash, md)

    html = markdown.markdown(
        stashed, extensions=["fenced_code", "tables", "toc"]
    )

    def _unstash(match: "re.Match[str]") -> str:
        idx = int(match.group(1))
        return spans[idx]

    return re.sub(r"@@MATH(\d+)@@", _unstash, html)
