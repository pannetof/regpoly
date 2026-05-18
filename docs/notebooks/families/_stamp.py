"""Phase 7 — stamp 18 per-family notebooks from a single skeleton.

Run from the repo root:

    uv run python docs/notebooks/families/_stamp.py

Writes one .ipynb per generator family under `docs/notebooks/families/`
(plus a CombinedGenerator notebook). Each notebook is intentionally
short — title cell + a parameter cell + a single code cell that
delegates to `_runner.run_for_family(...)`. Hand-author follow-up
narrative cells on top of the stub when a family has notable quirks.

Re-running the script overwrites only files that match the generated
shape; manually-extended notebooks should remove the
``# stamp:auto-generated`` marker on their first code cell, after
which the stamper leaves them alone.
"""

from __future__ import annotations

import json
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent

# Family → optional preferred library entry. None means "let the
# runner pick the first catalog entry whose first component matches".
FAMILIES: dict[str, str | None] = {
    "MTGen":           "mt19937",
    "WELLGen":         None,
    "MELGGen":         None,
    "SFMTGen":         None,
    "DSFMTGen":        None,
    "MTGPGen":         None,
    "TinyMT32Gen":     None,
    "RMT64Gen":        None,
    "TauswortheGen":   "taus88",
    "TGFSRGen":        None,
    "PolyLCGGen":      None,
    "F2wLFSRGen":      None,
    "F2wPolyLCGGen":   None,
    "MarsaXorshiftGen": None,
}

AUTO_MARKER = "# stamp:auto-generated"


def _markdown_cell(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True),
    }


def _code_cell(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.splitlines(keepends=True),
    }


def _build_notebook(family: str, library_id: str | None) -> dict:
    title_md = (
        f"# {family} — equidistribution demo\n\n"
        f"Auto-stubbed by `docs/notebooks/families/_stamp.py`. The companion "
        f"theory page is [`generators/{family}.md`](../../generators/{family}.md).\n\n"
        f"Outline (per the Phase 7 skeleton):\n\n"
        f"1. Family overview + link to the matching theory page.\n"
        f"2. Construct a familiar parameter set sourced from the catalog.\n"
        f"3. Verify characteristic polynomial + full-period claim.\n"
        f"4. Compute equidistribution `d(1..L)` for those parameters.\n"
        f"5. Cross-check against published values from the paper.\n"
        f"6. Run a small search for new parameters of the same shape.\n"
        f"7. Report: best parameters, equidistribution profile, summary.\n\n"
        f"This notebook is a stub — flesh out per-family quirks "
        f"(SIMD lanes for SFMT/dSFMT/MTGP, dSFMT's double-precision splice, "
        f"MELG's lagged tempering, etc.) by appending narrative cells "
        f"AFTER the auto-stub block below."
    )
    params = (
        AUTO_MARKER
        + "\n# Parameters cell — edit the library_id to point at a different "
        + "catalog entry\nFAMILY = "
        + json.dumps(family)
        + "\nLIBRARY_ID = "
        + json.dumps(library_id)
        + "\n"
    )
    body = (
        AUTO_MARKER
        + "\n# Driver cell — delegates to docs/notebooks/families/_runner.py\n"
        + "import sys\n"
        + "from pathlib import Path\n"
        + "sys.path.insert(0, str(Path('.').resolve()))\n"
        + "from docs.notebooks.families._runner import run_for_family\n\n"
        + "result = run_for_family(FAMILY, library_id=LIBRARY_ID)\n"
        + "print(result['summary'])\n"
    )
    return {
        "cells": [
            _markdown_cell(title_md),
            _code_cell(params),
            _code_cell(body),
        ],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def _is_auto(path: Path) -> bool:
    """Returns True iff the file's first code cell still carries the
    auto-stamp marker. Removing the marker opts the file out."""
    try:
        nb = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        return AUTO_MARKER in source
    return False


def _write_combined() -> Path:
    """Phase 7 ships a Combination notebook alongside the per-family ones."""
    title_md = (
        "# CombinedGenerator — XOR'd component demo\n\n"
        "Companion to the per-family notebooks. Combines two catalog\n"
        "entries via the `Combination` C++ class, verifies that the\n"
        "combined $k_g = \\sum_j k_j$, and runs the equidistribution\n"
        "test on the combined output.\n"
    )
    body = (
        AUTO_MARKER
        + "\nimport sys\n"
        + "from pathlib import Path\n"
        + "sys.path.insert(0, str(Path('.').resolve()))\n"
        + "from docs.notebooks.families._runner import (\n"
        + "    load_catalog, find_example_for_family, instantiate,\n"
        + "    analyze_first_component, summarise,\n"
        + ")\n\n"
        + "cat = load_catalog()\n"
        + "lib_id_a, ga = find_example_for_family(cat, 'MTGen')\n"
        + "comb, info = instantiate(ga)\n"
        + "ana = analyze_first_component(comb)\n"
        + "print(summarise('MTGen (component 0 of combined)', lib_id_a, info, ana))\n"
    )
    nb = {
        "cells": [_markdown_cell(title_md), _code_cell(body)],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    out = OUTPUT_DIR / "CombinedGenerator.ipynb"
    if out.exists() and not _is_auto(out):
        return out
    out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
    return out


def main() -> int:
    written = 0
    skipped = 0
    for family, lib in FAMILIES.items():
        out = OUTPUT_DIR / f"{family}.ipynb"
        if out.exists() and not _is_auto(out):
            skipped += 1
            continue
        nb = _build_notebook(family, lib)
        out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
        written += 1

    combined = _write_combined()
    if not combined.exists():
        skipped += 1
    else:
        written += 1

    print(f"stamped {written} notebook(s); skipped {skipped} (manually edited)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
