# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.
"""Sphinx configuration for the REGPOLY documentation site.

Stack: Sphinx + Doxygen (via Breathe + Exhale) + autodoc + Napoleon
       + MyST-Parser + MyST-NB + jupyter-cache + sphinx-design
       + PyData Sphinx Theme.

Build commands:

    uv sync --group docs
    uv run --group docs sphinx-build -b html docs site            # standard
    uv run --group docs sphinx-build -W -b html docs site         # strict (CI)
    uv run --group docs sphinx-build -b doctest docs site_doctest # check examples
    uv run --group docs sphinx-build -b linkcheck docs site_lc    # check links

Doxygen + graphviz must be installed on the build host (apt: doxygen graphviz).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Make every workspace src/ tree importable so autodoc can find Python modules.
ROOT = Path(__file__).resolve().parent.parent
for pkg in ("regpoly", "regpoly-cpp", "regpoly-legacy", "regpoly-web"):
    src = ROOT / "packages" / pkg / "src"
    if src.exists():
        sys.path.insert(0, str(src))

# -- Project information ------------------------------------------------------

project = "REGPOLY"
copyright = "2026, Francois Panneton"
author = "Francois Panneton"
release = "2.0.0"

# -- General configuration ----------------------------------------------------

extensions = [
    # Core docstring + cross-ref machinery
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    # Markdown + notebooks — myst_nb loads myst_parser internally; do
    # NOT list both, or Sphinx errors on a duplicate-transform removal.
    "myst_nb",
    # C++ via Doxygen
    "breathe",
    "exhale",
    # UI extensions
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    # Redirects
    "sphinx_reredirects",
]

# myst_nb auto-registers .md (via myst_parser) and .ipynb. Don't override.
source_suffix = {
    ".rst":   "restructuredtext",
    ".md":    "myst-nb",
    ".ipynb": "myst-nb",
}

master_doc = "index"
exclude_patterns = [
    # Build-time generated catalog fragment, pulled into papers/index.md via
    # `{include}` — don't parse it as a standalone document.
    "papers/_auto_catalog.md",
    # Notebooks we DO NOT render
    "notebooks/research/**",
    # Markdown sources of notebooks the user doesn't want included
    "notebooks/families/_runner.py",
    "notebooks/families/_stamp.py",
    # Build outputs
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

primary_domain = "py"
# Default fenced code blocks to plain text — many narrative pages have
# pseudocode / ASCII art with Unicode arrows (→, —, ·) that fail
# Pygments' python3 lexer. Labeled fences (```python, ```cpp, etc.)
# still highlight correctly.
highlight_language = "text"

# -- HTML output --------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = "REGPOLY"
html_static_path = ["_static"]
html_css_files: list[str] = []
html_baseurl = "https://pannetof.github.io/regpoly/"

html_theme_options = {
    "show_nav_level": 2,
    "navigation_depth": 4,
    "show_toc_level": 2,
    "header_links_before_dropdown": 7,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/pannetof/regpoly",
            "icon": "fa-brands fa-github",
        },
    ],
    "use_edit_page_button": False,
    "navbar_align": "left",
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version", "theme-version"],
}

# Notebook pages get a "Download notebook" button in the sidebar so
# readers can grab the underlying .ipynb and run it locally. The
# custom `sourcelink` template (docs/_templates/sourcelink.html) adds
# the HTML5 `download` attribute so the browser actually downloads
# the file instead of rendering it as plain text.
html_sidebars = {
    "notebooks/**": ["sidebar-nav-bs", "page-toc", "sourcelink"],
}

# Drop the default `.txt` suffix Sphinx appends to source copies, so
# the linked file ends in `.ipynb` (preserves the right extension when
# the browser saves it).
html_sourcelink_suffix = ""

templates_path = ["_templates"]

# -- MyST options -------------------------------------------------------------

myst_enable_extensions = [
    "amsmath",       # \begin{align*} ... \end{align*}
    "colon_fence",   # :::{note} ...
    "deflist",       # def lists
    "dollarmath",    # $...$ and $$...$$
    "smartquotes",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3

# -- MyST-NB (notebooks) ------------------------------------------------------

# Hash-based notebook execution: first build on a fresh checkout
# executes every included notebook and stores the outputs in
# `.jupyter_cache/` (gitignored). Subsequent builds reuse the cache
# unless the notebook source changes.
#
# CI note: the docs workflow (Step 11) must wire `actions/cache@v4`
# against `.jupyter_cache/` keyed on a hash of:
#   - every .ipynb file under docs/notebooks/
#   - the compiled pybind11 extension (.venv/.../*_regpoly_cpp*.so)
# so that a C++ code change invalidates notebook outputs.
nb_execution_mode = "cache"

# Alias every kernel reference in a notebook (e.g. the "regpoly" custom
# kernel some pre-v2 notebooks were authored with) to the local venv's
# python3 kernel. Avoids any dependency on the user's global Jupyter
# kernel registry.
nb_kernel_rgx_aliases = {r".*": "python3"}

# Don't abort the docs build when one notebook errors during execution —
# the error cell renders in the output so the reader sees it. This is
# safer than -W on the slow lane because (a) pre-v2 reference notebooks
# may depend on external URLs or fixtures that have moved, and (b) we
# don't want a single broken notebook to break the entire docs deploy.
nb_execution_raise_on_error = False
nb_execution_timeout = 120          # per-cell seconds; was 600 (overkill)
nb_execution_allow_errors = True    # render error tracebacks instead of failing

# Notebooks excluded from execution (still rendered with their committed
# outputs). Add entries here whenever a notebook depends on network /
# moved fixtures / hardware not available in CI.
nb_execution_excludepatterns = [
    "**/research/**",                            # 2 CA research notebooks
    "**/reference/primitive.ipynb",              # fetches external URL
    # genf2w uses the pre-v2 `coeffs:` YAML key format; the new
    # factory expects `coeff` (singular, flat uint_vec) and the
    # YAML preprocessor's coeffs→coeff translation is a known
    # regression. The committed cell outputs render as-is so the
    # notebook's narrative is preserved until the loader is fixed.
    "**/reference/genf2w_equidistribution.ipynb",
]
nb_execution_cache_path = ".jupyter_cache"    # gitignored
nb_merge_streams = True                        # cleaner stdout/stderr in output cells

# -- Autodoc + Napoleon (NumPy style) -----------------------------------------

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_typehints_format = "short"
autodoc_class_signature = "separated"

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = True
napoleon_attr_annotations = True

# -- Intersphinx --------------------------------------------------------------

intersphinx_mapping = {
    "python":     ("https://docs.python.org/3", None),
    "numpy":      ("https://numpy.org/doc/stable/", None),
    "pandas":     ("https://pandas.pydata.org/docs/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "scipy":      ("https://docs.scipy.org/doc/scipy/", None),
}

# -- Breathe + Exhale (C++) ---------------------------------------------------

DOXY_OUT = str(ROOT / "docs" / "_doxygen" / "xml")
breathe_projects = {"regpoly_cpp": DOXY_OUT}
breathe_default_project = "regpoly_cpp"

exhale_args = {
    "containmentFolder":     "./api/cpp/library",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "C++ API",
    "doxygenStripFromPath":  str(ROOT / "packages" / "regpoly-cpp" / "src" / "include"),
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": "\n".join([
        f"INPUT                = {ROOT / 'packages' / 'regpoly-cpp' / 'src' / 'include'}",
        "FILE_PATTERNS        = *.h",
        "RECURSIVE            = YES",
        "GENERATE_HTML        = NO",
        "GENERATE_XML         = YES",
        "EXTRACT_ALL          = NO",
        "EXTRACT_PRIVATE      = NO",
        "EXTRACT_STATIC       = NO",
        "QUIET                = YES",
        # NOT YET enabled — promoted in Step 12 once backfill lands.
        "WARN_IF_UNDOCUMENTED = NO",
        # Don't fail the docs build on Doxygen warnings during the migration window.
        # Step 12 enables WARN_AS_ERROR = YES.
        "WARN_AS_ERROR        = NO",
        "PREDEFINED           = "
            'PYBIND11_NAMESPACE=pybind11 '
            'PYBIND11_MODULE(x,y)=void __pybind_module() '
            'REGPOLY_PUBLIC= '
            'REGPOLY_INTERNAL= ',
        # Skip internal-only headers by convention (underscore-prefixed basename).
        "EXCLUDE_PATTERNS     = */_*.h",
    ]),
}

# -- sphinx-reredirects (preserve old MkDocs URLs) ----------------------------
# Add concrete entries here when content moves; the redirect file lists them
# all explicitly so a docs reader's bookmark continues to work.
redirects: dict[str, str] = {
    # Examples (will fill in as URLs change):
    # "api/python/regpoly.library/": "api/python/regpoly.library.html",
}

# -- Suppress benign warnings -------------------------------------------------

# -- builder-inited hook: auto-stub per-paper pages from docs/library/*.yaml -

import yaml


def _format_authors(authors):
    if not authors:
        return ""
    names = [
        f"{a.get('given', '').strip()} {a.get('family', '').strip()}".strip()
        for a in authors
    ]
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return " and ".join(names)
    if len(names) <= 3:
        return ", ".join(names[:-1]) + ", and " + names[-1]
    return f"{names[0]} et al."


def _format_citation(meta):
    parts = []
    authors = _format_authors(meta.get("authors") or [])
    year = meta.get("year")
    title = meta.get("title")
    venue = meta.get("venue")
    volume = meta.get("volume")
    issue = meta.get("issue")
    pages = meta.get("pages")
    if authors:
        parts.append(f"{authors}, {year}." if year else f"{authors}.")
    if title:
        parts.append(f"*{title}.*")
    venue_bits = []
    if venue:
        venue_bits.append(str(venue))
    if volume:
        vol = f"**{volume}**"
        if issue:
            vol += f"({issue})"
        venue_bits.append(vol)
    if pages:
        venue_bits.append(f"pp. {pages}")
    if venue_bits:
        parts.append(", ".join(venue_bits) + ".")
    return " ".join(parts)


def _generate_paper_pages(app=None, config=None):
    """Sphinx `builder-inited` hook — auto-stub per-paper pages."""
    library = ROOT / "docs" / "library"
    papers_out = ROOT / "docs" / "papers"
    papers_out.mkdir(parents=True, exist_ok=True)

    paper_ids = []
    for yaml_path in sorted(library.glob("*.yaml")):
        if yaml_path.name.endswith("_params.yaml"):
            continue
        if yaml_path.name in {"index.md", "about.md"}:
            continue
        try:
            meta = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
        except yaml.YAMLError:
            continue
        if not isinstance(meta, dict) or "id" not in meta:
            continue

        pid = meta["id"]
        title = meta.get("title", pid)
        out = papers_out / f"{pid}.md"

        lines = [f"# {title}", ""]
        citation = _format_citation(meta)
        if citation:
            lines += [citation, ""]
        if meta.get("doi"):
            lines += [f"[DOI]({meta['doi']})", ""]
        if meta.get("bibkey"):
            lines += [f"BibTeX key: `{meta['bibkey']}`", ""]
        if meta.get("tags"):
            lines += ["Tags: " + " ".join(f"`{t}`" for t in meta["tags"]), ""]
        if meta.get("deferred"):
            lines += [
                "> **Deferred catalog entry.** This paper is referenced for "
                "citation but no instantiable generator entries are loaded — "
                "typically because MTToolBox does not ship a per-mexp "
                "`calc_equidist` binary for the family. The web UI shows the "
                "paper but does not surface runnable parameter sets.",
                "",
            ]
        if meta.get("abstract"):
            lines += ["## Abstract", "", str(meta["abstract"]).strip(), ""]

        generators = meta.get("generators") or []
        if generators:
            lines += [
                "## Generators in this paper",
                "",
                "| Library id | Display | Family | Combined | $L_{\\max}$ |",
                "|---|---|---|---|---|",
            ]
            for g in generators:
                gid = g.get("id", "")
                disp = g.get("display", gid)
                fam = g.get("family", "")
                combined = "yes" if g.get("combined") else "no"
                lmax = g.get("Lmax", "")
                lines.append(
                    f"| `{gid}` | {disp} | `{fam}` | {combined} | {lmax} |"
                )
            lines += [
                "",
                "Each entry is loaded by `regpoly.library.Catalog` and "
                "surfaced on the web UI under `/library/{paper_id}/{gen_id}`.",
                "",
            ]

        lines += [
            "---",
            "",
            f"_Source YAML:_ "
            f"[`docs/library/{pid}.yaml`](https://github.com/pannetof/"
            f"regpoly/blob/master/docs/library/{pid}.yaml). "
            f"This page is auto-generated at build time by the "
            f"`_generate_paper_pages` hook in `docs/conf.py`; edit the YAML "
            f"to change the page.",
            "",
        ]
        out.write_text("\n".join(lines), encoding="utf-8")
        paper_ids.append(pid)

    # Write the auto-generated "Catalog-backed paper pages" section to
    # `docs/papers/_auto_catalog.md`. The hand-authored `index.md` pulls
    # it in via a MyST `{include}` directive — that way the source
    # `index.md` stays committed-clean and the auto-file is gitignored.
    auto_path = ROOT / "docs" / "papers" / "_auto_catalog.md"
    if paper_ids:
        body = ["## Catalog-backed paper pages", "",
                "Auto-generated from `docs/library/*.yaml` at build time:", ""]
        for pid in sorted(paper_ids):
            body.append(f"- [`{pid}`]({pid}.md)")
        body.append("")
        auto_path.write_text("\n".join(body), encoding="utf-8")


def setup(app):
    app.connect("builder-inited", _generate_paper_pages)
    return {"parallel_read_safe": True, "parallel_write_safe": True}


# -- Suppress benign warnings -------------------------------------------------

suppress_warnings = [
    # Exhale registers some nested types both as members AND as standalone
    # pages → known-benign duplicate C++ declaration warning. Sphinx 9 uses
    # the "cpp" domain; "cpp.duplicate_declaration" works in newer versions.
    # See https://github.com/svenevs/exhale/issues/77.
    "cpp.duplicate_declaration",
    "ref.duplicate_declaration",
    # The Exhale "collapsible-lists" static asset includes a LICENSE.md that
    # Sphinx flags as not-in-toctree. Benign.
    "toc.not_included",
    # Notebooks rendered from .ipynb committed-outputs (nb_execution_mode="off").
    # Suppress the "no cache entry" notices.
    "mystnb.exec",
]
