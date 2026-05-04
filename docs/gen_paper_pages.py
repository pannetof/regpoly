"""Phase 6.4 — auto-stub one rendered page per paper from docs/library/.

Wired into mkdocs via the `mkdocs-gen-files` plugin (see mkdocs.yml).
For every paper-organised catalog YAML under `docs/library/` (the
ones that carry `id`, `authors`, `title`, …), this script writes a
markdown page at `papers/{id}.md` summarising the paper and its
catalog generator entries.

`*_params.yaml` files (per-family cross-check fixtures, e.g.
`sfmt_params.yaml`) are skipped — they have no `authors:` / `title:`
metadata and are not papers.

The generated pages live in mkdocs-gen-files' virtual filesystem;
they are NOT written under `docs/papers/` on disk. mkdocs-material
sees them at build time and threads them into the site, but git
sees nothing — no source-of-truth duplication, and editing the
catalog YAML is the only way to change them.
"""

from __future__ import annotations

from pathlib import Path

import mkdocs_gen_files
import yaml

LIBRARY_DIR = Path("docs/library")


def _format_authors(authors: list[dict]) -> str:
    """Render a YAML authors list as 'Given Family, Given Family' or
    'Given Family et al.' for >3 entries."""
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


def _format_citation(meta: dict) -> str:
    """Single-line citation: Authors, Year. *Title.* Venue volume(issue), pages."""
    parts: list[str] = []
    authors = _format_authors(meta.get("authors") or [])
    year = meta.get("year")
    title = meta.get("title")
    venue = meta.get("venue")
    volume = meta.get("volume")
    issue = meta.get("issue")
    pages = meta.get("pages")
    if authors:
        if year:
            parts.append(f"{authors}, {year}.")
        else:
            parts.append(f"{authors}.")
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


def _write_paper_page(meta: dict) -> None:
    pid = meta.get("id")
    if not pid:
        return
    title = meta.get("title", pid)
    out_path = f"papers/{pid}.md"

    lines: list[str] = []
    lines.append(f"# {title}")
    lines.append("")

    citation = _format_citation(meta)
    if citation:
        lines.append(citation)
        lines.append("")

    if meta.get("doi"):
        lines.append(f"[DOI]({meta['doi']})")
        lines.append("")

    if meta.get("bibkey"):
        lines.append(f"BibTeX key: `{meta['bibkey']}`")
        lines.append("")

    if meta.get("tags"):
        tag_strs = [f"`{t}`" for t in meta["tags"]]
        lines.append("Tags: " + " ".join(tag_strs))
        lines.append("")

    if meta.get("deferred"):
        lines.append(
            "> **Deferred catalog entry.** This paper is referenced for "
            "citation but no instantiable generator entries are loaded — "
            "typically because MTToolBox does not ship a per-mexp "
            "`calc_equidist` binary for the family. The web UI shows the "
            "paper but does not surface runnable parameter sets."
        )
        lines.append("")

    if meta.get("abstract"):
        lines.append("## Abstract")
        lines.append("")
        lines.append(str(meta["abstract"]).strip())
        lines.append("")

    generators = meta.get("generators") or []
    if generators:
        lines.append("## Generators in this paper")
        lines.append("")
        lines.append("| Library id | Display | Family | Combined | $L_{\\max}$ |")
        lines.append("|---|---|---|---|---|")
        for g in generators:
            gid = g.get("id", "")
            disp = g.get("display", gid)
            fam = g.get("family", "")
            combined = "yes" if g.get("combined") else "no"
            lmax = g.get("Lmax", "")
            lines.append(
                f"| `{gid}` | {disp} | `{fam}` | {combined} | {lmax} |"
            )
        lines.append("")
        lines.append(
            "Each entry is loaded by `regpoly.library.Catalog` and "
            "surfaced on the web UI under `/library/{paper_id}/{gen_id}`."
        )
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append(
        f"_Source YAML: [`docs/library/{pid}.yaml`](https://github.com/"
        f"pannetof/regpoly/blob/master/docs/library/{pid}.yaml)._ "
        f"This page is auto-generated at build time by "
        f"`docs/gen_paper_pages.py` from the catalog YAML; edit the YAML "
        f"to change the page."
    )
    lines.append("")

    with mkdocs_gen_files.open(out_path, "w") as fp:
        fp.write("\n".join(lines))
    mkdocs_gen_files.set_edit_path(out_path, f"library/{pid}.yaml")


def _is_paper_yaml(path: Path) -> bool:
    """Skip *_params.yaml cross-check fixtures and the on-disk index."""
    if path.name.endswith("_params.yaml"):
        return False
    if path.name in {"index.md", "about.md"}:
        return False
    return True


def main() -> None:
    paper_ids: list[str] = []
    for yaml_path in sorted(LIBRARY_DIR.glob("*.yaml")):
        if not _is_paper_yaml(yaml_path):
            continue
        with yaml_path.open("r", encoding="utf-8") as fp:
            try:
                meta = yaml.safe_load(fp)
            except yaml.YAMLError:
                continue
        if not isinstance(meta, dict) or "id" not in meta:
            continue
        _write_paper_page(meta)
        paper_ids.append(meta["id"])

    # Augment the existing papers/index.md with an auto-generated
    # listing of every catalog-backed page. The on-disk index.md is
    # left intact; we append a Catalog section at the end.
    index_path = "papers/index.md"
    src_index = Path("docs/papers/index.md")
    if src_index.exists():
        with mkdocs_gen_files.open(index_path, "w") as fp:
            fp.write(src_index.read_text(encoding="utf-8"))
            if paper_ids:
                fp.write("\n\n## Catalog-backed paper pages\n\n")
                fp.write(
                    "Auto-generated from `docs/library/*.yaml` at build "
                    "time:\n\n"
                )
                for pid in sorted(paper_ids):
                    fp.write(f"- [`{pid}`]({pid}.md)\n")


main()
