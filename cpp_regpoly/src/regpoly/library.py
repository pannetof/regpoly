"""Published-generators library — paper-centric catalog.

Each ``docs/library/*.yaml`` file describes one *paper* (the unit the
reader thinks in) containing the bibliographic metadata, an abstract,
and a list of *generators* the paper proposed.  A generator is a
parametrised instance of an existing C++ generator family (no new
family code is introduced).

URL shape:
- ``/library``                         list of papers
- ``/library/<paper_id>``              paper detail + generator list
- ``/library/<paper_id>/<gen_id>``     individual generator + actions
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

_log = logging.getLogger(__name__)

_ID_RE = re.compile(r"^[a-z0-9][a-z0-9\-]*$")


@dataclass
class Author:
    family: str          # surname
    given: str = ""      # given names / initials

    def display(self) -> str:
        return f"{self.family}, {self.given}" if self.given else self.family

    def short(self) -> str:
        """Last-name form used in the paper pill."""
        return self.family


@dataclass
class Generator:
    """One generator proposed in a paper."""

    id: str              # globally unique slug (``lfsr113``)
    display: str
    family: str          # C++ family name
    target: str
    combined: bool
    Lmax: int
    components: list[dict]
    notes_md: str = ""
    starred: bool = False    # per-generator flag; the paper's own
                             # `starred` is stored on Paper directly
    errors: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


@dataclass
class Paper:
    """One bibliographic entry with any number of generators."""

    id: str              # slug, matches filename stem
    authors: list[Author]
    year: int
    title: str
    venue: str
    volume: str = ""
    issue: str = ""
    pages: str = ""
    doi: str = ""        # full https://doi.org/… URL
    pdf: str = ""        # repo-relative "papers/foo.pdf", optional
    bibkey: str = ""     # optional BibTeX key
    abstract_md: str = ""
    notes_md: str = ""
    tags: list[str] = field(default_factory=list)
    starred: bool = False
    generators: list[Generator] = field(default_factory=list)
    source_path: Path = field(default_factory=lambda: Path("."))
    source_mtime: float = 0.0
    errors: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors and all(g.valid for g in self.generators)

    # ── Presentation helpers ──────────────────────────────────────────

    def author_list_short(self) -> str:
        """E.g. ``Matsumoto and Nishimura`` or ``L'Ecuyer``."""
        n = len(self.authors)
        if n == 0:
            return ""
        if n == 1:
            return self.authors[0].short()
        if n == 2:
            return (f"{self.authors[0].short()} and "
                    f"{self.authors[1].short()}")
        return f"{self.authors[0].short()} et al."

    def display(self) -> str:
        """Sidebar label: ``Matsumoto & Nishimura 1998``."""
        n = len(self.authors)
        if n == 0:
            who = self.id
        elif n == 1:
            who = self.authors[0].short()
        elif n == 2:
            who = (f"{self.authors[0].short()} & "
                   f"{self.authors[1].short()}")
        else:
            who = f"{self.authors[0].short()} et al."
        return f"{who} {self.year}"

    def acmtrans_citation(self) -> str:
        """Render the reference in ACM Transactions style.

        Example::
           L'Ecuyer, P. 1999. Tables of Maximally Equidistributed
           Combined LFSR Generators. Math. Comp. 68, 225 (1999), 261–269.
           https://doi.org/10.1090/S0025-5718-99-00996-5
        """
        parts: list[str] = []
        if self.authors:
            block = _format_authors_acm(self.authors).rstrip(".")
            parts.append(f"{block}.")
        parts.append(f"{self.year}.")
        if self.title:
            parts.append(f"{self.title.rstrip('.')}.")
        venue_bits: list[str] = []
        if self.venue:
            venue_bits.append(self.venue)
        if self.volume:
            if self.issue:
                venue_bits.append(f"{self.volume}, {self.issue}")
            else:
                venue_bits.append(str(self.volume))
        venue_str = " ".join(venue_bits).strip()
        tail = f"({self.year})"
        if self.pages:
            tail += f", {self.pages}"
        tail += "."
        if venue_str:
            parts.append(f"{venue_str} {tail}")
        if self.doi:
            parts.append(self.doi)
        return " ".join(p for p in parts if p)


def _format_authors_acm(authors: list[Author]) -> str:
    """Format ``Family, G.`` joined with ``and``."""
    def _one(a: Author) -> str:
        if a.given:
            initials = "".join(
                (tok[0] + ".") for tok in re.split(r"\s+", a.given.strip())
                if tok)
            return f"{a.family}, {initials}"
        return a.family
    formatted = [_one(a) for a in authors]
    if len(formatted) == 1:
        return formatted[0]
    if len(formatted) == 2:
        return f"{formatted[0]} and {formatted[1]}"
    return ", ".join(formatted[:-1]) + f", and {formatted[-1]}"


def config_hash(
    family: str, params: dict, tempering: list[dict],
) -> str:
    """Stable short hash of a single component configuration."""
    payload = json.dumps(
        {"family": family, "params": params, "tempering": tempering},
        sort_keys=True, default=str,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


# ═══════════════════════════════════════════════════════════════════════
# Catalog
# ═══════════════════════════════════════════════════════════════════════


class Catalog:
    """Registry of papers and the generators they contain."""

    def __init__(self, library_dir: Optional[Path]) -> None:
        self.library_dir = library_dir
        self._papers: dict[str, Paper] = {}
        self._generator_index: dict[str, tuple[str, int]] = {}
        # ^ generator id -> (paper_id, index within paper.generators)

    def load(self) -> None:
        self._papers = {}
        self._generator_index = {}
        if self.library_dir is None or not self.library_dir.is_dir():
            return
        for path in sorted(self.library_dir.glob("*.yaml")):
            try:
                paper = _parse_paper(path)
            except Exception as exc:
                _log.error(
                    "library: failed to parse %s: %s", path, exc)
                continue
            self._insert(paper)

    def reload_if_stale(self) -> None:
        if self.library_dir is None or not self.library_dir.is_dir():
            self._papers = {}
            self._generator_index = {}
            return
        on_disk = {p.stem: p for p in self.library_dir.glob("*.yaml")}
        for pid in list(self._papers):
            if pid not in on_disk:
                self._drop(pid)
        for stem, path in on_disk.items():
            mtime = path.stat().st_mtime
            existing = self._papers.get(stem)
            if existing is not None and existing.source_mtime == mtime:
                continue
            try:
                paper = _parse_paper(path)
            except Exception as exc:
                _log.error(
                    "library: reload %s failed: %s", path, exc)
                continue
            self._drop(stem)
            self._insert(paper)

    def _insert(self, paper: Paper) -> None:
        if paper.id in self._papers:
            _log.error(
                "library: duplicate paper id %r (keeping first)", paper.id)
            return
        self._papers[paper.id] = paper
        for i, g in enumerate(paper.generators):
            if g.id in self._generator_index:
                paper.errors.append(
                    f"generator id {g.id!r} already exists in paper "
                    f"{self._generator_index[g.id][0]!r}")
                continue
            self._generator_index[g.id] = (paper.id, i)

    def _drop(self, paper_id: str) -> None:
        paper = self._papers.pop(paper_id, None)
        if paper is None:
            return
        for g in paper.generators:
            if self._generator_index.get(g.id, ("", -1))[0] == paper_id:
                self._generator_index.pop(g.id, None)

    # ── Views ─────────────────────────────────────────────────────────

    def papers(
        self, *,
        starred: bool | None = None,
        tag: str | None = None,
        include_invalid: bool = False,
    ) -> list[Paper]:
        out = list(self._papers.values())
        if not include_invalid:
            out = [p for p in out if p.valid]
        if starred is not None:
            out = [p for p in out if p.starred == starred]
        if tag:
            out = [p for p in out if tag in p.tags]
        return sorted(out, key=lambda p: (p.year, p.id))

    def paper(self, paper_id: str) -> Paper | None:
        return self._papers.get(paper_id)

    def generator(
        self, gen_id: str,
    ) -> tuple[Paper, Generator] | None:
        loc = self._generator_index.get(gen_id)
        if loc is None:
            return None
        pid, idx = loc
        paper = self._papers.get(pid)
        if paper is None or idx >= len(paper.generators):
            return None
        return paper, paper.generators[idx]

    def all_generators(
        self, *, family: str | None = None,
    ) -> list[tuple[Paper, Generator]]:
        out = []
        for p in sorted(self._papers.values(),
                        key=lambda p: (p.year, p.id)):
            for g in p.generators:
                if family and g.family != family:
                    continue
                out.append((p, g))
        return out


# ═══════════════════════════════════════════════════════════════════════
# YAML parsing
# ═══════════════════════════════════════════════════════════════════════


_KNOWN_TARGETS = ("tested_generator", "primitive_generator")


def _parse_paper(path: Path) -> Paper:
    import yaml

    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, dict):
        raise ValueError("top-level YAML must be a mapping")

    errors: list[str] = []
    pid = str(raw.get("id", "")).strip()
    if not _ID_RE.match(pid):
        errors.append(f"id must match [a-z0-9][a-z0-9-]*, got {pid!r}")
    elif pid != path.stem:
        errors.append(
            f"id {pid!r} must equal filename stem {path.stem!r}")

    authors = _parse_authors(raw.get("authors"), errors)
    year = raw.get("year")
    if not isinstance(year, int):
        errors.append("year must be an integer")
        year = 0

    title = str(raw.get("title") or "").strip()
    venue = str(raw.get("venue") or "").strip()
    if not title:
        errors.append("title is required")
    if not venue:
        errors.append("venue is required")

    doi = str(raw.get("doi") or "").strip()
    if not doi:
        errors.append("doi is required")
    elif not (doi.startswith("https://doi.org/")
              or doi.startswith("https://www.jstor.org/")
              or doi.startswith("http://www.jstor.org/")):
        errors.append(
            "doi must be a full https URL "
            "(doi.org or jstor.org stable link)")

    pdf = str(raw.get("pdf") or "").strip()
    bibkey = str(raw.get("bibkey") or "").strip()

    volume = str(raw.get("volume") or "").strip()
    issue = str(raw.get("issue") or "").strip()
    pages = str(raw.get("pages") or "").strip()

    abstract_md = str(raw.get("abstract") or "")
    notes_md = str(raw.get("notes") or "")
    tags_raw = raw.get("tags") or []
    tags = [str(t) for t in tags_raw] if isinstance(tags_raw, list) else []
    starred = bool(raw.get("starred", False))

    generators = _parse_generators(raw.get("generators"), errors)

    paper = Paper(
        id=pid,
        authors=authors,
        year=year,
        title=title,
        venue=venue,
        volume=volume,
        issue=issue,
        pages=pages,
        doi=doi,
        pdf=pdf,
        bibkey=bibkey,
        abstract_md=abstract_md,
        notes_md=notes_md,
        tags=tags,
        starred=starred,
        generators=generators,
        source_path=path,
        source_mtime=path.stat().st_mtime,
        errors=errors,
    )

    # Cross-checks that only matter when the structure is intact.
    if not errors:
        _validate_pdf_exists(paper)
        for g in paper.generators:
            if g.valid:
                _dry_build_generator(g, Lmax_default=paper_default_Lmax(paper))

    return paper


def _parse_authors(raw, errors: list[str]) -> list[Author]:
    if raw is None:
        errors.append("authors list is required")
        return []
    if not isinstance(raw, list) or not raw:
        errors.append("authors must be a non-empty list")
        return []
    out: list[Author] = []
    for i, a in enumerate(raw):
        if isinstance(a, str):
            # Accept plain-string authors: "Pierre L'Ecuyer" → split at
            # last whitespace for surname (simple heuristic, override
            # with the structured form when it matters).
            tokens = a.strip().rsplit(None, 1)
            if len(tokens) == 1:
                out.append(Author(family=tokens[0]))
            else:
                out.append(Author(given=tokens[0], family=tokens[1]))
            continue
        if not isinstance(a, dict):
            errors.append(f"authors[{i}] must be a string or mapping")
            continue
        family = str(a.get("family") or "").strip()
        if not family:
            errors.append(f"authors[{i}].family is required")
            continue
        out.append(Author(
            family=family,
            given=str(a.get("given") or "").strip(),
        ))
    return out


def _parse_generators(raw, errors: list[str]) -> list[Generator]:
    if raw is None or not isinstance(raw, list) or not raw:
        errors.append("generators must be a non-empty list")
        return []
    out: list[Generator] = []
    for i, g in enumerate(raw):
        if not isinstance(g, dict):
            errors.append(f"generators[{i}] must be a mapping")
            continue
        gerrors: list[str] = []
        gid = str(g.get("id", "")).strip()
        if not _ID_RE.match(gid):
            gerrors.append(
                f"id must match [a-z0-9][a-z0-9-]*, got {gid!r}")
        display = str(g.get("display") or gid).strip() or gid
        family = str(g.get("family") or "").strip()
        target = str(g.get("target") or "tested_generator").strip()
        combined = bool(g.get("combined", False))
        Lmax = int(g.get("Lmax") or 0)
        if target not in _KNOWN_TARGETS:
            gerrors.append(
                f"target must be one of {_KNOWN_TARGETS}, got {target!r}")
        if Lmax <= 0:
            gerrors.append(
                f"Lmax must be a positive int, got {Lmax!r}")
        if not family:
            gerrors.append("family is required")

        components = _normalize_components(g.get("components"), gerrors)
        if target == "primitive_generator":
            if combined:
                gerrors.append(
                    "primitive_generator target cannot be combined")
            if len(components) != 1:
                gerrors.append(
                    "primitive_generator target must have exactly "
                    "one component")
            elif components[0]["tempering"]:
                gerrors.append(
                    "primitive_generator target must have empty tempering")

        out.append(Generator(
            id=gid,
            display=display,
            family=family,
            target=target,
            combined=combined,
            Lmax=Lmax,
            components=components,
            notes_md=str(g.get("notes") or ""),
            starred=bool(g.get("starred", False)),
            errors=gerrors,
        ))
    return out


def _normalize_components(raw, errors: list[str]) -> list[dict]:
    if not isinstance(raw, list) or not raw:
        errors.append("components must be a non-empty list")
        return []
    out: list[dict] = []
    for idx, comp in enumerate(raw):
        if not isinstance(comp, dict):
            errors.append(f"components[{idx}] must be a mapping")
            continue
        family = str(comp.get("family", "")).strip()
        L = int(comp.get("L", 0))
        params = comp.get("params") or {}
        tempering = comp.get("tempering") or []
        if not family:
            errors.append(f"components[{idx}].family is required")
        if L <= 0:
            errors.append(f"components[{idx}].L must be a positive int")
        if not isinstance(params, dict):
            errors.append(f"components[{idx}].params must be a mapping")
            params = {}
        if not isinstance(tempering, list):
            errors.append(
                f"components[{idx}].tempering must be a list")
            tempering = []
        norm_temper: list[dict] = []
        for tj, t in enumerate(tempering):
            if not isinstance(t, dict) or "type" not in t:
                errors.append(
                    f"components[{idx}].tempering[{tj}] needs a "
                    "'type' key")
                continue
            norm_temper.append(
                {"type": str(t["type"]),
                 **{k: v for k, v in t.items() if k != "type"}})
        out.append({
            "family": family,
            "L": L,
            "params": dict(params),
            "tempering": norm_temper,
        })
    return out


def paper_default_Lmax(paper: Paper) -> int:
    """Pick a default Lmax for dry-builds when a generator forgot it."""
    for g in paper.generators:
        if g.Lmax:
            return g.Lmax
    return 32


def _dry_build_generator(g: Generator, Lmax_default: int) -> None:
    from regpoly.generateur import Generateur
    from regpoly.transformation import Transformation

    for idx, comp in enumerate(g.components):
        try:
            gen = Generateur.create(
                comp["family"], comp["L"], **comp["params"])
        except Exception as exc:
            g.errors.append(
                f"components[{idx}]: Generateur.create({comp['family']}, "
                f"L={comp['L']}, …) failed: {exc}")
            continue
        w_default = gen.params.get("w")
        for tj, step in enumerate(comp["tempering"]):
            t = dict(step)
            type_name = t.pop("type")
            if "w" not in t and w_default is not None:
                t["w"] = w_default
            try:
                Transformation.create(type_name, **t)
            except Exception as exc:
                g.errors.append(
                    f"components[{idx}].tempering[{tj}] "
                    f"({type_name}): {exc}")


def _validate_pdf_exists(paper: Paper) -> None:
    if not paper.pdf:
        return
    here = paper.source_path.resolve().parent
    for _ in range(5):
        if (here / "docs").is_dir() and here.name != "docs":
            pdf_path = here / paper.pdf
            if not pdf_path.is_file():
                paper.errors.append(
                    f"pdf points to missing file: {paper.pdf}")
            return
        here = here.parent


__all__ = [
    "Catalog",
    "Paper",
    "Generator",
    "Author",
    "config_hash",
]
