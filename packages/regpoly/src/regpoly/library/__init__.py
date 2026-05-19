# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Published-generators library — paper-centric catalog.

The catalog parser, dataclasses, and YAML loader all live in C++
(`regpoly_cpp.Catalog` / `Paper` / `CatalogGenerator` / `Author` /
`config_hash`). This module is a thin Python shim that:

- accepts both `str` and `pathlib.Path` arguments to `Catalog(...)`,
- exposes `source_path` as a `pathlib.Path` on returned `Paper` objects
  (the web app type-checks against `Path`),
- re-exports the C++ types under the names existing callers use.

Each ``docs/library/*.yaml`` file describes one *paper* (the unit the
reader thinks in) containing the bibliographic metadata, an abstract,
and a list of *generators* the paper proposed. A generator here is a
parametrised instance of an existing C++ generator family — no new
family code is introduced.

URL shape served by the web layer:

| URL                                   | Page                                |
|---------------------------------------|-------------------------------------|
| `/library`                            | list of papers                      |
| `/library/<paper_id>`                 | paper detail + generator list       |
| `/library/<paper_id>/<gen_id>`        | individual generator + actions      |

Examples
--------
Build a catalog from the workspace docs tree, list papers, and
look up one generator:

>>> from regpoly.library import Catalog
>>> cat = Catalog("docs/library")
>>> cat.load()                                   # doctest: +SKIP
>>> for paper in cat.papers():                   # doctest: +SKIP
...     print(paper.id, paper.title)
>>> hit = cat.generator("mt19937")               # doctest: +SKIP
>>> if hit is not None:                          # doctest: +SKIP
...     paper, gen = hit
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from regpoly_cpp import _regpoly_cpp as _cpp

#: Bibliographic author record (re-exported from the C++ catalog).
Author = _cpp.Author

#: A single parametrised generator entry inside a paper. Aliases
#: `regpoly_cpp.CatalogGenerator`. Distinct from the runtime
#: `regpoly.core.generator.Generator` (which is the *family* class —
#: this is a *parameter set*).
Generator = _cpp.CatalogGenerator

#: Paper-level record: bibliographic metadata + abstract + the list of
#: generators the paper proposes. Re-exported from the C++ catalog.
Paper = _cpp.Paper

#: Compute a content hash of a parameter dict. Used by the web app to
#: detect when two `Generator` entries describe identical parameters.
#: Signature: ``config_hash(params: dict) -> str``.
config_hash = _cpp.config_hash


class Catalog:
    """Registry of papers and the generators they contain.

    Thin Python wrapper around `regpoly_cpp.Catalog`. The C++ side
    walks `library_dir/*.yaml`, parses each file as a paper, validates
    that every generator entry resolves to a known C++ family, and
    indexes everything by paper / generator id.

    Attributes
    ----------
    library_dir
        `Path` to the directory the catalog reads from
        (`None` when constructed with no path).

    See Also
    --------
    :cpp:class:`regpoly::library::Catalog` : the C++ implementation this wraps
        (header: ``packages/regpoly-cpp/src/include/library/catalog.h``).
    """

    def __init__(self, library_dir: Optional[Path | str]) -> None:
        """Construct a catalog rooted at ``library_dir``.

        The constructor does not load anything from disk — call
        :meth:`regpoly.library.Catalog.load` (or
        :meth:`regpoly.library.Catalog.reload_if_stale`)
        explicitly.

        Parameters
        ----------
        library_dir
            Path to the catalog directory (containing
            paper YAML files). Either a `str`, a `pathlib.Path`,
            or `None`. When `None`, the underlying C++ catalog is
            constructed with an empty path and will fail any
            subsequent `load()` — use this only for tests that
            replace `self._impl` afterwards.
        """
        if library_dir is None:
            self._impl = _cpp.Catalog("")
        else:
            self._impl = _cpp.Catalog(str(library_dir))
        self.library_dir = (
            Path(library_dir) if library_dir is not None else None)

    def load(self) -> None:
        """Read every paper YAML file under `library_dir`.

        Builds the in-memory index of papers and generators. Idempotent:
        a second call rebuilds the index from scratch.

        Raises
        ------
        FileNotFoundError
            If `library_dir` does not exist.
        ValueError
            If any paper YAML fails to parse, or if two
            papers / generators share the same id.
        """
        self._impl.load()

    def reload_if_stale(self) -> None:
        """Re-`load` only if at least one source YAML has changed on disk.

        Cheap to call from request handlers: returns immediately when
        the catalog is up to date.
        """
        self._impl.reload_if_stale()

    def papers(
        self, *,
        starred: bool | None = None,
        tag: str | None = None,
        include_invalid: bool = False,
    ) -> list[Paper]:
        """List papers matching the given filters.

        Parameters
        ----------
        starred
            If `True`, return only papers flagged
            ``starred: true`` in their YAML. If `False`, return
            only un-starred. If `None`, no filter.
        tag
            If set, return only papers whose ``tags:`` list
            contains this value.
        include_invalid
            If `True`, also return papers whose
            generator entries failed validation. Default `False`
            (invalid papers are dropped silently).

        Returns
        -------
        list of Paper
            List of matching `Paper` objects, in catalog order.
        """
        return self._impl.papers(
            starred=starred, tag=tag, include_invalid=include_invalid)

    def paper(self, paper_id: str) -> Paper | None:
        """Look up one paper by id.

        Parameters
        ----------
        paper_id
            The paper's id (the YAML file's stem, e.g.
            ``"matsumoto-nishimura-1998"``).

        Returns
        -------
        Paper or None
            The matching `Paper`, or `None` if no paper has this id.
        """
        return self._impl.paper(paper_id)

    def generator(
        self, gen_id: str,
    ) -> tuple[Paper, Generator] | None:
        """Look up one generator by id, plus its containing paper.

        Parameters
        ----------
        gen_id
            The generator's id (unique across the entire
            catalog, e.g. ``"mt19937"``).

        Returns
        -------
        tuple of (Paper, Generator) or None
            A `(paper, generator)` tuple, or `None` if no generator
            has this id.
        """
        return self._impl.generator(gen_id)

    def all_generators(
        self, *, family: str | None = None,
    ) -> list[tuple[Paper, Generator]]:
        """List every generator in the catalog, optionally filtered by family.

        Parameters
        ----------
        family
            If set, return only generators whose ``family``
            field matches this string (canonical `-Gen` class
            name; legacy aliases are *not* expanded). If `None`,
            return all generators.

        Returns
        -------
        list of (Paper, Generator)
            List of `(paper, generator)` pairs in catalog order.
        """
        return self._impl.all_generators(family=family)


__all__ = [
    "Catalog",
    "Paper",
    "Generator",
    "Author",
    "config_hash",
]
