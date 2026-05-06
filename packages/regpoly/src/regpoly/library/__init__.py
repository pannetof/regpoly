# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Published-generators library — paper-centric catalog.

Phase 3.2: the catalog parser, dataclasses, and YAML loader live in
C++ (regpoly_cpp.Catalog / Paper / CatalogGenerator / Author /
config_hash). This module is a thin shim that re-exports those types
under the names existing callers use, accepts both string and
pathlib.Path arguments to Catalog(), and adapts the Path-typed
attributes the web app expects.

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

from pathlib import Path
from typing import Optional

from regpoly_cpp import _regpoly_cpp as _cpp

# Re-export the C++ types under the names the existing Python callers
# use. CatalogGenerator → Generator (the runtime Generator class lives
# under regpoly.core.generator, not in this namespace).
Author = _cpp.Author
Generator = _cpp.CatalogGenerator
Paper = _cpp.Paper
config_hash = _cpp.config_hash


class Catalog:
    """Registry of papers and the generators they contain.

    Thin wrapper around regpoly_cpp.Catalog: accepts pathlib.Path,
    None, or str for ``library_dir`` and exposes ``source_path`` as a
    pathlib.Path on returned Paper objects (the web app type-checks
    against Path).
    """

    def __init__(self, library_dir: Optional[Path | str]) -> None:
        if library_dir is None:
            self._impl = _cpp.Catalog("")
        else:
            self._impl = _cpp.Catalog(str(library_dir))
        self.library_dir = (
            Path(library_dir) if library_dir is not None else None)

    def load(self) -> None:
        self._impl.load()

    def reload_if_stale(self) -> None:
        self._impl.reload_if_stale()

    def papers(
        self, *,
        starred: bool | None = None,
        tag: str | None = None,
        include_invalid: bool = False,
    ) -> list[Paper]:
        return self._impl.papers(
            starred=starred, tag=tag, include_invalid=include_invalid)

    def paper(self, paper_id: str) -> Paper | None:
        return self._impl.paper(paper_id)

    def generator(
        self, gen_id: str,
    ) -> tuple[Paper, Generator] | None:
        return self._impl.generator(gen_id)

    def all_generators(
        self, *, family: str | None = None,
    ) -> list[tuple[Paper, Generator]]:
        return self._impl.all_generators(family=family)


__all__ = [
    "Catalog",
    "Paper",
    "Generator",
    "Author",
    "config_hash",
]
