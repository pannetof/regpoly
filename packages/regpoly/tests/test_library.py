"""Smoke tests for the paper-centric published-generators catalog."""

from __future__ import annotations

import os
import time
from pathlib import Path

from regpoly.library import Catalog, Paper, config_hash
from regpoly.web.config import find_library_dir


def _repo_catalog() -> Catalog:
    d = find_library_dir()
    assert d is not None, "docs/library dir not found"
    c = Catalog(d)
    c.load()
    return c


def test_catalog_loads_papers():
    c = _repo_catalog()
    ids = {p.id for p in c.papers(include_invalid=True)}
    assert {"matsumoto-nishimura-1998",
            "lecuyer-1996", "lecuyer-1999"} <= ids


def test_papers_are_valid():
    c = _repo_catalog()
    invalid = [p for p in c.papers(include_invalid=True) if not p.valid]
    if invalid:
        lines = []
        for p in invalid:
            lines.append(f"{p.id}: {p.errors}")
            for g in p.generators:
                if not g.valid:
                    lines.append(f"    {g.id}: {g.errors}")
        raise AssertionError("Invalid papers:\n" + "\n".join(lines))


def test_generator_lookup_is_unique():
    c = _repo_catalog()
    gen_ids = {"mt19937", "taus88", "lfsr113"}
    for gid in gen_ids:
        loc = c.generator(gid)
        assert loc is not None, f"{gid} missing"
        paper, gen = loc
        assert gen.id == gid


def test_paper_display_format():
    c = _repo_catalog()
    p = c.paper("matsumoto-nishimura-1998")
    assert p.display() == "Matsumoto & Nishimura 1998"
    p = c.paper("lecuyer-1999")
    assert p.display() == "L'Ecuyer 1999"


def test_acmtrans_citation_has_all_pieces():
    c = _repo_catalog()
    cite = c.paper("lecuyer-1999").acmtrans_citation()
    # Must contain author surname, year, title opening, venue, pages, DOI.
    assert "L'Ecuyer" in cite
    assert "1999" in cite
    assert "Tables of Maximally" in cite
    assert "Mathematics of Computation" in cite
    assert "261–269" in cite
    assert "https://www.jstor.org/stable/2585109" in cite


def test_config_hash_stable_under_key_reorder():
    params_a = {"k": 31, "nb_terms": 3, "quicktaus": True,
                "poly": [0, 6, 31], "s": 18}
    params_b = {"s": 18, "poly": [0, 6, 31], "nb_terms": 3,
                "k": 31, "quicktaus": True}
    assert config_hash("Tausworthe", params_a, []) == \
           config_hash("Tausworthe", params_b, [])


def test_reload_if_stale(tmp_path: Path):
    src = find_library_dir() / "matsumoto-nishimura-1998.yaml"
    tgt = tmp_path / "matsumoto-nishimura-1998.yaml"
    tgt.write_text(src.read_text())
    c = Catalog(tmp_path)
    c.load()
    before = c.paper("matsumoto-nishimura-1998").source_mtime
    time.sleep(0.02)
    os.utime(tgt, None)
    c.reload_if_stale()
    after = c.paper("matsumoto-nishimura-1998").source_mtime
    assert after >= before
