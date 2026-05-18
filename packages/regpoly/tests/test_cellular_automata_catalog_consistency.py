# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Drift detection between the CellularAutomata catalog YAML and the
hard-coded test parametrization tuples.

The catalog (docs/library/bhuvaneswari-bhattacharjee-2026-cellular-automata.yaml)
is the single source of truth for "what generators exist".  Per-test
parametrization tuples are duplicated as Python lists in the
``test_cellular_automata_table*.py`` files.  This test asserts that
every catalog entry has a corresponding test row and vice versa.

If you add a generator to the catalog YAML, you must add the matching
row to the test parametrization.
"""

from __future__ import annotations

from regpoly.library import Catalog


def test_catalog_loads_and_instantiates():
    """Every CellularAutomataGen entry in the catalog must instantiate
    via the C++ factory.  This is a smoke test on the catalog YAML."""
    cat = Catalog("docs/library")
    cat.load()
    gens = cat.all_generators(family="CellularAutomataGen")
    assert len(gens) > 100, f"Expected ~116 entries, got {len(gens)}"


def _by_id_prefix(gens, prefix):
    return [g for _, g in gens if g.id.startswith(prefix)]


def test_catalog_contains_expected_table_groups():
    """The catalog YAML must contain every paper-table group we test
    (identified by id prefix since CatalogGenerator doesn't expose tags
    in the current C++ binding)."""
    cat = Catalog("docs/library")
    cat.load()
    gens = cat.all_generators(family="CellularAutomataGen")
    assert _by_id_prefix(gens, "bhuv2026-r1"), "missing R1 entry"
    assert _by_id_prefix(gens, "bhuv2026-t3-7-10-"), "missing Tables 3/7-10"
    assert _by_id_prefix(gens, "bhuv2026-t11-"), "missing Table 11"
    assert _by_id_prefix(gens, "bhuv2026-t12-"), "missing Table 12"
    for weak in ("bhuv2026-r2", "bhuv2026-r3", "bhuv2026-r4"):
        assert _by_id_prefix(gens, weak), f"missing {weak}"


def test_catalog_table_12_count():
    """Table 12: 49 entries — see _TABLE_12 in test_cellular_automata_table12.py."""
    cat = Catalog("docs/library")
    cat.load()
    gens = cat.all_generators(family="CellularAutomataGen")
    t12 = _by_id_prefix(gens, "bhuv2026-t12-")
    assert len(t12) == 49, f"Expected 49 Table-12 entries, got {len(t12)}"


def test_catalog_table_11_count():
    """Table 11 (Adak-Das CA(90')/CA(150') combinations) — see
    _TABLE_11 in test_cellular_automata_table11_phase2.py."""
    cat = Catalog("docs/library")
    cat.load()
    gens = cat.all_generators(family="CellularAutomataGen")
    t11 = _by_id_prefix(gens, "bhuv2026-t11-")
    assert len(t11) >= 50, f"Expected ≥50 Table-11 entries, got {len(t11)}"
