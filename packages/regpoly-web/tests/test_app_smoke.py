# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 5.3 smoke suite for regpoly-web.

Drives the FastAPI app through a TestClient against an ephemeral
SQLite DB. Phase 5.4+ will add Playwright golden paths.
"""

from __future__ import annotations


def test_package_imports() -> None:
    """Basic import sanity — kept from the Phase 0 stub so collection
    still works even if the FastAPI app fails to construct."""
    import regpoly_web

    assert regpoly_web is not None


def test_app_constructs(client) -> None:
    """The TestClient fixture exercises the lifespan startup. If this
    test sets up the fixture without raising, the app + catalog +
    worker pool all initialised cleanly."""
    assert client is not None


def test_families_index(client) -> None:
    r = client.get("/api/families")
    assert r.status_code == 200
    families = r.json()
    assert isinstance(families, list)
    assert len(families) > 0
    names = {f["name"] for f in families}
    # Every shipped family must be visible. Sanity-check a representative
    # sample (full list is in routes/families.py).
    assert {"MTGen", "TauswortheGen", "TGFSRGen", "MELGGen"} <= names


def test_family_detail_roundtrip(client) -> None:
    r = client.get("/api/families/MTGen")
    assert r.status_code == 200
    body = r.json()
    assert body["name"] == "MTGen"
    assert isinstance(body["params"], list)
    assert len(body["params"]) > 0
    # MT family is not enumerable.
    assert "enumerable" in body


def test_family_detail_unknown_404(client) -> None:
    r = client.get("/api/families/NoSuchFamily")
    assert r.status_code == 404


def test_transformations_index(client) -> None:
    r = client.get("/api/transformations")
    assert r.status_code == 200
    body = r.json()
    names = {t["name"] for t in body}
    assert {"tempMK", "tempMK2", "permut", "laggedTempering"} <= names


def test_transformation_detail(client) -> None:
    r = client.get("/api/transformations/tempMK")
    assert r.status_code == 200
    body = r.json()
    assert body["name"] == "tempMK"
    assert isinstance(body["params"], list)
    assert len(body["params"]) > 0


def test_transformation_detail_unknown_404(client) -> None:
    r = client.get("/api/transformations/no-such-trans")
    assert r.status_code == 404


def test_tests_index(client) -> None:
    r = client.get("/api/tests")
    assert r.status_code == 200
    body = r.json()
    names = {t["name"] for t in body}
    assert {"equidistribution", "collision_free", "tuplets"} == names


def test_library_lists_panneton_lecuyer_2005_xorshift(client) -> None:
    """``docs/library/panneton-lecuyer-2005-xorshift.yaml`` registers the
    Panneton & L'Ecuyer 2005 paper with all 37 Table III / IV entries.
    The catalog endpoint must surface the paper, and a sample generator
    detail must round-trip with the expected MarsaXorshiftGen family."""
    r = client.get("/api/library/papers")
    assert r.status_code == 200
    papers = r.json()
    match = [p for p in papers if p.get("id") == "panneton-lecuyer-2005-xorshift"]
    assert match, "Panneton & L'Ecuyer 2005 paper not in /api/library/papers"
    assert match[0]["year"] == 2005

    r = client.get("/api/library/papers/panneton-lecuyer-2005-xorshift")
    assert r.status_code == 200
    body = r.json()
    gens = body.get("generators", [])
    assert len(gens) == 37, f"expected 37 generators, got {len(gens)}"
    gen_ids = {g["id"] for g in gens}
    assert "panneton-lecuyer-2005-tiii-01" in gen_ids
    assert "panneton-lecuyer-2005-tiii-22" in gen_ids
    assert "panneton-lecuyer-2005-tiv-38" in gen_ids

    r = client.get("/api/library/generators/panneton-lecuyer-2005-tiii-22")
    assert r.status_code == 200
    g = r.json()
    assert g["family"] == "MarsaXorshiftGen"


def test_primitive_search_detail_has_delete_button(client) -> None:
    """The detail page exposes a Delete button that calls the existing
    DELETE /api/primitive-searches/{id}/generators endpoint (which
    drops every found generator + the run row, then redirects to
    /searches). Render-only check; the API endpoint already has its
    own coverage."""
    r = client.get("/primitive-search/1")
    assert r.status_code == 200
    html = r.text
    assert "del()" in html
    assert "busyDelete" in html
    assert "btn-outline-danger" in html


def test_primitive_search_form_f2w_hides_nocoeff(client) -> None:
    """For F2w families the form renders the paper-notation layout:
    `nocoeff` is hidden (it's derived server-side from {m, t, q}), and
    {t, q, modM, coeff} are randomized per iteration by default."""
    r = client.get("/primitive-search?family=F2wLFSRGen")
    assert r.status_code == 200
    html = r.text
    # The F2w branch sets up an `isF2wFamily()` helper that drives the
    # `nocoeff` filter and the buildBody skip.
    assert "isF2wFamily" in html
    # The new rand_type and the paper-form params surface in the
    # rendered spec list.
    assert "irreducible_gf2" in html or "F2wLFSRGen" in html


def test_f2w_family_detail_paper_notation_specs(client) -> None:
    """The F2wLFSRGen ParamSpec list exposes paper-notation entry:
    new {nb_terms, t, q} params, modM samplable via `irreducible_gf2`."""
    r = client.get("/api/families/F2wLFSRGen")
    assert r.status_code == 200
    body = r.json()
    names = {p["name"] for p in body["params"]}
    assert {"nb_terms", "t", "q", "modM", "nocoeff", "coeff"} <= names
    assert "m" not in names
    specs = {p["name"]: p for p in body["params"]}
    assert specs["modM"]["rand_type"] == "irreducible_gf2"
    assert specs["t"]["rand_type"] == "range"
    assert specs["q"]["rand_type"] == "range"
    # nocoeff is derived from m/t/q in from_params, so the search loop
    # must skip it: has_default=True with no rand_type.
    assert specs["nocoeff"]["has_default"] is True
    assert specs["nocoeff"]["rand_type"] in ("", "none")


def test_primitive_search_detail_has_se_column_helpers(client) -> None:
    """The detail page exposes a conditional SE column rendered only
    when the run enabled the equidistribution post-filter."""
    r = client.get("/primitive-search/1")
    assert r.status_code == 200
    html = r.text
    # Alpine helpers that drive the column visibility.
    assert "hasSeFilter()" in html
    assert "formatSe(g)" in html


def test_generator_detail_has_gap_table(client) -> None:
    """The /generators/{id} page renders both the SVG bar chart AND a
    numeric dimension-gap table (mirrors /tested-generators/{id}).
    The table was missing earlier — equidStripes() was defined but
    never invoked in the markup."""
    r = client.get("/generators/1")
    assert r.status_code == 200
    html = r.text
    # Helper must be wired into the markup, not just defined.
    assert "equidStripes(gen.pis_gaps)" in html
    # Table chrome: matches the existing tested-generators rendering.
    assert "equidist-grid" in html
    assert "Dimension-gap table" in html


def test_searches_summary_surfaces_max_se(client) -> None:
    """The /searches table renders ``s.max_se`` in the Summary column
    when the equidistribution post-filter was used — same column users
    rely on to spot active filters at a glance."""
    r = client.get("/searches")
    assert r.status_code == 200
    html = r.text
    # Alpine wiring: the summary() helper threads max_se onto the
    # primitive-search summary string.
    assert "s.max_se" in html
    # The rendered prefix users actually see.
    assert "SE ≤" in html or "SE ≤ " in html or "SE \\u2264" in html \
        or "SE &le;" in html or "`SE ≤ " in html or "SE ≤" in html


def test_tools_page_has_workers_tab(client) -> None:
    """The Tools page exposes a Workers tab that polls /api/workers/status."""
    r = client.get("/tools")
    assert r.status_code == 200
    html = r.text
    # Tab pill is registered alongside Import / Export / Imports.
    assert 'data-tools-tab="workers"' in html
    assert ">Workers</a>" in html
    # Alpine wiring for the polling/render.
    assert "loadWorkers()" in html
    assert "startWorkersPolling()" in html
    assert "/api/workers/status" in html


def test_primitive_search_detail_columns_are_sortable(client) -> None:
    """Generator-table headers are click-to-sort: each sortable column
    wires `toggleSort(...)` and renders a `sortIndicator(...)` arrow."""
    r = client.get("/primitive-search/1")
    assert r.status_code == 200
    html = r.text
    # Helpers
    assert "toggleSort(col)" in html
    assert "sortIndicator(col)" in html
    # ID, k, and (conditional) SE are always sortable.
    for col in ("id", "k", "pis_se"):
        assert f"toggleSort('{col}')" in html, col
        assert f"sortIndicator('{col}')" in html, col
    # WELL m1/m2/m3 columns are sortable too.
    for col in ("m1", "m2", "m3"):
        assert f"toggleSort('{col}')" in html, col


def test_primitive_search_form_well_is_max_cost_only(client) -> None:
    """For WELL families the form exposes only the `max_cost` lever:
    the matrices editor is absent (server samples per iteration),
    paper presets are not advertised, and the max_cost input renders
    with a `required` attribute."""
    r = client.get("/primitive-search?family=WELLGen")
    assert r.status_code == 200
    html = r.text

    # max_cost is the (required) WELL search lever.
    assert "max_cost" in html
    assert "wellMaxCostIsValid" in html
    assert ">Max cost</label>" in html

    # The old per-slot editor and preset machinery are gone.
    assert "matricesSlots" not in html
    assert "applyPreset" not in html
    assert "materializeMatrices" not in html
    # Paper preset names must not leak into the rendered form.
    assert "WELL512a" not in html
    assert "WELL19937a" not in html
    assert "WELL44497a" not in html
