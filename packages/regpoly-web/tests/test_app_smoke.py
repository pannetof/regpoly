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
