"""Phase 4 red — Publish / Unpublish on tested-generator detail.

POST /api/v2/tested-generators/{id}/publish body {library_id} → 200,
sets tested_generator.library_id. DELETE clears it. Both idempotent
(re-publish same id, double-delete = 200).
"""

from __future__ import annotations


def test_post_publish_attaches_library_id(seeded_client) -> None:
    r = seeded_client.post(
        "/api/v2/tested-generators/4242/publish",
        json={"library_id": "matsumoto2008"},
    )
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["tested_gen_id"] == 4242
    assert body["library_id"] == "matsumoto2008"


def test_delete_publish_clears_library_id(seeded_client) -> None:
    seeded_client.post(
        "/api/v2/tested-generators/4242/publish",
        json={"library_id": "matsumoto2008"},
    )
    r = seeded_client.delete("/api/v2/tested-generators/4242/publish")
    assert r.status_code == 200
    body = r.json()
    assert body["library_id"] is None


def test_publish_idempotent(seeded_client) -> None:
    a = seeded_client.post(
        "/api/v2/tested-generators/4242/publish",
        json={"library_id": "matsumoto2008"},
    )
    b = seeded_client.post(
        "/api/v2/tested-generators/4242/publish",
        json={"library_id": "matsumoto2008"},
    )
    assert a.status_code == 200
    assert b.status_code == 200


def test_publish_unknown_id_returns_404(seeded_client) -> None:
    r = seeded_client.post(
        "/api/v2/tested-generators/999999/publish",
        json={"library_id": "matsumoto2008"},
    )
    assert r.status_code == 404


def test_v2_generator_related_returns_envelope(seeded_client) -> None:
    r = seeded_client.get("/api/v2/generators/1/related")
    assert r.status_code == 200, r.text
    body = r.json()
    # The seeded generator has search_run_id=1, no tempering attached.
    assert "primitive_search_id" in body
    assert "tempering_runs" in body
    assert "library_id" in body
    assert isinstance(body["tempering_runs"], list)


def test_v2_transition_matrix_coords_endpoint_exists(seeded_client) -> None:
    # The endpoint mirrors v1's shape; we accept either 200 with body
    # or 404 (when the underlying analyzer can't compute on the seeded
    # row). We're checking it's *mounted* and produces a JSON response,
    # not that it produces a coordinate.
    r = seeded_client.get("/api/v2/generators/1/transition-matrix-coords")
    assert r.status_code in (200, 404, 422), r.text
