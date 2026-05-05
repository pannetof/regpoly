"""Phase 5 red — Dry-run import preview."""

from __future__ import annotations


def test_preview_returns_envelope(seeded_client, tmp_path) -> None:
    yml = tmp_path / "g.yml"
    yml.write_text(
        "generators:\n"
        "  - family: MTGen\n"
        "    L: 64\n"
        "    structural_params: {w: 32, n: 2, m: 1, r: 31, u: 11}\n"
        "    search_params: {a: '0xdeadbeef'}\n"
    )
    files = {"file": ("g.yml", yml.read_bytes(), "text/yaml")}
    r = seeded_client.post(
        "/api/v2/import/generators/preview",
        files=files,
    )
    assert r.status_code == 200, r.text
    body = r.json()
    assert "would_add" in body
    assert "would_skip" in body
    assert "conflicts" in body
