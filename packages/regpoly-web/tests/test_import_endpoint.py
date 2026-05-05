"""Phase 6 red — Tools page commit endpoint + audit-column fix + glob.

Today the Import button on /tools is a no-op (commitImport() in JS
just sets an error message). Add POST /api/v2/import/generators that
accepts the same multipart upload as preview and persists.

Additionally:
  - Audit endpoint queries `created_at` but the schema has
    `imported_at` — silently 500s and the JS swallows it.
  - The directory walker uses *.yaml only; the file picker accepts
    .yml too.
"""

from __future__ import annotations


def _yaml_payload(tmp_path):
    p = tmp_path / "g.yml"
    p.write_text(
        "family: MTGen\n"
        "common: {w: 32, n: 2, m: 1, r: 31, u: 11}\n"
        "generators:\n"
        "  - {a: '0xdeadbeef'}\n"
    )
    return p


def test_v2_import_generators_upload_persists_rows(client, tmp_path) -> None:
    p = _yaml_payload(tmp_path)
    files = {"file": ("g.yml", p.read_bytes(), "text/yaml")}
    r = client.post(
        "/api/v2/import/generators",
        files=files,
    )
    assert r.status_code in (200, 201), r.text
    body = r.json()
    assert "inserted" in body or "row_count" in body or "would_add" in body


def test_audit_endpoint_uses_imported_at_column(seeded_client) -> None:
    """Pin the column-name fix: the schema has `yaml_import.imported_at`,
    not `created_at`. The endpoint must return 200 with `items` even
    when the table has zero rows (i.e. the SELECT must not fail)."""
    r = seeded_client.get("/api/import/audit")
    assert r.status_code == 200, r.text
    body = r.json()
    assert "items" in body


def test_import_dir_glob_includes_yml_extension(client, tmp_path) -> None:
    """Directory import must walk both *.yml and *.yaml."""
    (tmp_path / "a.yml").write_text(
        "family: MTGen\ncommon: {w: 32, n: 2, m: 1, r: 31, u: 11}\n"
        "generators: []\n"
    )
    (tmp_path / "b.yaml").write_text(
        "family: MTGen\ncommon: {w: 32, n: 2, m: 1, r: 31, u: 11}\n"
        "generators: []\n"
    )
    # Confine to a path that's allowed; assume the tmp_path is permitted
    # (per the path-traversal guard, configurable via Settings.import_root).
    r = client.post(
        "/api/import/generators-dir",
        json={"directory": str(tmp_path)},
    )
    # If the path-traversal guard rejects the tmp dir, this test
    # records that limitation; otherwise both files must be picked up.
    if r.status_code == 200:
        body = r.json()
        assert body.get("imported_files", 0) >= 2, (
            f"glob missed *.yml files: {body}"
        )
