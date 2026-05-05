"""Phase 5 red — Tools page (Import / Export / Imports)."""

from __future__ import annotations


def test_sidebar_has_tools_link(seeded_client) -> None:
    r = seeded_client.get("/")
    assert r.status_code == 200
    body = r.text
    assert 'href="/tools"' in body, "Sidebar must link to /tools"


def test_tools_page_renders_three_tabs(seeded_client) -> None:
    r = seeded_client.get("/tools")
    assert r.status_code == 200
    body = r.text
    for tab in ("Import", "Export", "Imports"):
        assert tab in body, f"Tools page must mention {tab} tab"


def test_imports_tab_lists_yaml_imports(seeded_client) -> None:
    r = seeded_client.get("/tools?tab=imports")
    assert r.status_code == 200
