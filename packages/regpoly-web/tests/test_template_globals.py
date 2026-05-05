"""Phase 1 red set — verify the redesigned base shell drops the
families globals and registers nav_items.

The pre-redesign base.html pinned the family hierarchy + paper list in
the sidebar via two template globals (`families_primary`,
`families_other`) and `library_papers`. The redesigned shell ships a
flat sidebar driven by a single `nav_items` global.
"""

from __future__ import annotations


def test_families_globals_removed(client) -> None:
    """The Jinja env exposed by the FastAPI app must NOT carry the
    pre-redesign family-hierarchy globals."""
    app = client.app
    templates = getattr(app.state, "templates", None)
    assert templates is not None, (
        "app.state.templates should be the Jinja2Templates instance"
    )
    globals_ = templates.env.globals
    assert "families_primary" not in globals_, (
        "families_primary must be dropped — sidebar no longer renders families"
    )
    assert "families_other" not in globals_, (
        "families_other must be dropped"
    )


def test_nav_items_global_registered(client) -> None:
    """The flat sidebar reads from one source of truth."""
    app = client.app
    templates = app.state.templates
    nav_items = templates.env.globals.get("nav_items")
    assert nav_items, "nav_items global must be registered for the flat sidebar"
    labels = [item.get("label") if isinstance(item, dict) else item
              for item in nav_items]
    expected = {"Dashboard", "Generators", "Tested generators",
                "Searches", "Library", "Tools"}
    assert set(labels) >= expected, (
        f"nav_items must include {expected}, got {labels}"
    )
