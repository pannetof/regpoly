"""Phase 0 smoke test for regpoly-web.

Asserts the package is importable. Phase 5 replaces this with the full FastAPI
TestClient suite per the v2.0 plan.
"""

from __future__ import annotations


def test_package_imports() -> None:
    import regpoly_web

    assert regpoly_web is not None
