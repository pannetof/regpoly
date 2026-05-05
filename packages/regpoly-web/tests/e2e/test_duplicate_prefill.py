"""Phase 6 red — Duplicate prefill end-to-end."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_primitive_form_consumes_s_param_query(page, base_url) -> None:
    page.goto(
        f"{base_url}/primitive-search?family=MTGen&L=64"
        "&s_w=32&s_n=2&s_m=1&s_r=31&s_u=11&max_tries=10",
        wait_until="networkidle",
    )
    val_w = page.locator("input[name='w'], [data-param='w'] input")
    if val_w.count():
        v = val_w.first.input_value()
        assert v == "32", f"primitive form ignored s_w URL param, got {v!r}"


def test_tempering_form_consumes_duplicate_from_session_storage(
    page, base_url,
) -> None:
    page.goto(base_url, wait_until="networkidle")
    page.evaluate(
        "() => sessionStorage.setItem("
        "  'regpoly.tempering.duplicate',"
        "  JSON.stringify({"
        "    cloned_from: 1, nb_tries: 5,"
        "    test_config: {L: 32}, components: []"
        "  })"
        ")"
    )
    page.goto(f"{base_url}/tempering-search?duplicate_from=1",
              wait_until="networkidle")
    nb_tries = page.locator("input[name='nb_tries'], [data-param='nb_tries']")
    if nb_tries.count():
        v = nb_tries.first.input_value()
        assert v == "5", (
            f"tempering form ignored duplicate_from sessionStorage, got {v!r}"
        )
