"""F2w-family `_effective_L` is min(k, 32) when the form leaves L blank."""
from regpoly_web.models import PrimitiveSearchCreate
from regpoly_web.routes.primitive_search import _effective_L


def _f2w(family, w, r, **extras) -> PrimitiveSearchCreate:
    return PrimitiveSearchCreate(
        family=family, L=0,
        structural_params={"w": w, "r": r,
                           "nb_terms": 3, "normal_basis": False,
                           "step": 1, **extras},
        fixed_params={},
    )


def test_f2w_lfsr_caps_L_at_32_for_paper_sized_state() -> None:
    # w=32, r=3 ⇒ k=96; L should be 32, not 96.
    body = _f2w("F2wLFSRGen", 32, 3)
    assert _effective_L(body, k=96) == 32


def test_f2w_polylcg_caps_L_at_32() -> None:
    body = _f2w("F2wPolyLCGGen", 32, 25)
    assert _effective_L(body, k=800) == 32


def test_f2w_small_state_uses_full_k() -> None:
    # w=8, r=3 ⇒ k=24 < 32; L should be 24.
    body = _f2w("F2wLFSRGen", 8, 3)
    assert _effective_L(body, k=24) == 24


def test_f2w_explicit_L_wins() -> None:
    body = _f2w("F2wLFSRGen", 32, 3)
    body.L = 16
    assert _effective_L(body, k=96) == 16
