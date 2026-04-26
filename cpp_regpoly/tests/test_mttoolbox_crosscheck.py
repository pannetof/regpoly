"""Cross-check: regpoly d(1..5) ≡ MTToolBox calc_equidist d(1..5).

Loads the per-family parameter catalogs from
``docs/library/*_params.yaml``, builds the corresponding regpoly
``Generateur``, computes d(1..5) via the appropriate equidistribution
method, and asserts equality with the reference values captured by
``tests/_gen_mttoolbox_reference.py`` from MTToolBox's binaries.

Marker conventions:
- ``@pytest.mark.slow`` for entries with ``mexp > slow_threshold_mexp``.
  Hidden from the default ``pytest`` run by ``pyproject.toml``'s
  ``addopts = "-m 'not slow'"``; opt in with ``pytest -m slow``.
- Catalogs marked ``deferred: true`` are skipped at collection time
  with the reason exposed via ``pytest.skip``.

Run:
    pytest tests/test_mttoolbox_crosscheck.py             # fast lane
    pytest -m slow tests/test_mttoolbox_crosscheck.py     # slow lane
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

from regpoly.generateur import Generateur
from regpoly.combinaison import Combinaison
from regpoly.transformation import Transformation
from regpoly.analyses.equidistribution_test import (
    EquidistributionTest,
    METHOD_NOTPRIMITIVE,
    METHOD_SIMD_NOTPRIMITIVE,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_LIB = REPO_ROOT / "docs" / "library"

MAX_V = 5  # cross-check granularity (per user spec)


# Map family-name → which equidistribution method drives the test.
_METHOD_FOR_FAMILY = {
    "SFMT": METHOD_SIMD_NOTPRIMITIVE,
    "dSFMT": METHOD_SIMD_NOTPRIMITIVE,
    "MTGP": METHOD_SIMD_NOTPRIMITIVE,
    # Default for everything else: METHOD_NOTPRIMITIVE.
}


def _method_for(family: str) -> int:
    return _METHOD_FOR_FAMILY.get(family, METHOD_NOTPRIMITIVE)


def _family_for_create(family: str) -> str:
    """Map yaml `family` field to the string Generateur.create() accepts."""
    aliases = {
        "MersenneTwister": "MT",
        "XORSHIFT128": "xorshift",
        "TinyMT32": "tinymt",
        "RMT64": "rmt",
        "dSFMT": "dsfmt",
    }
    return aliases.get(family, family)


def _load_reference(catalog_stem: str) -> dict[str, list[int]]:
    """Load the d(1..5) reference dict for a family.

    catalog_stem = "sfmt" → tests/sfmt_mttoolbox_d5.py → SFMT_MTTOOLBOX_D5.
    Returns an empty dict if the reference module hasn't been generated
    yet (the cross-check entries will pytest.skip in that case).
    """
    ref_path = REPO_ROOT / "tests" / f"{catalog_stem}_mttoolbox_d5.py"
    if not ref_path.exists():
        return {}
    ns: dict = {}
    exec(ref_path.read_text(), ns)
    for key in ns:
        if key.endswith("_MTTOOLBOX_D5"):
            return ns[key]
    return {}


def _collect_cases() -> list[tuple]:
    """Scan all docs/library/*_params.yaml and yield one case per generator."""
    cases: list[tuple] = []
    for catalog_path in sorted(DOCS_LIB.glob("*_params.yaml")):
        with catalog_path.open() as f:
            catalog = yaml.safe_load(f)
        family = catalog["family"]
        deferred = catalog.get("deferred", False)
        deferred_reason = catalog.get("deferred_reason", "deferred")
        threshold = catalog.get("slow_threshold_mexp", 20000)
        catalog_stem = catalog_path.stem.replace("_params", "")
        ref = _load_reference(catalog_stem)
        for entry in catalog.get("generators", []):
            cases.append({
                "family": family,
                "catalog_stem": catalog_stem,
                "id": entry["id"],
                "mexp": entry["mexp"],
                "regpoly_params": entry.get("regpoly", {}),
                "expected_d5": ref.get(entry["id"]),
                "deferred": deferred,
                "deferred_reason": deferred_reason,
                "slow": entry["mexp"] > threshold,
                "xfail_reason": entry.get("cross_check_xfail"),
                "tempering": entry.get("tempering", []),
            })
    return cases


def _make_param(case: dict) -> pytest.param:
    marks = []
    if case["slow"]:
        marks.append(pytest.mark.slow)
    if case.get("xfail_reason"):
        marks.append(pytest.mark.xfail(reason=case["xfail_reason"], strict=False))
    return pytest.param(case, id=f"{case['catalog_stem']}-{case['id']}",
                        marks=marks)


_CASES = _collect_cases()


@pytest.mark.parametrize("case", [_make_param(c) for c in _CASES])
def test_regpoly_matches_mttoolbox_d5(case: dict) -> None:
    if case["deferred"]:
        pytest.skip(f"deferred: {case['deferred_reason']}")
    if case["expected_d5"] is None:
        pytest.skip(
            f"no reference data for {case['family']} {case['id']} — "
            f"run `python tests/_gen_mttoolbox_reference.py "
            f"{case['catalog_stem']}` to capture"
        )

    family_alias = _family_for_create(case["family"])
    L = 32  # 32-bit accuracy matches calc_equidist for SFMT/MT/XORSHIFT
    if case["family"] == "dSFMT":
        L = 52
    elif case["family"] == "RMT64":
        L = 64
    # Pass `mexp` from the entry top-level into params (most generators
    # take it as a structural parameter; xorshift doesn't care).
    full_params = dict(case["regpoly_params"])
    if case["mexp"] is not None and "mexp" not in full_params and \
       case["family"] not in ("XORSHIFT128", "TinyMT32"):
        full_params["mexp"] = case["mexp"]
    gen = Generateur.create(family_alias, L=L, **full_params)
    # Build any tempering chain declared in the catalog (e.g. MT19937's
    # standard 4-stage temper).  Each entry is {type: <str>, params: <dict>}.
    temper_chain = [
        Transformation.create(t["type"], **t.get("params", {}))
        for t in case.get("tempering", [])
    ]
    C = Combinaison.CreateFromFiles([[gen]], Lmax=L, temperings=[temper_chain])
    C.reset()

    method = _method_for(case["family"])
    test = EquidistributionTest(
        L=MAX_V,
        delta=[sys.maxsize] * (MAX_V + 1),
        mse=sys.maxsize,
        method=method,
    )
    result = test.run(C)
    actual_d5 = list(result.ecart[1:MAX_V + 1])
    expected_d5 = case["expected_d5"]
    assert actual_d5 == expected_d5, (
        f"{case['family']} {case['id']} (mexp={case['mexp']}): "
        f"d(1..5) regpoly={actual_d5} vs MTToolBox={expected_d5}"
    )
