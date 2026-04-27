"""Display formatting for generator / transformation parameters.

All bitmask-like fields (anything the C++ ``ParamSpec`` flags as
``rand_type="bitmask"``, ``rand_type="bitmask_vec"``, or
``type="uint_vec"``) are rendered as hexadecimal strings so every
template shows them in the same ``0xHEX`` form without needing its
own heuristic.

The formatter caches the per-family / per-transformation key sets so
it can be called on every API payload without paying a C++ spec
lookup each time.
"""

from __future__ import annotations

from functools import lru_cache

import regpoly._regpoly_cpp as _cpp


HEX_MASK_64 = (1 << 64) - 1


def _hex(value: int) -> str:
    """Canonical hex rendering for a bitmask-shaped integer.

    Treats the stored value as unsigned (masks with 64 bits) so that
    negative ``int64_t`` values stored from a ``uint64_t`` bitmask
    round-trip to the same hex representation.  Uppercase,
    ``0x``-prefixed, no zero-padding.
    """
    if not isinstance(value, int) or isinstance(value, bool):
        return str(value)
    return f"0x{value & HEX_MASK_64:X}"


@lru_cache(maxsize=64)
def _gen_bitmask_keys(family: str) -> tuple[frozenset[str], frozenset[str]]:
    """Return ``(scalar_hex_keys, vector_hex_keys)`` for ``family``."""
    try:
        specs = _cpp.get_gen_param_specs(family)
    except Exception:
        return frozenset(), frozenset()
    scalar, vec = set(), set()
    for s in specs:
        rt = s.get("rand_type", "") or ""
        ty = s.get("type", "") or ""
        if rt == "bitmask":
            scalar.add(s["name"])
        elif rt == "bitmask_vec" or ty == "uint_vec":
            vec.add(s["name"])
    return frozenset(scalar), frozenset(vec)


@lru_cache(maxsize=32)
def _trans_bitmask_keys(type_name: str) -> tuple[frozenset[str], frozenset[str]]:
    try:
        specs = _cpp.get_trans_param_specs(type_name)
    except Exception:
        return frozenset(), frozenset()
    scalar, vec = set(), set()
    for s in specs:
        rt = s.get("rand_type", "") or ""
        ty = s.get("type", "") or ""
        if rt == "bitmask":
            scalar.add(s["name"])
        elif rt == "bitmask_vec" or ty == "uint_vec":
            vec.add(s["name"])
    return frozenset(scalar), frozenset(vec)


def _apply(params: dict, scalar: frozenset[str],
           vec: frozenset[str]) -> dict:
    if not params:
        return params
    out = dict(params)
    for k in scalar:
        v = out.get(k)
        if isinstance(v, int) and not isinstance(v, bool):
            out[k] = _hex(v)
    for k in vec:
        v = out.get(k)
        if isinstance(v, list):
            out[k] = [_hex(x) if isinstance(x, int)
                      and not isinstance(x, bool) else x
                      for x in v]
    return out


def format_gen_params(family: str, params: dict | None) -> dict:
    """Return a shallow copy of ``params`` with bitmask scalars /
    vectors rendered as hex strings, using the ``family``'s ParamSpec
    as the source of truth for which keys are bitmask-shaped."""
    if not params:
        return params or {}
    scalar, vec = _gen_bitmask_keys(family)
    return _apply(params, scalar, vec)


def format_trans_params(type_name: str,
                        params: dict | None) -> dict:
    """Same as :func:`format_gen_params` but for transformation
    params.  ``params`` is expected to include the ``type`` key; it
    is preserved unchanged."""
    if not params:
        return params or {}
    scalar, vec = _trans_bitmask_keys(type_name)
    return _apply(params, scalar, vec)


def format_tempering_list(tempering: list[dict] | None) -> list[dict]:
    """Format a list of ``{type, **params}`` tempering steps."""
    if not tempering:
        return []
    return [format_trans_params(t.get("type", ""), t) for t in tempering]


__all__ = [
    "format_gen_params",
    "format_trans_params",
    "format_tempering_list",
]
