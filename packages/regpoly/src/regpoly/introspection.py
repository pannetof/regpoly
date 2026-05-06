# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
regpoly.introspection — read-only metadata helpers about generator
families and tempering types.

Phase 5.1: thin Python shim over the C++ factory specs. Lets the web
app and other consumers query "what params does family X take?" /
"can family Y be enumerated?" without importing _cpp directly.

Returns plain dicts so callers don't have to know about the underlying
ParamSpec C++ struct.
"""

from __future__ import annotations

from typing import Any

from regpoly_cpp import _regpoly_cpp as _cpp


def get_gen_param_specs(family: str) -> list[dict[str, Any]]:
    """Return the parameter specs for a generator family.

    Each spec is a dict with keys:
        name, type, structural, has_default, default,
        rand_type, rand_args, optimizable.
    """
    return _cpp.get_gen_param_specs(family)


def get_trans_param_specs(type_name: str) -> list[dict[str, Any]]:
    """Return the parameter specs for a transformation (tempering)
    type. Same dict shape as get_gen_param_specs."""
    return _cpp.get_trans_param_specs(type_name)


def family_is_enumerable(family: str) -> bool:
    """True iff the family has a finite enumerator (the search CLI's
    enumerate-mode supports it).
    """
    return bool(_cpp.family_is_enumerable(family))


__all__ = [
    "get_gen_param_specs",
    "get_trans_param_specs",
    "family_is_enumerable",
]
