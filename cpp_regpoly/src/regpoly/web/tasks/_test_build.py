"""Shared factory for the three test types used by the tempering-search
worker and the published-generators library runner."""

from __future__ import annotations


def build_test(test_config: dict, Lmax: int):
    """Instantiate the test class matching ``test_config['type']``."""
    from regpoly.analyses.equidistribution_test import EquidistributionTest
    from regpoly.analyses.collision_free_test import CollisionFreeTest
    from regpoly.analyses.tuplets_test import TupletsTest

    dispatch = {
        "equidistribution": EquidistributionTest,
        "collision_free":   CollisionFreeTest,
        "tuplets":          TupletsTest,
    }
    ttype = test_config.get("type", "equidistribution")
    cls = dispatch.get(ttype)
    if cls is None:
        raise ValueError(f"Unknown test type: {ttype}")
    return cls._from_params(test_config, Lmax)


__all__ = ["build_test"]
