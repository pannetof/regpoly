# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — tempering worker heartbeat.

Today `tasks/tempering.py:_search_one_combo` only writes a progress
row when a NEW best is found inside the inner loop. A long-grinding
combo with no improvement looks frozen on the live UI for minutes.

Pin the contract: the worker must emit a progress row at least every
N seconds even with no new best. We test by importing the helper that
decides when to write — it should accept a `last_heartbeat_at` arg
and return True after `heartbeat_interval` seconds elapse.
"""

from __future__ import annotations


def test_tempering_heartbeat_helper_exists() -> None:
    from regpoly_web.tasks.tempering import should_emit_heartbeat

    assert callable(should_emit_heartbeat)


def test_tempering_heartbeat_fires_after_interval() -> None:
    from regpoly_web.tasks.tempering import should_emit_heartbeat

    assert should_emit_heartbeat(now=1000.0, last_emit=900.0,
                                 interval=5.0) is True


def test_tempering_heartbeat_holds_within_interval() -> None:
    from regpoly_web.tasks.tempering import should_emit_heartbeat

    assert should_emit_heartbeat(now=1002.0, last_emit=1001.0,
                                 interval=5.0) is False
