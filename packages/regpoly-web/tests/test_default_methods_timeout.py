# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Regression: ``_default_test_methods_cached`` must not hang the web.

``is_full_period()`` in C++ trial-divides 2^K-1 for non-Mersenne-prime
K. For K=800 (paper-form F2w) the call can take tens of seconds to
minutes — long enough to make the FastAPI event loop unresponsive
and the library page render as "nothing at all".

We wrap the call in ``asyncio.wait_for(timeout=N)`` so on overflow we
cache a ``None`` fallback. This test injects a slow stub for
``_default_test_methods_for`` and asserts:

  1. The async wrapper returns within the timeout budget.
  2. The result is the documented all-``None`` fallback.
  3. The fallback is cached, so a second call returns instantly.
"""

from __future__ import annotations

import asyncio
import time
import types

import pytest


@pytest.mark.asyncio
async def test_default_methods_times_out_and_caches_fallback(monkeypatch):
    from regpoly_web.routes import library as lib

    # The wrapper reads REGPOLY_DEFAULT_METHODS_TIMEOUT at import time
    # into the module-level constant. Override the constant directly so
    # this test can run in seconds, not five.
    monkeypatch.setattr(lib, "_DEFAULT_METHODS_TIMEOUT_SECONDS", 0.5)

    def slow_default_methods(g):  # pragma: no cover — runs to timeout
        # Sleep with a busy loop so the threadpool worker doesn't
        # cooperate with cancellation — mirrors the real C++ call
        # holding the GIL with no checkpoints.
        end = time.monotonic() + 5.0
        while time.monotonic() < end:
            pass
        return {"equidistribution": "harase"}

    monkeypatch.setattr(lib, "_default_test_methods_for", slow_default_methods)

    fake_request = types.SimpleNamespace(
        app=types.SimpleNamespace(state=types.SimpleNamespace())
    )
    fake_paper = types.SimpleNamespace(source_mtime=42.0)
    fake_gen = types.SimpleNamespace(id="poison-pill", Lmax=800)

    t0 = time.monotonic()
    result1 = await lib._default_test_methods_cached(
        fake_request, fake_paper, fake_gen,
    )
    elapsed1 = time.monotonic() - t0
    assert elapsed1 < 2.0, f"timeout fallback didn't fire (took {elapsed1:.2f}s)"
    assert all(v is None for v in result1.values()), result1

    # Negative result is cached → second call is instant and doesn't
    # spawn another doomed thread.
    t1 = time.monotonic()
    result2 = await lib._default_test_methods_cached(
        fake_request, fake_paper, fake_gen,
    )
    elapsed2 = time.monotonic() - t1
    assert elapsed2 < 0.05, f"cache miss on second call ({elapsed2:.3f}s)"
    assert result2 == result1
