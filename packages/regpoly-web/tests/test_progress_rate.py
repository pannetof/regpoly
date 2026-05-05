"""Phase 3 red — Rolling-rate helper, keyed by run_id.

The helper is a 5-second rolling window over (timestamp, tries) samples
exposed via `tasks/_progress_rate.py::RollingRate`. Each progress write
calls `RollingRate.for_run(run_id).observe(now, tries)` then reads
`.rate()`. Different runs must NOT cross-contaminate.
"""

from __future__ import annotations


def test_rolling_rate_window_basic() -> None:
    from regpoly_web.tasks._progress_rate import RollingRate

    r = RollingRate(window_sec=5.0)
    r.observe(t=0.0, tries=0)
    r.observe(t=1.0, tries=100)
    r.observe(t=2.0, tries=200)
    # 200 tries in 2 s → 100/s.
    assert 90.0 <= r.rate(now=2.0) <= 110.0


def test_rolling_rate_drops_old_samples() -> None:
    from regpoly_web.tasks._progress_rate import RollingRate

    r = RollingRate(window_sec=5.0)
    r.observe(t=0.0, tries=0)
    r.observe(t=10.0, tries=10_000)  # only this sample is within window
    r.observe(t=11.0, tries=10_100)
    rate = r.rate(now=11.0)
    # Old sample at t=0 dropped; the window covers [6, 11] s.
    # Effective rate computed over (10→11) ≈ 100/s, not 1000/s.
    assert 90.0 <= rate <= 200.0


def test_rolling_rate_keyed_by_run_id_no_cross_contamination() -> None:
    from regpoly_web.tasks._progress_rate import for_run, drop_run

    a = for_run(101, window_sec=5.0)
    b = for_run(102, window_sec=5.0)
    assert a is not b
    a.observe(t=0.0, tries=0)
    a.observe(t=1.0, tries=1000)
    b.observe(t=0.0, tries=0)
    b.observe(t=1.0, tries=10)
    assert a.rate(now=1.0) > 100.0
    assert b.rate(now=1.0) < 100.0
    drop_run(101)
    drop_run(102)
