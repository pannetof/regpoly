# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — three sparkline endpoints converge on a single helper.

Today there are three different bucket-downsampling implementations:
  - routes/v2/dashboard.py:_sparkline_for_run (uses i*n//points indexing)
  - routes/v2/searches.py:_history_points_for_run (same)
  - routes/v2/searches_history.py:_bucket_downsample (last-value-per-bucket)

Three slightly different math choices. Promote the
last-value-per-bucket variant to a shared helper and have all three
call sites import it.
"""

from __future__ import annotations


def test_shared_downsampler_module_exists() -> None:
    from regpoly_web.routes.v2._downsample import bucket_downsample

    assert callable(bucket_downsample)


def test_shared_downsampler_returns_exact_K_points() -> None:
    from regpoly_web.routes.v2._downsample import bucket_downsample

    arr = list(range(100))
    out = bucket_downsample(arr, 10)
    assert len(out) == 10


def test_shared_downsampler_returns_input_when_smaller() -> None:
    from regpoly_web.routes.v2._downsample import bucket_downsample

    arr = [1, 2, 3]
    out = bucket_downsample(arr, 10)
    assert out == arr


def test_shared_downsampler_uses_last_value_per_bucket() -> None:
    from regpoly_web.routes.v2._downsample import bucket_downsample

    # 10 input → 5 buckets of 2 → take the last value of each pair.
    arr = list(range(10))  # 0..9
    out = bucket_downsample(arr, 5)
    assert len(out) == 5
    # Last value per bucket: arr[1], arr[3], arr[5], arr[7], arr[9]
    assert out == [1, 3, 5, 7, 9]
