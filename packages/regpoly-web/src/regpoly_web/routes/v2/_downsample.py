# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Shared bucket-downsampling helper.

Three v2 history endpoints (`dashboard.py`, `searches.py`,
`searches_history.py`) used to roll their own slightly-different
implementations. This module is the single canonical version: bucket
the input by index into K groups, take the last value of each bucket.

Properties:
- output length == min(len(values), k)
- preserves the final value (the bucket containing it picks it last)
- monotonically-decreasing series stay monotonically-decreasing
- O(k) work after the initial pass
"""

from __future__ import annotations

from typing import Sequence, TypeVar

T = TypeVar("T")


def bucket_downsample(values: Sequence[T], k: int) -> list[T]:
    if k <= 0 or not values:
        return []
    if len(values) <= k:
        return list(values)
    out: list[T] = []
    n = len(values)
    for i in range(k):
        start = (i * n) // k
        end = ((i + 1) * n) // k
        if end <= start:
            end = start + 1
        out.append(values[end - 1])
    return out
