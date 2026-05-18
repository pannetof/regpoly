# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""regpoly_legacy — optional pure-Python add-on for reading pre-v2 .dat
parameter files into regpoly Generator / Transformation objects.
"""

from regpoly_legacy.reader import LegacyReader
from regpoly_legacy.seek_factory import seek_from_legacy

__all__ = ["LegacyReader", "seek_from_legacy"]
