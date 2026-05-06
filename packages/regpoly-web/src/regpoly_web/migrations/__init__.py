# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Schema migration package for regpoly-web.

Each module bumps the DB by one major version. Migrations are
invoked from `database.init_sync` when an existing DB file is opened
at an older version, so the web app upgrades user data in-place on
startup. Each migration is also runnable as a standalone CLI from
`packages/regpoly-web/scripts/migrate_<v>.py`.
"""
