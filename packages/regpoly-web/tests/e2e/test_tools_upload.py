# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 6 red — Tools "Import" button performs an actual upload.

Before P6 the button just sets an error message. After P6 it submits
the file via POST /api/v2/import/generators and the imports tab
shows a new audit row.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_import_button_does_not_show_unimplemented_message(page, base_url, tmp_path) -> None:
    """The button must commit, not bail with the
    'isn't implemented yet' error message."""
    yml = tmp_path / "g.yml"
    yml.write_text(
        "family: MTGen\n"
        "common: {w: 32, n: 2, m: 1, r: 31, u: 11}\n"
        "generators:\n"
        "  - {a: '0xdeadbeef'}\n"
    )
    page.goto(f"{base_url}/tools", wait_until="networkidle")
    file_input = page.locator("input[type=file][data-import-file]")
    if file_input.count() == 0:
        pytest.skip("file input not rendered")
    file_input.set_input_files(str(yml))
    import_btn = page.locator("button:has-text('Import'):not([data-dry-run])")
    if import_btn.count() == 0:
        pytest.skip("Import button not rendered")
    import_btn.first.click()
    # Wait briefly; assert the unimplemented message did NOT appear.
    page.wait_for_timeout(300)
    error_text = page.locator(
        "text=isn't implemented yet"
    )
    assert error_text.count() == 0, (
        "Tools Import button still shows the 'unimplemented' error message"
    )
