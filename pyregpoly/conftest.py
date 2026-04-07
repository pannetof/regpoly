"""
conftest.py — shared pytest configuration.

Adds --sizes command-line option so both test files pick up the sizes to test.

Usage:
    pytest test_bitvect.py test_matrix.py --sizes 8,16,32,64,94
"""


def pytest_addoption(parser):
    parser.addoption(
        "--sizes",
        default="8",
        help="Comma-separated bit-vector/matrix sizes to test, e.g. 8,16,32,64,94",
    )


def pytest_generate_tests(metafunc):
    if "size" in metafunc.fixturenames:
        raw = metafunc.config.getoption("--sizes")
        sizes = [int(s.strip()) for s in raw.split(",")]
        metafunc.parametrize("size", sizes)
