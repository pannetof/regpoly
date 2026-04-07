"""
setup.py — Build Cython extensions if Cython is available.

Generated .c files go into build/cython/ to keep source directories clean.
The .so files are placed by setuptools into the source tree for editable installs.

When Cython is not installed, the package still works as pure Python.
"""

from setuptools import setup

CYTHON_SOURCES = [
    "src/regpoly/bitvect.py",
    "src/regpoly/generateur.py",
    "src/regpoly/generators/polylcg.py",
    "src/regpoly/generators/tausworthe.py",
    "src/regpoly/generators/tgfsr.py",
    "src/regpoly/generators/mt.py",
    "src/regpoly/generators/genf2w.py",
    "src/regpoly/gauss_matrix.py",
]

try:
    from Cython.Build import cythonize
    ext_modules = cythonize(
        CYTHON_SOURCES,
        build_dir="build/cython",
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
        },
    )
except ImportError:
    ext_modules = []

setup(ext_modules=ext_modules)
