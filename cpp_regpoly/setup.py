"""
setup.py — Build the C++ extension module via pybind11.
"""

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "regpoly._regpoly_cpp",
        sources=[
            "src/cpp/src/generateur.cpp",
            "src/cpp/src/gen_polylcg.cpp",
            "src/cpp/src/gen_tausworthe.cpp",
            "src/cpp/src/gen_tgfsr.cpp",
            "src/cpp/src/gen_mt.cpp",
            "src/cpp/src/gen_f2w_base.cpp",
            "src/cpp/src/gen_f2w_polylcg.cpp",
            "src/cpp/src/gen_f2w_lfsr.cpp",
            "src/cpp/src/gf2w_arith.cpp",
            "src/cpp/src/bitvect.cpp",
            "src/cpp/src/gauss.cpp",
            "src/cpp/src/params.cpp",
            "src/cpp/src/gen_matsumoto.cpp",
            "src/cpp/src/gen_marsaxorshift.cpp",
            "src/cpp/src/gen_ac1d.cpp",
            "src/cpp/src/gen_carry2.cpp",
            "src/cpp/src/factory.cpp",
            "src/cpp/src/trans_permutation.cpp",
            "src/cpp/src/trans_temper_mk.cpp",
            "src/cpp/src/poly_lattice.cpp",
            "src/cpp/src/lattice_polys.cpp",
            "src/cpp/src/bindings.cpp",
        ],
        include_dirs=["src/cpp/include"],
        libraries=["ntl", "pthread"],
        cxx_std=17,
        extra_compile_args=["-O3", "-march=native"],
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
