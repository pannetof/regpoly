"""
setup.py — Build the C++ extension module via pybind11.
"""

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "regpoly._regpoly_cpp",
        sources=[
            # core
            "src/cpp/core/generator.cpp",
            "src/cpp/core/params.cpp",
            "src/cpp/core/factory.cpp",
            "src/cpp/core/gen_enumerator.cpp",
            # algebra
            "src/cpp/algebra/bitvect.cpp",
            "src/cpp/algebra/gauss.cpp",
            "src/cpp/algebra/gf2w_arith.cpp",
            "src/cpp/algebra/bm.cpp",
            # generators
            "src/cpp/generators/polylcg.cpp",
            "src/cpp/generators/tausworthe.cpp",
            "src/cpp/generators/tgfsr.cpp",
            "src/cpp/generators/mt.cpp",
            "src/cpp/generators/f2w_base.cpp",
            "src/cpp/generators/f2w_polylcg.cpp",
            "src/cpp/generators/f2w_lfsr.cpp",
            "src/cpp/generators/matsumoto.cpp",
            "src/cpp/generators/marsaxorshift.cpp",
            "src/cpp/generators/ac1d.cpp",
            "src/cpp/generators/well.cpp",
            "src/cpp/generators/melg.cpp",
            "src/cpp/generators/sfmt.cpp",
            "src/cpp/generators/dsfmt.cpp",
            "src/cpp/generators/mtgp.cpp",
            "src/cpp/generators/xorshift128.cpp",
            "src/cpp/generators/tinymt32.cpp",
            "src/cpp/generators/rmt64.cpp",
            # transforms
            "src/cpp/transforms/permutation.cpp",
            "src/cpp/transforms/temper_mk.cpp",
            "src/cpp/transforms/lag_mask.cpp",
            # lattice
            "src/cpp/lattice/dual_lattice.cpp",
            "src/cpp/lattice/me_helpers.cpp",
            "src/cpp/lattice/me_harase.cpp",
            "src/cpp/lattice/me_notprimitive.cpp",
            "src/cpp/lattice/me_notprimitive_simd.cpp",
            "src/cpp/lattice/temper_optimizer.cpp",
            # pybind glue
            "src/cpp/bindings/bindings.cpp",
        ],
        include_dirs=[
            "src/cpp/include/core",
            "src/cpp/include/algebra",
            "src/cpp/include/generators",
            "src/cpp/include/transforms",
            "src/cpp/include/lattice",
        ],
        libraries=["ntl", "pthread"],
        cxx_std=17,
        extra_compile_args=["-O3", "-march=native"],
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
