"""
test_base.py — Abstract base class for ME/CF tests.

Provides the shared _prepare_mat() and _prepare_mat_packed() static methods
used by EquidistributionTest, CollisionFreeTest, and TupletsTest.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from regpoly.bitvect import BitVect
from regpoly.matrix import Matrix
from regpoly.gauss_matrix import PyIntGaussMatrix, PackedGaussMatrix
from regpoly.tests.test_results_base import AbstractTestResults

if TYPE_CHECKING:
    from regpoly.combinaison import Combinaison

# Default backend: "packed" for uint64 word arrays, "pyint" for Python ints
_DEFAULT_BACKEND = "packed"


class AbstractTest(ABC):
    """
    Abstract base for all statistical test classes.

    A test object holds the configuration parameters.  Calling run()
    executes the test against the current combination and returns an
    AbstractTestResults object.
    """

    @abstractmethod
    def run(self, C: "Combinaison", *args, **kwargs) -> AbstractTestResults:
        """Execute the test on C and return a results object."""

    @staticmethod
    def _prepare_mat(C: "Combinaison", indice_max: int) -> Matrix:
        """
        PrepareMat: build the k_g × indice_max generator matrix as a Matrix.

        Used by transition_matrix and other callers that need the Matrix object.
        """
        L   = C.L
        mat = Matrix(C.k_g, L, indice_max)
        ml  = 0
        for j, comp in enumerate(C.components):
            gen = C[j]
            k   = gen.k
            bc  = BitVect.zeros(k)
            bc.put_bit(0, 1)
            for _gl in range(k):
                state    = gen.initialize_state(bc)
                bc     >>= 1
                mat[ml, 0] = comp.transform(state)
                for mc in range(1, indice_max):
                    state      = next(gen)
                    mat[ml, mc] = comp.transform(state)
                ml += 1
        return mat

    @staticmethod
    def _prepare_mat_packed(
        C: "Combinaison",
        indice_max: int,
        backend: str = _DEFAULT_BACKEND,
    ):
        """
        PrepareMat (packed): build the generator matrix for Gaussian elimination.

        Returns a GaussMatrix (PyIntGaussMatrix or PackedGaussMatrix depending
        on the backend parameter).
        """
        L          = C.L
        total_cols = indice_max * L
        L_mask     = (1 << L) - 1

        if backend == "packed":
            mat = PackedGaussMatrix(C.k_g, total_cols)
        else:
            mat = PyIntGaussMatrix(C.k_g, total_cols)

        ml = 0
        for j, comp in enumerate(C.components):
            gen    = C[j]
            k      = gen.k
            bc     = BitVect.zeros(k)
            bc.put_bit(0, 1)
            has_transform = len(comp.trans) > 0
            for _gl in range(k):
                state = gen.initialize_state(bc)
                bc  >>= 1
                # Pack first column group (mc=0)
                if has_transform:
                    val = comp.transform(state)._val & L_mask
                else:
                    val = state._val & L_mask
                row = val << (total_cols - L)
                # Pack remaining column groups
                for mc in range(1, indice_max):
                    state = next(gen)
                    if has_transform:
                        val = comp.transform(state)._val & L_mask
                    else:
                        val = state._val & L_mask
                    shift = total_cols - (mc + 1) * L
                    row |= val << shift

                if backend == "packed":
                    mat.set_row_from_int(ml, row)
                else:
                    mat.rows[ml] = row
                ml += 1
        return mat

    @staticmethod
    def _prepare_mat_words(C: "Combinaison", indice_max: int) -> list:
        """
        PrepareMat (word-level): build a list of k_g lists, each of length
        indice_max, holding raw generator state words (including 64-bit
        overflow from C's ulong arithmetic) for L <= 32.
        """
        kg  = C.k_g
        arr = [[0] * indice_max for _ in range(kg)]
        ml  = 0
        for j, comp in enumerate(C.components):
            gen = C[j]
            k   = gen.k
            bc  = BitVect.zeros(k)
            bc.put_bit(0, 1)
            for _gl in range(k):
                gen.initialize_state(bc)
                bc >>= 1
                arr[ml][0] = gen.gen_state._val & 0xFFFFFFFFFFFFFFFF
                for mc in range(1, indice_max):
                    next(gen)
                    arr[ml][mc] = gen.gen_state._val & 0xFFFFFFFFFFFFFFFF
                ml += 1
        return arr
