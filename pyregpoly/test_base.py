"""
test_base.py — Abstract base class for ME/CF tests.

Provides the shared _prepare_mat() static method used by both
EquidistributionTest and CollisionFreeTest.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

from bitvect import BitVect
from matrix import Matrix
from test_results_base import AbstractTestResults

if TYPE_CHECKING:
    from combinaison import Combinaison


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
        PrepareMat: build the k_g × indice_max generator matrix.

        Row ml, column mc holds the L-bit tempered output at step mc when
        the generator is initialised with canonical basis vector ml (the
        unit vector with public bit ml set).

        Used by both EquidistributionTest and CollisionFreeTest.
        """
        L   = C.L
        mat = Matrix(C.k_g, L, indice_max)
        ml  = 0
        for j, comp in enumerate(C.components):
            gen = C[j]
            k   = gen.k
            bc  = BitVect.zeros(k)
            bc.put_bit(0, 1)        # canonic(0): public bit 0 (MSB) = 1
            for _gl in range(k):
                state    = gen.initialize_state(bc)
                bc     >>= 1        # next canonical basis vector
                mat[ml, 0] = comp.transform(state)
                for mc in range(1, indice_max):
                    state      = next(gen)
                    mat[ml, mc] = comp.transform(state)
                ml += 1
        return mat

    @staticmethod
    def _prepare_mat_words(C: "Combinaison", indice_max: int) -> np.ndarray:
        """
        PrepareMat (word-level): build a (k_g, indice_max) int64 array.

        Each entry stores the raw generator state word (including 64-bit
        overflow from C's ulong arithmetic) for L <= 32.  This replicates
        C's behaviour where Transform copies vect[0] (with overflow)
        directly into the matrix.
        """
        kg  = C.k_g
        arr = np.zeros((kg, indice_max), dtype=np.uint64)
        ml  = 0
        for j, comp in enumerate(C.components):
            gen = C[j]
            k   = gen.k
            bc  = BitVect.zeros(k)
            bc.put_bit(0, 1)
            for _gl in range(k):
                gen.initialize_state(bc)
                bc >>= 1
                # Store raw word (with overflow) — Transform just copies vect[0]
                # Mask to 64 bits to match C's ulong width
                arr[ml, 0] = gen.gen_state._val & 0xFFFFFFFFFFFFFFFF
                for mc in range(1, indice_max):
                    next(gen)
                    arr[ml, mc] = gen.gen_state._val & 0xFFFFFFFFFFFFFFFF
                ml += 1
        return arr
