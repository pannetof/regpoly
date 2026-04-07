"""
equidistribution_test.py — ME equidistribution test (METHOD_MATRICIAL).

METHOD_DUALLATTICE raises NotImplementedError if requested.
"""

from __future__ import annotations

import math
import sys
from typing import TYPE_CHECKING

import numpy as np

from test_base import AbstractTest
from equidistribution_results import EquidistributionResults

if TYPE_CHECKING:
    from combinaison import Combinaison

METHOD_MATRICIAL   = 0
METHOD_DUALLATTICE = 1
METHOD_NOTHING     = 2


class EquidistributionTest(AbstractTest):
    """
    ME test configuration and algorithm.

    Holds the test parameters; calling run(C) executes the test and
    returns an EquidistributionResults object.

    Attributes
    ----------
    L       : int       — maximum resolution to test
    delta   : list[int] — per-resolution gap bound (delta[l]); delta[0] unused
    mse     : int       — maximum allowed sum of gaps (quasi-ME threshold)
    meverif : bool      — whether the test is enabled
    method  : int       — METHOD_MATRICIAL or METHOD_NOTHING
    """

    def __init__(
        self,
        L: int,
        delta: list[int],
        mse: int,
        meverif: bool = True,
        method: int = METHOD_MATRICIAL,
    ) -> None:
        if method == METHOD_DUALLATTICE:
            raise NotImplementedError("METHOD_DUALLATTICE is not implemented")
        if method not in (METHOD_MATRICIAL, METHOD_NOTHING):
            raise ValueError(f"Unknown method: {method}")
        self.L       = L
        self.delta   = list(delta)      # indexed 0..L; delta[0] unused
        self.mse     = mse
        self.meverif = meverif
        self.method  = method

    # -- AbstractTest interface -------------------------------------------

    def run(self, C: "Combinaison", *args, **kwargs) -> EquidistributionResults:
        """
        TestME (METHOD_MATRICIAL): compute dimension gaps Delta_l for all l
        in Psi_12.

        Builds the k_g × smax generator matrix, then for each resolution l
        in Psi_12 runs Gaussian elimination to determine t_l (the actual
        equidistribution dimension at resolution l).  Short-circuits when a
        gap exceeds delta[l] or the running sum of gaps exceeds mse.

        Returns an EquidistributionResults object with the computed gaps.
        """
        if not self.meverif or self.method == METHOD_NOTHING:
            return EquidistributionResults(
                L=self.L, ecart=[0] * (self.L + 1),
                psi12=[False] * (self.L + 1), se=0,
                verified=False, mse=self.mse,
                meverif=self.meverif, delta=self.delta,
            )

        if self.L < C.L:
            raise ValueError(
                f"TestME: EquidistributionTest.L ({self.L}) < generator L ({C.L})"
            )

        ecart = [-1] * (self.L + 1)
        se    = 0
        psi12 = self._compute_psi12(C)

        indice_max = min(C.k_g, C.smax)
        mat        = self._prepare_mat(C, indice_max)
        raw_full   = np.asarray(mat._mat, dtype=np.int64)  # (k_g, indice_max * L)

        verif = False
        maxl  = self.L
        l     = 1
        while l <= self.L:
            if ecart[l] == -1 and (psi12[l] or verif):
                t   = min(C.k_g // l, C.smax)
                raw = raw_full[:C.k_g, :t * C.L].copy()
                t_l = self._dimension_equid(raw, C.k_g, l, C.L, C.smax)
                ecart[l] = t - t_l
                se       += ecart[l]

                if ecart[l] > self.delta[l] or se > self.mse:
                    maxl = l
                    break

                if ecart[l] != 0:
                    verif = True
                    if l != 1:
                        l -= 2      # step back; outer l += 1 will give l-1
                else:
                    verif = False
            l += 1

        # Finalise: uncomputed entries within maxl are 0; beyond are maxsize.
        se = 0
        for l in range(1, maxl + 1):
            if ecart[l] == -1:
                ecart[l] = 0
            se += ecart[l]
        for l in range(maxl + 1, self.L + 1):
            if ecart[l] == -1:
                ecart[l] = sys.maxsize

        return EquidistributionResults(
            L=self.L, ecart=ecart, psi12=psi12, se=se,
            verified=True, mse=self.mse,
            meverif=self.meverif, delta=self.delta,
        )

    # -- Private helpers --------------------------------------------------

    def _compute_psi12(self, C: "Combinaison") -> list[bool]:
        """SetPsi12: compute Psi_12 — resolutions l whose gap must be tested."""
        psi = [False] * (self.L + 1)
        # PSI 1: l = 1 .. floor(sqrt(k_g))
        r = int(math.isqrt(C.k_g))
        for l in range(1, min(r, self.L) + 1):
            psi[l] = True
        # PSI 2: floor(k_g / t) for t = m .. floor(sqrt(k_g - 1))
        r2 = int(math.isqrt(C.k_g - 1))
        m  = C.k_g // self.L
        if m < 2:
            m = 2
        for t in range(m, r2 + 1):
            psi[min(C.k_g // t, self.L)] = True
        return psi

    @staticmethod
    def _dimension_equid(
        raw: np.ndarray,
        kg: int,
        l: int,
        L: int,
        smax: int,
    ) -> int:
        """
        DimensionEquid: Gaussian elimination to compute t_l.

        raw  : (kg, t*L) int64 array — modified in place; caller must pass a copy
        L    : bits per column group (= generator output resolution)
        l    : resolution being tested (bits per group used for pivoting)
        Returns t_l — number of complete l-bit columns eliminated (= actual
        equidistribution dimension at resolution l).
        """
        t    = min(kg // l, smax)
        rang = 0
        for j in range(t):
            for cl in range(l):
                col = j * L + cl        # bit cl (MSB-first) of column group j
                i   = rang
                while i < kg and raw[i, col] == 0:
                    i += 1
                if i < kg:              # pivot found
                    if i != rang:
                        raw[[rang, i]] = raw[[i, rang]]
                    lo = j * L
                    hi = t * L
                    for ii in range(rang + 1, kg):
                        if raw[ii, col]:
                            raw[ii, lo:hi] ^= raw[rang, lo:hi]
                    rang += 1
                else:
                    return j            # no pivot: j complete columns processed
        return t

    @staticmethod
    def _resolution_equid(
        raw: np.ndarray,
        kg: int,
        t: int,
        L: int,
        indices: list,
    ) -> int:
        """
        ResolutionEquid: Gaussian elimination to compute l_t.

        Outer loop iterates over bit levels cl, inner over the t selected
        column groups.  XOR operations span all columns in raw (0 to raw.shape[1]).

        raw     : (kg, ncols) int64 array — modified in place; share across calls
                  for the same matrix (as in the C code) to accumulate reduction.
        kg      : number of rows
        t       : dimension (number of column groups to test)
        L       : bits per column group in raw
        indices : list of t column-group indices to use

        Returns l_t — number of complete bit levels eliminated (= actual
        equidistribution resolution for dimension t).
        """
        l    = min(kg // t, L)
        rang = 0
        for cl in range(l):
            for j in range(t):
                col = indices[j] * L + cl  # bit cl of column group indices[j]
                i   = rang
                while i < kg and raw[i, col] == 0:
                    i += 1
                if i < kg:              # pivot found
                    if i != rang:
                        raw[[rang, i]] = raw[[i, rang]]
                    for ii in range(rang + 1, kg):
                        if raw[ii, col]:
                            raw[ii, :] ^= raw[rang, :]   # XOR all columns
                    rang += 1
                else:
                    return cl           # no pivot at this bit level
        return l

    @staticmethod
    def _resolution_equid_words(
        words: np.ndarray,
        kg: int,
        t: int,
        L: int,
        indices: list,
    ) -> int:
        """
        ResolutionEquid (word-level): Gaussian elimination on a (kg, ncols)
        int64 array where each entry is a raw 32-bit word (potentially with
        64-bit overflow from C's ulong arithmetic).

        Replicates C's pivot search: ``vect[cldiv] < MMC[clmod]`` where
        MMC[clmod] = 0x80000000 >> clmod.  The XOR spans all columns.
        """
        l    = min(kg // t, L)
        rang = 0
        for cl in range(l):
            mmc = 0x80000000 >> cl
            for j in range(t):
                col = indices[j]
                i   = rang
                while i < kg and words[i, col] < mmc:
                    i += 1
                if i < kg:
                    if i != rang:
                        words[[rang, i]] = words[[i, rang]]
                    for ii in range(rang + 1, kg):
                        if words[ii, col] & mmc:
                            words[ii, :] ^= words[rang, :]
                    rang += 1
                else:
                    return cl
        return l
