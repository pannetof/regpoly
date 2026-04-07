"""
generateur.py — Abstract base class for a single PRNG component.
"""

from __future__ import annotations

import sys
from abc import ABC, abstractmethod

from bitvect import BitVect


class Generateur(ABC):
    """
    Abstract base for a single PRNG component.

    Generateur is a Python iterator: call next(gen) to advance one step,
    or use it directly in a for-loop / itertools pipeline.

    Attributes
    ----------
    gen_state : BitVect  — k-bit state vector
    smax      : int      — max statistical dimension (sys.maxsize = INT_MAX sentinel)
    k         : int      — degree (total state bits)
    step      : int      — step size
    L         : int      — output resolution in bits
    """

    def __init__(self, k: int) -> None:
        self.gen_state: BitVect = BitVect.zeros(k)
        self.smax: int = sys.maxsize
        self.k: int = k
        self.step: int = 1
        self.L: int = 0

    @abstractmethod
    def display(self) -> None:
        """Print generator parameters to stdout."""

    @classmethod
    @abstractmethod
    def name(cls) -> str:
        """Return the generator type name string."""

    @abstractmethod
    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """Initialise state from init_bv; return the first output word."""

    @abstractmethod
    def __next__(self) -> BitVect:
        """Advance by one step and return the output word."""

    def __iter__(self) -> "Generateur":
        return self

    @abstractmethod
    def copy(self) -> "Generateur":
        """Return a deep copy of self."""

    def char_poly(self) -> BitVect:
        """
        polychar: compute the minimal polynomial via Berlekamp-Massey.

        Initialises a copy of this generator with a canonical state (public
        bit 0 = 1), collects 2*k output bits (MSB of each output word), and
        applies the Berlekamp-Massey algorithm (SequenceMinimalPolynomial).

        Returns a k-bit BitVect where public bit j holds the coefficient of
        x^(k-j) for j = 0..k-1 (the leading x^k term is not stored).
        """
        K = self.k
        gen = self.copy()

        init_bv = BitVect.zeros(K)
        init_bv.put_bit(0, 1)
        gen.initialize_state(init_bv)

        seq = [0] * (2 * K)
        for i in range(2 * K):
            out = next(gen)
            seq[i] = out.get_bit(0)

        # Berlekamp-Massey (SequenceMinimalPolynomial from polynomials.c)
        # C[j] = coefficient of x^j in the connection polynomial
        C = [0] * (K + 1)
        B = [0] * (K + 1)
        C[0] = B[0] = 1   # BVCanonic: public bit 0 = 1
        L = 0
        x = 1
        for N in range(2 * K):
            d = seq[N]
            for j in range(1, L + 1):
                if C[j]:
                    d ^= seq[N - j]
            d &= 1
            if d == 0:
                x += 1
            elif 2 * L > N:
                temp = [0] * (K + 1)
                for i in range(K + 1 - x):
                    temp[i + x] = B[i]
                for i in range(K + 1):
                    C[i] ^= temp[i]
                x += 1
            else:
                T = C[:]
                temp = [0] * (K + 1)
                for i in range(K + 1 - x):
                    temp[i + x] = B[i]
                for i in range(K + 1):
                    C[i] ^= temp[i]
                L = N + 1 - L
                B = T
                x = 1

        # BVPoly[k] = C[K-k] for k = 0..K-1  (bit K = C[0] is out of range)
        bv = BitVect.zeros(K)
        for k in range(K):
            bv.put_bit(k, C[K - k])
        return bv

    def transition_matrix(self) -> "Matrix":
        """
        TransitionMatrix: compute and return the K×K state-transition matrix.

        For each canonical basis vector eᵢ (i = 0..K-1), initialises a copy of
        this generator with eᵢ, advances one step, and reads the K output bits
        into column i of the matrix.

        Returns
        -------
        A : Matrix — K×K transition matrix (l=K, t=1)

        Raises
        ------
        ValueError if L < k (output resolution smaller than degree).
        """
        from matrix import Matrix

        K = self.k
        if self.L < K:
            raise ValueError(
                "Cannot compute transition matrix: "
                "resolution of the generator smaller than degree (L<k)"
            )

        gen = self.copy()
        bc  = BitVect.zeros(K)
        bc.put_bit(0, 1)          # BVCanonic(&BC, 0)
        A   = Matrix(K, K, 1)

        for i in range(K):
            gen.initialize_state(bc)   # INITGEN
            bc >>= 1                   # BVRShift(&BC, &BC, 1)
            state = next(gen)          # ITERATION

            for row in range(K):
                bv = A.get_bitvect(row, 0)
                bv.put_bit(i, state.get_bit(row))
                A.set_bitvect(row, 0, bv)

        return A

    @classmethod
    @abstractmethod
    def CreateListFromFile(cls, filename: str, L: int) -> list["Generateur"]:
        """Read generator parameters from filename and return them as a list.

        Each subclass parses its own file format and constructs its generators.
        L is the maximum output resolution in bits.  The caller is responsible
        for handling the 'same' flag and populating the Component.
        """
