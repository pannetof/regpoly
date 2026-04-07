"""
generateur.py — Abstract base class for a single PRNG component.

Cannot be a @cython.cclass because it inherits from ABC.  The concrete
methods (char_poly, transition_matrix) still benefit from Cython compilation
of the module as a whole.
"""


import cython
import sys
from abc import ABC, abstractmethod

from regpoly.bitvect import BitVect


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
    L         : int      — output resolution in bits
    """

    def __init__(self, k: cython.int) -> None:
        self.gen_state: BitVect = BitVect.zeros(k)
        self.smax: cython.int = sys.maxsize
        self.k: cython.int = k
        self.L: cython.int = 0

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
        K: cython.int = self.k
        gen = self.copy()

        init_bv: BitVect = BitVect.zeros(K)
        init_bv.put_bit(0, 1)
        gen.initialize_state(init_bv)

        seq: list = [0] * (2 * K)
        i: cython.int
        for i in range(2 * K):
            out: BitVect = next(gen)
            seq[i] = out.get_bit(0)

        # Berlekamp-Massey (SequenceMinimalPolynomial from polynomials.c)
        C: list = [0] * (K + 1)
        B: list = [0] * (K + 1)
        C[0] = B[0] = 1
        L: cython.int = 0
        x: cython.int = 1
        N: cython.int
        j: cython.int
        d: cython.int
        for N in range(2 * K):
            d = seq[N]
            for j in range(1, L + 1):
                if C[j]:
                    d ^= seq[N - j]
            d &= 1
            if d == 0:
                x += 1
            elif 2 * L > N:
                temp: list = [0] * (K + 1)
                for i in range(K + 1 - x):
                    temp[i + x] = B[i]
                for i in range(K + 1):
                    C[i] ^= temp[i]
                x += 1
            else:
                T: list = C[:]
                temp = [0] * (K + 1)
                for i in range(K + 1 - x):
                    temp[i + x] = B[i]
                for i in range(K + 1):
                    C[i] ^= temp[i]
                L = N + 1 - L
                B = T
                x = 1

        bv: BitVect = BitVect.zeros(K)
        k: cython.int
        for k in range(K):
            bv.put_bit(k, C[K - k])
        return bv

    def _transition_state(self) -> object:
        """
        Return the k-bit state as an integer for transition matrix computation.

        Default: return gen_state._val directly (assumes gen_state.n == k).
        Override in subclasses where the state size differs from k
        (e.g., Tausworthe where gen_state.n == L, or MT with circular
        buffer rotation).
        """
        return self.gen_state._val

    def transition_matrix(self, transposed: bool = False) -> "Matrix":
        """
        TransitionMatrix: compute and return the K×K state-transition matrix.

        For each canonical basis vector eᵢ (i = 0..K-1), initialises a copy of
        this generator with eᵢ, advances one step, and reads the K state bits
        into column i of the matrix.

        The matrix is always K×K and does not depend on L.

        Parameters
        ----------
        transposed : bool — if True, return Aᵀ (row i = column i of A).

        Returns
        -------
        A : Matrix — K×K transition matrix (l=K, t=1), or its transpose.
        """
        from regpoly.matrix import Matrix

        K: object = self.k
        gen = self.copy()
        bc: BitVect = BitVect.zeros(K)
        bc.put_bit(0, 1)

        # Always build the transpose first: row i = column i of A
        # = the state after one step with basis vector eᵢ.
        # This is O(K) — one int assignment per column, no inner loop.
        At: Matrix = Matrix(K, K, 1)
        i: cython.int
        for i in range(K):
            gen.initialize_state(bc)
            bc >>= 1
            next(gen)
            At._rows[i] = gen._transition_state()

        if transposed:
            return At

        # Transpose Aᵀ → A: scatter each row of Aᵀ (= column of A)
        # into the rows of A.
        A: Matrix = Matrix(K, K, 1)
        col_bit: object
        row: cython.int
        for i in range(K):
            sv: object = At._rows[i]
            col_bit = 1 << (K - 1 - i)
            for row in range(K):
                if sv >> (K - 1 - row) & 1:
                    A._rows[row] |= col_bit
        return A

    @classmethod
    @abstractmethod
    def CreateListFromFile(cls, filename: str, L: int) -> list["Generateur"]:
        """Read generator parameters from filename and return them as a list."""
