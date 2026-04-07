"""
matrix.py — Matrix class implementing all matrix operations from vectorsF2.c.

Internal storage
----------------
A galois.GF(2) array of shape (nblignes, t * l), where:
    nblignes  : number of rows
    t         : number of column groups per row
    l         : bits per column group  (= BitVect width)

Column group j occupies flat columns [j*l : (j+1)*l].
Within each group, column 0 is the MSB (matching BitVect convention).

Forward-compatibility
---------------------
Every field operation goes through the module-level GF2 object.
To move to GF(q), replace:
    GF2 = galois.GF(2)
with:
    GF2 = galois.GF(q)
and adapt the BitVect ↔ array conversion helpers.
"""

import numpy as np
import galois

from bitvect import BitVect

GF2 = galois.GF(2)


class Matrix:
    __slots__ = ("nblignes", "l", "t", "_mat")

    # ------------------------------------------------------------------ #
    #  Construction                                                        #
    # ------------------------------------------------------------------ #

    def __init__(self, nblignes: int, l: int, t: int):
        """
        AllocMat: allocate a zero-filled matrix.

        nblignes : number of rows
        l        : bits per column group (BitVect width)
        t        : number of column groups per row
        """
        self.nblignes: int = nblignes
        self.l: int = l
        self.t: int = t
        self._mat: galois.FieldArray = GF2(np.zeros((nblignes, t * l), dtype=np.int64))

    # ------------------------------------------------------------------ #
    #  Internal helpers                                                   #
    # ------------------------------------------------------------------ #

    def _col_slice(self, j: int) -> slice:
        """Column slice in _mat for column group j."""
        return slice(j * self.l, (j + 1) * self.l)

    @staticmethod
    def _bv_to_np(bv: BitVect, n: int) -> np.ndarray:
        """Convert the first n bits of a BitVect to a 1-D int64 array (MSB first)."""
        return np.array([(bv._val >> (bv.n - 1 - i)) & 1 for i in range(n)],
                        dtype=np.int64)

    @staticmethod
    def _np_to_bv(arr, n: int) -> BitVect:
        """Convert a 1-D array of 0/1 values (MSB first) to a BitVect of width n."""
        val = 0
        a = np.asarray(arr, dtype=np.int64)
        for i in range(n):
            if a[i]:
                val |= 1 << (n - 1 - i)
        return BitVect(n, val)

    def _raw(self) -> np.ndarray:
        """Return a writable int64 view of the internal galois array."""
        return np.asarray(self._mat, dtype=np.int64)

    def _commit(self, raw: np.ndarray) -> None:
        """Write a modified int64 array back into _mat."""
        self._mat = GF2(raw)

    # ------------------------------------------------------------------ #
    #  Element access (BitVect ↔ matrix cell)                            #
    # ------------------------------------------------------------------ #

    def get_bitvect(self, row: int, j: int) -> BitVect:
        """Return a copy of the BitVect at (row, column group j)."""
        return self._np_to_bv(self._mat[row, self._col_slice(j)], self.l)

    def set_bitvect(self, row: int, j: int, bv: BitVect) -> None:
        """Set the BitVect at (row, column group j)."""
        self._mat[row, self._col_slice(j)] = GF2(self._bv_to_np(bv, self.l))

    def __getitem__(self, key) -> BitVect:
        """m[row, j] — shorthand for get_bitvect."""
        row, j = key
        return self.get_bitvect(row, j)

    def __setitem__(self, key, bv: BitVect) -> None:
        """m[row, j] = bv — shorthand for set_bitvect."""
        row, j = key
        self.set_bitvect(row, j, bv)

    # ------------------------------------------------------------------ #
    #  Copy                                                               #
    # ------------------------------------------------------------------ #

    def copy(self) -> "Matrix":
        """Return a full deep copy of this matrix."""
        m = Matrix(self.nblignes, self.l, self.t)
        m._mat = GF2(np.copy(self._raw()))
        return m

    def copy_from(self, src: "Matrix", nl: int, t: int) -> None:
        """
        CopyMat: copy nl rows and t column groups from src into self.
        Both matrices must have the same l.
        """
        if src.l != self.l:
            raise ValueError(f"CopyMat: l mismatch ({self.l} vs {src.l})")
        if nl > src.nblignes or t > src.t:
            raise ValueError("CopyMat: source too small")
        if nl > self.nblignes or t > self.t:
            raise ValueError("CopyMat: destination too small")
        self._mat[:nl, : t * self.l] = src._mat[:nl, : t * self.l]

    # ------------------------------------------------------------------ #
    #  Row operations                                                     #
    # ------------------------------------------------------------------ #

    def exchange_rows(self, i: int, j: int) -> None:
        """ExchangeVect: swap rows i and j in place."""
        if i != j:
            self._mat[[i, j]] = self._mat[[j, i]]

    def xor_rows(self, r: int, s: int, min_j: int = 0, max_j: int = None) -> None:
        """
        XorVect: row r ^= row s for column groups min_j .. max_j-1.
        Defaults to all column groups (0 .. t-1).
        """
        if max_j is None:
            max_j = self.t
        lo = min_j * self.l
        hi = max_j * self.l
        self._mat[r, lo:hi] += self._mat[s, lo:hi]

    # ------------------------------------------------------------------ #
    #  Gaussian elimination                                               #
    # ------------------------------------------------------------------ #

    def complete_elimination(self,
                              nblignes: int = None,
                              l: int = None,
                              t: int = None) -> int:
        """
        CompleteElimination: reduced row echelon form (RREF) over GF(2).
        Modifies self in place.  Returns the rank.
        """
        n  = nblignes if nblignes is not None else self.nblignes
        ll = l        if l        is not None else self.l
        tt = t        if t        is not None else self.t
        ncols = tt * ll
        rref = self._mat[:n, :ncols].row_reduce()
        self._mat[:n, :ncols] = rref
        return int(np.asarray(rref, dtype=np.int64).any(axis=1).sum())

    @property
    def rank(self) -> int:
        """Rank of the matrix over GF(2). Does not modify self."""
        return int(np.linalg.matrix_rank(self._mat))

    # ------------------------------------------------------------------ #
    #  Transpose                                                          #
    # ------------------------------------------------------------------ #

    def transpose(self) -> "Matrix":
        """
        TransposeMatrices: return the transpose as a new Matrix.

        For each column group s, the (nblignes × l) sub-matrix is transposed
        to (l × nblignes), producing a result with:
            nblignes_T = self.l
            l_T        = self.nblignes
            t_T        = self.t
        (The C source only handles t=1 and l ≤ 32; this implementation is
        general.)
        """
        T = Matrix(self.l, self.nblignes, self.t)
        raw_M = self._raw()
        raw_T = np.zeros((self.l, self.t * self.nblignes), dtype=np.int64)
        for s in range(self.t):
            M_s = raw_M[:, s * self.l : (s + 1) * self.l]            # (nblignes, l)
            raw_T[:, s * self.nblignes : (s + 1) * self.nblignes] = M_s.T
        T._mat = GF2(raw_T)
        return T

    # ------------------------------------------------------------------ #
    #  Inverse                                                            #
    # ------------------------------------------------------------------ #

    def inverse(self) -> tuple:
        """
        InverseMatrix: compute the GF(2) inverse of a square matrix.

        self must be square: t == 1 and nblignes == l.
        Does not modify self.

        Returns (True,  M_inv) if the matrix is invertible.
        Returns (False, None)  if the matrix is singular.
        """
        if self.t != 1 or self.nblignes != self.l:
            raise ValueError(
                "InverseMatrix: matrix must be square with t=1 and nblignes == l"
            )
        try:
            inv_raw = np.linalg.inv(self._mat)   # galois overrides np.linalg.inv
            M_inv = Matrix(self.nblignes, self.l, 1)
            M_inv._mat = GF2(np.asarray(inv_raw, dtype=np.int64))
            return True, M_inv
        except np.linalg.LinAlgError:
            return False, None

    # ------------------------------------------------------------------ #
    #  Display                                                            #
    # ------------------------------------------------------------------ #

    def display(self,
                t: int = None,
                l: int = None,
                kg: int = None) -> None:
        """
        DispMat: print the matrix to stdout.

        t, l, kg default to self.t, self.l, self.nblignes.
        """
        tt  = t  if t  is not None else self.t
        ll  = l  if l  is not None else self.l
        n   = kg if kg is not None else self.nblignes
        ncols = tt * ll

        header_100 = "    " + "".join(str((c // 100) % 10) for c in range(ncols))
        header_10  = "    " + "".join(str((c // 10)  % 10) for c in range(ncols))
        header_1   = "    " + "".join(str(c % 10)          for c in range(ncols))
        print(header_100)
        print(header_10)
        print(header_1)
        for i in range(n):
            row = "".join(str(int(self._mat[i, c])) for c in range(ncols))
            print(f"{i:3d}[{row}]")
        print()

    def __repr__(self) -> str:
        return (
            f"Matrix(nblignes={self.nblignes}, l={self.l}, t={self.t},"
            f" rank={self.rank})"
        )
