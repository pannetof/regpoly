"""
matrix.py — Matrix class implementing all matrix operations from vectorsF2.c.

Internal storage
----------------
A list of `nblignes` Python ints, each representing a packed row of
`t * l` bits.  Bit 0 (MSB) is at position `(t*l - 1)` of the int.

Column group j occupies bits [j*l : (j+1)*l].
Within each group, column 0 is the MSB (matching BitVect convention).
"""

from regpoly.bitvect import BitVect


class Matrix:
    __slots__ = ("nblignes", "l", "t", "_rows")

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
        self._rows: list[int] = [0] * nblignes

    # ------------------------------------------------------------------ #
    #  Internal helpers                                                   #
    # ------------------------------------------------------------------ #

    @property
    def _total_cols(self) -> int:
        return self.t * self.l

    def _get_bit(self, row: int, col: int) -> int:
        """Get bit at (row, col). col 0 = MSB."""
        return (self._rows[row] >> (self._total_cols - 1 - col)) & 1

    def _put_bit(self, row: int, col: int, val: int) -> None:
        """Set bit at (row, col). col 0 = MSB."""
        pos = self._total_cols - 1 - col
        if val:
            self._rows[row] |= (1 << pos)
        else:
            self._rows[row] &= ~(1 << pos)

    # ------------------------------------------------------------------ #
    #  Element access (BitVect ↔ matrix cell)                            #
    # ------------------------------------------------------------------ #

    def get_bitvect(self, row: int, j: int) -> BitVect:
        """Return a copy of the BitVect at (row, column group j)."""
        tc = self._total_cols
        start = j * self.l
        # Extract l bits starting at position start
        val = (self._rows[row] >> (tc - start - self.l)) & ((1 << self.l) - 1)
        return BitVect(self.l, val)

    def set_bitvect(self, row: int, j: int, bv: BitVect) -> None:
        """Set the BitVect at (row, column group j)."""
        tc = self._total_cols
        start = j * self.l
        shift = tc - start - self.l
        mask = ((1 << self.l) - 1) << shift
        self._rows[row] = (self._rows[row] & ~mask) | ((bv._val & ((1 << self.l) - 1)) << shift)

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
        m._rows = self._rows[:]
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
        src_tc = src._total_cols
        dst_tc = self._total_cols
        nbits = t * self.l
        for i in range(nl):
            # Extract top nbits from src row, place into top nbits of dst row
            src_val = (src._rows[i] >> (src_tc - nbits)) & ((1 << nbits) - 1)
            mask = ((1 << nbits) - 1) << (dst_tc - nbits)
            self._rows[i] = (self._rows[i] & ~mask) | (src_val << (dst_tc - nbits))

    # ------------------------------------------------------------------ #
    #  Row operations                                                     #
    # ------------------------------------------------------------------ #

    def exchange_rows(self, i: int, j: int) -> None:
        """ExchangeVect: swap rows i and j in place."""
        if i != j:
            self._rows[i], self._rows[j] = self._rows[j], self._rows[i]

    def xor_rows(self, r: int, s: int, min_j: int = 0, max_j: int = None) -> None:
        """
        XorVect: row r ^= row s for column groups min_j .. max_j-1.
        Defaults to all column groups (0 .. t-1).
        """
        if max_j is None:
            self._rows[r] ^= self._rows[s]
        else:
            tc = self._total_cols
            lo = min_j * self.l
            hi = max_j * self.l
            nbits = hi - lo
            shift = tc - hi
            mask = ((1 << nbits) - 1) << shift
            self._rows[r] ^= self._rows[s] & mask

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
        n     = nblignes if nblignes is not None else self.nblignes
        ll    = l        if l        is not None else self.l
        tt    = t        if t        is not None else self.t
        ncols = tt * ll
        tc    = self._total_cols

        rang = 0
        for col in range(ncols):
            pos = tc - 1 - col
            # Find pivot
            i = rang
            while i < n and not (self._rows[i] >> pos & 1):
                i += 1
            if i >= n:
                continue
            if i != rang:
                self._rows[rang], self._rows[i] = self._rows[i], self._rows[rang]
            # Eliminate all other rows (both above and below for RREF)
            for ii in range(n):
                if ii != rang and (self._rows[ii] >> pos & 1):
                    self._rows[ii] ^= self._rows[rang]
            rang += 1
        return rang

    @property
    def rank(self) -> int:
        """Rank of the matrix over GF(2). Does not modify self."""
        # Work on a copy
        rows = self._rows[:]
        tc = self._total_cols
        n = self.nblignes
        rang = 0
        for col in range(tc):
            pos = tc - 1 - col
            i = rang
            while i < n and not (rows[i] >> pos & 1):
                i += 1
            if i >= n:
                continue
            if i != rang:
                rows[rang], rows[i] = rows[i], rows[rang]
            for ii in range(rang + 1, n):
                if rows[ii] >> pos & 1:
                    rows[ii] ^= rows[rang]
            rang += 1
        return rang

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
        """
        T = Matrix(self.l, self.nblignes, self.t)
        for s in range(self.t):
            for i in range(self.nblignes):
                for j in range(self.l):
                    bit = self._get_bit(i, s * self.l + j)
                    if bit:
                        T._put_bit(j, s * self.nblignes + i, 1)
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
        n = self.nblignes
        # Augmented matrix [A | I] using 2n-bit packed rows
        aug = [0] * n
        for i in range(n):
            aug[i] = (self._rows[i] << n) | (1 << (n - 1 - i))

        rang = 0
        for col in range(n):
            pos = 2 * n - 1 - col
            i = rang
            while i < n and not (aug[i] >> pos & 1):
                i += 1
            if i >= n:
                return False, None
            if i != rang:
                aug[rang], aug[i] = aug[i], aug[rang]
            for ii in range(n):
                if ii != rang and (aug[ii] >> pos & 1):
                    aug[ii] ^= aug[rang]
            rang += 1

        M_inv = Matrix(n, n, 1)
        mask_n = (1 << n) - 1
        for i in range(n):
            M_inv._rows[i] = aug[i] & mask_n
        return True, M_inv

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
        tt    = t  if t  is not None else self.t
        ll    = l  if l  is not None else self.l
        n     = kg if kg is not None else self.nblignes
        ncols = tt * ll

        header_100 = "    " + "".join(str((c // 100) % 10) for c in range(ncols))
        header_10  = "    " + "".join(str((c // 10)  % 10) for c in range(ncols))
        header_1   = "    " + "".join(str(c % 10)          for c in range(ncols))
        print(header_100)
        print(header_10)
        print(header_1)
        tc = self._total_cols
        for i in range(n):
            row = "".join(
                str((self._rows[i] >> (tc - 1 - c)) & 1)
                for c in range(ncols)
            )
            print(f"{i:3d}[{row}]")
        print()

    def display_transposed(self) -> None:
        """
        Display this matrix transposed: print column j of self as row j.

        Reads bit j from each row i to build output row j.
        Useful when self stores Aᵀ and we want to display A.
        """
        nrows = self.nblignes
        tc = self._total_cols

        header_100 = "    " + "".join(str((c // 100) % 10) for c in range(nrows))
        header_10  = "    " + "".join(str((c // 10)  % 10) for c in range(nrows))
        header_1   = "    " + "".join(str(c % 10)          for c in range(nrows))
        print(header_100)
        print(header_10)
        print(header_1)
        for j in range(tc):
            bit_pos = tc - 1 - j
            row = "".join(
                str((self._rows[i] >> bit_pos) & 1)
                for i in range(nrows)
            )
            print(f"{j:3d}[{row}]")
        print()

    def __repr__(self) -> str:
        return (
            f"Matrix(nblignes={self.nblignes}, l={self.l}, t={self.t},"
            f" rank={self.rank})"
        )
