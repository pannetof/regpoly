"""
gauss_matrix.py — Matrix backends for Gaussian elimination.

Two implementations with identical interfaces:
  PyIntGaussMatrix  — rows stored as Python arbitrary-precision ints (current)
  PackedGaussMatrix — rows stored as arrays of uint64 words (native C speed)

Both provide row_xor, bit_test, first_nonzero_in_col, swap_rows, and
mask_xor for the Gaussian elimination routines.
"""

import cython
import array as _array


# ── Python int backend (original) ──────────────────────────────────────

class PyIntGaussMatrix:
    """
    Gaussian elimination matrix stored as list[int].

    Each row is a single Python arbitrary-precision int of total_cols bits,
    with bit 0 (MSB) at position total_cols-1 of the int.
    """
    __slots__ = ("nrows", "total_cols", "rows")

    def __init__(self, nrows: int, total_cols: int) -> None:
        self.nrows: int = nrows
        self.total_cols: int = total_cols
        self.rows: list = [0] * nrows

    def copy(self) -> "PyIntGaussMatrix":
        m = PyIntGaussMatrix.__new__(PyIntGaussMatrix)
        m.nrows = self.nrows
        m.total_cols = self.total_cols
        m.rows = self.rows[:]
        return m

    def bit_test(self, row: int, col: int) -> int:
        """Return 1 if bit at (row, col) is set, 0 otherwise. col 0 = MSB."""
        return (self.rows[row] >> (self.total_cols - 1 - col)) & 1

    def swap_rows(self, a: int, b: int) -> None:
        self.rows[a], self.rows[b] = self.rows[b], self.rows[a]

    def row_xor(self, dst: int, src: int) -> None:
        """rows[dst] ^= rows[src]"""
        self.rows[dst] ^= self.rows[src]

    def masked_row_xor(self, dst: int, src: int, col_lo: int, col_hi: int) -> None:
        """rows[dst] ^= (rows[src] masked to columns [col_lo, col_hi))."""
        nbits: int = col_hi - col_lo
        shift: int = self.total_cols - col_hi
        mask: object = ((1 << nbits) - 1) << shift
        self.rows[dst] ^= self.rows[src] & mask

    def find_pivot(self, col: int, start_row: int) -> int:
        """Return first row >= start_row with bit set at col, or -1."""
        bit_pos: object = self.total_cols - 1 - col
        i: int = start_row
        while i < self.nrows:
            if (self.rows[i] >> bit_pos) & 1:
                return i
            i += 1
        return -1

    def eliminate_column(self, pivot_row: int, col: int, start_row: int) -> None:
        """XOR pivot_row into all rows [start_row..nrows) that have col set."""
        bit_pos: object = self.total_cols - 1 - col
        pivot: object = self.rows[pivot_row]
        i: int
        for i in range(start_row, self.nrows):
            if (self.rows[i] >> bit_pos) & 1:
                self.rows[i] ^= pivot

    def eliminate_column_masked(
        self, pivot_row: int, col: int, start_row: int,
        col_lo: int, col_hi: int,
    ) -> None:
        """XOR masked pivot_row into rows [start_row..nrows) that have col set."""
        bit_pos: object = self.total_cols - 1 - col
        nbits: int = col_hi - col_lo
        shift: int = self.total_cols - col_hi
        mask: object = ((1 << nbits) - 1) << shift
        pivot_masked: object = self.rows[pivot_row] & mask
        i: int
        for i in range(start_row, self.nrows):
            if (self.rows[i] >> bit_pos) & 1:
                self.rows[i] ^= pivot_masked


# ── Packed uint64 backend ──────────────────────────────────────────────

_WL: cython.int = 64


class PackedGaussMatrix:
    """
    Gaussian elimination matrix stored as flat array of uint64 words.

    Each row occupies nwords consecutive uint64 entries in a flat array.
    Bit 0 (MSB of the logical row) is at the MSB of word 0.
    """
    __slots__ = ("nrows", "total_cols", "nwords", "data")

    def __init__(self, nrows: int, total_cols: int) -> None:
        self.nrows: cython.int = nrows
        self.total_cols: cython.int = total_cols
        nw: cython.int = (total_cols + _WL - 1) // _WL
        self.nwords: cython.int = nw
        self.data: list = [0] * (nrows * nw)

    def copy(self) -> "PackedGaussMatrix":
        m = PackedGaussMatrix.__new__(PackedGaussMatrix)
        m.nrows = self.nrows
        m.total_cols = self.total_cols
        m.nwords = self.nwords
        m.data = self.data[:]
        return m

    def set_row_from_int(self, row: cython.int, val: object) -> None:
        """Set row from a Python int (MSB = bit 0). Used during matrix build."""
        base: cython.int = row * self.nwords
        nw: cython.int = self.nwords
        w: cython.int
        for w in range(nw):
            shift: object = (nw - 1 - w) * _WL
            self.data[base + w] = (val >> shift) & 0xFFFFFFFFFFFFFFFF

    def bit_test(self, row: cython.int, col: cython.int) -> cython.int:
        """Return 1 if bit at (row, col) is set. col 0 = MSB."""
        base: cython.int = row * self.nwords
        w: cython.int = col // _WL
        b: cython.int = _WL - 1 - (col % _WL)
        return (self.data[base + w] >> b) & 1

    def swap_rows(self, a: cython.int, b: cython.int) -> None:
        base_a: cython.int = a * self.nwords
        base_b: cython.int = b * self.nwords
        nw: cython.int = self.nwords
        w: cython.int
        for w in range(nw):
            self.data[base_a + w], self.data[base_b + w] = (
                self.data[base_b + w], self.data[base_a + w]
            )

    def row_xor(self, dst: cython.int, src: cython.int) -> None:
        """data[dst] ^= data[src], word by word."""
        base_d: cython.int = dst * self.nwords
        base_s: cython.int = src * self.nwords
        nw: cython.int = self.nwords
        w: cython.int
        for w in range(nw):
            self.data[base_d + w] ^= self.data[base_s + w]

    def masked_row_xor(
        self, dst: cython.int, src: cython.int,
        col_lo: cython.int, col_hi: cython.int,
    ) -> None:
        """data[dst] ^= data[src] for bits in columns [col_lo, col_hi)."""
        base_d: cython.int = dst * self.nwords
        base_s: cython.int = src * self.nwords
        w_lo: cython.int = col_lo // _WL
        w_hi: cython.int = (col_hi - 1) // _WL  # inclusive last word
        w: cython.int
        # Full words in the middle
        for w in range(w_lo, w_hi + 1):
            mask: object = 0xFFFFFFFFFFFFFFFF
            # Mask out bits before col_lo in the first word
            if w == w_lo:
                bit_in_word: cython.int = col_lo % _WL
                if bit_in_word != 0:
                    mask &= (1 << (_WL - bit_in_word)) - 1
            # Mask out bits after col_hi in the last word
            if w == w_hi:
                tail: cython.int = (w_hi + 1) * _WL - col_hi
                if tail != 0:
                    mask &= ~((1 << tail) - 1)
            self.data[base_d + w] ^= self.data[base_s + w] & mask

    def find_pivot(self, col: cython.int, start_row: cython.int) -> cython.int:
        """Return first row >= start_row with bit set at col, or -1."""
        w: cython.int = col // _WL
        b_mask: object = 1 << (_WL - 1 - (col % _WL))
        nw: cython.int = self.nwords
        i: cython.int
        for i in range(start_row, self.nrows):
            if self.data[i * nw + w] & b_mask:
                return i
        return -1

    def eliminate_column(self, pivot_row: cython.int, col: cython.int, start_row: cython.int) -> None:
        """XOR pivot_row into all rows [start_row..nrows) that have col set."""
        nw: cython.int = self.nwords
        base_p: cython.int = pivot_row * nw
        w_col: cython.int = col // _WL
        b_mask: object = 1 << (_WL - 1 - (col % _WL))
        i: cython.int
        w: cython.int
        for i in range(start_row, self.nrows):
            base_i: cython.int = i * nw
            if self.data[base_i + w_col] & b_mask:
                for w in range(nw):
                    self.data[base_i + w] ^= self.data[base_p + w]

    def eliminate_column_masked(
        self, pivot_row: cython.int, col: cython.int, start_row: cython.int,
        col_lo: cython.int, col_hi: cython.int,
    ) -> None:
        """XOR masked pivot_row into rows [start_row..nrows) that have col set."""
        nw: cython.int = self.nwords
        base_p: cython.int = pivot_row * nw
        w_col: cython.int = col // _WL
        b_mask: object = 1 << (_WL - 1 - (col % _WL))
        w_lo: cython.int = col_lo // _WL
        w_hi: cython.int = (col_hi - 1) // _WL

        # Precompute word masks
        masks: list = []
        w: cython.int
        for w in range(w_lo, w_hi + 1):
            m: object = 0xFFFFFFFFFFFFFFFF
            if w == w_lo:
                bit_in_word: cython.int = col_lo % _WL
                if bit_in_word != 0:
                    m &= (1 << (_WL - bit_in_word)) - 1
            if w == w_hi:
                tail: cython.int = (w_hi + 1) * _WL - col_hi
                if tail != 0:
                    m &= ~((1 << tail) - 1)
            masks.append(m)

        # Precompute masked pivot words
        pwords: list = []
        for idx in range(len(masks)):
            pwords.append(self.data[base_p + w_lo + idx] & masks[idx])

        i: cython.int
        for i in range(start_row, self.nrows):
            base_i: cython.int = i * nw
            if self.data[base_i + w_col] & b_mask:
                for idx in range(len(masks)):
                    self.data[base_i + w_lo + idx] ^= pwords[idx]
