"""
test_matrix.py — Complete test battery for Matrix, parametrised by size.

The `size` fixture drives both the bit-vector width (l) and the number of
rows (nblignes) so that square matrices of any requested dimension are tested.

Run with default size (8):
    pytest test_matrix.py -v

Run with multiple sizes:
    pytest test_matrix.py -v --sizes 8,16,32,64,94
"""

import re
import numpy as np
import pytest
from bitvect import BitVect
from matrix import Matrix, GF2


# ------------------------------------------------------------------ #
#  Helpers                                                            #
# ------------------------------------------------------------------ #

def identity_matrix(n):
    """Return an n×n identity Matrix (t=1, l=n)."""
    m = Matrix(n, n, 1)
    for i in range(n):
        m[i, 0] = BitVect(n, 1 << (n - 1 - i))
    return m


def alt_val(n):
    """n-bit alternating pattern MSB=1."""
    v = 0
    for i in range(0, n, 2):
        v |= 1 << (n - 1 - i)
    return v


def all_ones(n):
    return (1 << n) - 1


# ================================================================== #
#  Construction                                                        #
# ================================================================== #

class TestConstruction:
    def test_alloc_shape(self, size):
        m = Matrix(size, size, 2)
        assert m.nblignes == size
        assert m.l == size
        assert m.t == 2
        assert np.asarray(m._mat).shape == (size, 2 * size)

    def test_alloc_zero_filled(self, size):
        m = Matrix(size, size, 1)
        assert np.all(np.asarray(m._mat) == 0)

    def test_alloc_single_row(self, size):
        m = Matrix(1, size, 1)
        assert np.asarray(m._mat).shape == (1, size)


# ================================================================== #
#  Element access                                                     #
# ================================================================== #

class TestElementAccess:
    def test_set_get_bitvect_group0(self, size):
        m = Matrix(size, size, 2)
        bv = BitVect(size, alt_val(size))
        m.set_bitvect(0, 0, bv)
        result = m.get_bitvect(0, 0)
        assert result._val == alt_val(size)
        assert result.n == size

    def test_set_get_bitvect_group1(self, size):
        m = Matrix(size, size, 2)
        bv = BitVect(size, all_ones(size))
        m.set_bitvect(size - 1, 1, bv)
        assert m.get_bitvect(size - 1, 1)._val == all_ones(size)
        assert m.get_bitvect(size - 1, 0)._val == 0

    def test_getitem_setitem(self, size):
        m = Matrix(size, size, 1)
        bv = BitVect(size, alt_val(size))
        m[0, 0] = bv
        assert m[0, 0]._val == alt_val(size)

    def test_get_returns_copy(self, size):
        """get_bitvect returns a copy, not a live reference."""
        m = Matrix(size, size, 1)
        m[0, 0] = BitVect(size, alt_val(size))
        bv = m[0, 0]
        bv._val = 0
        assert m[0, 0]._val == alt_val(size)


# ================================================================== #
#  Copy                                                               #
# ================================================================== #

class TestCopy:
    def test_deep_copy(self, size):
        m = identity_matrix(size)
        m2 = m.copy()
        raw_m = np.asarray(m._mat, dtype=np.int64)
        raw_m2 = np.asarray(m2._mat, dtype=np.int64)
        assert np.array_equal(raw_m, raw_m2)
        m2[0, 0] = BitVect(size, 0)
        assert m[0, 0]._val == 1 << (size - 1)  # original unchanged

    def test_copy_from_partial(self, size):
        half = size // 2
        src = Matrix(size, size, 2)
        src[0, 0] = BitVect(size, all_ones(size))
        src[1, 1] = BitVect(size, alt_val(size))
        dst = Matrix(size + 4, size, 2)
        dst.copy_from(src, nl=2, t=2)
        assert dst[0, 0]._val == all_ones(size)
        assert dst[1, 1]._val == alt_val(size)

    def test_copy_from_l_mismatch(self, size):
        src = Matrix(2, size, 1)
        dst = Matrix(2, size * 2, 1)
        with pytest.raises(ValueError, match="l mismatch"):
            dst.copy_from(src, 2, 1)

    def test_copy_from_source_too_small(self, size):
        src = Matrix(2, size, 1)
        dst = Matrix(4, size, 1)
        with pytest.raises(ValueError, match="source too small"):
            dst.copy_from(src, 3, 1)

    def test_copy_from_dest_too_small(self, size):
        src = Matrix(size, size, 2)
        dst = Matrix(2, size, 1)
        with pytest.raises(ValueError, match="destination too small"):
            dst.copy_from(src, size, 2)


# ================================================================== #
#  Row operations                                                     #
# ================================================================== #

class TestRowOperations:
    def test_exchange_rows(self, size):
        m = Matrix(size, size, 1)
        m[0, 0] = BitVect(size, all_ones(size))
        m[size - 1, 0] = BitVect(size, 0)
        m.exchange_rows(0, size - 1)
        assert m[0, 0]._val == 0
        assert m[size - 1, 0]._val == all_ones(size)

    def test_exchange_same_row(self, size):
        m = Matrix(size, size, 1)
        m[0, 0] = BitVect(size, alt_val(size))
        m.exchange_rows(0, 0)
        assert m[0, 0]._val == alt_val(size)

    def test_xor_rows_full(self, size):
        m = Matrix(size, size, 1)
        a = alt_val(size)
        b = all_ones(size)
        m[0, 0] = BitVect(size, a)
        m[1, 0] = BitVect(size, b)
        m.xor_rows(0, 1)
        assert m[0, 0]._val == a ^ b
        assert m[1, 0]._val == b  # source unchanged

    def test_xor_rows_partial(self, size):
        m = Matrix(size, size, 2)
        m[0, 0] = BitVect(size, all_ones(size))
        m[0, 1] = BitVect(size, all_ones(size))
        m[1, 0] = BitVect(size, alt_val(size))
        m[1, 1] = BitVect(size, alt_val(size))
        # XOR only column group 1
        m.xor_rows(0, 1, min_j=1, max_j=2)
        assert m[0, 0]._val == all_ones(size)   # group 0 unchanged
        assert m[0, 1]._val == all_ones(size) ^ alt_val(size)


# ================================================================== #
#  Gaussian elimination                                               #
# ================================================================== #

class TestCompleteElimination:
    def test_identity_full_rank(self, size):
        m = identity_matrix(size)
        assert m.complete_elimination() == size

    def test_duplicate_row_rank(self, size):
        m = Matrix(size, size, 1)
        m[0, 0] = BitVect(size, 1 << (size - 1))
        m[1, 0] = BitVect(size, 1)
        m[2, 0] = BitVect(size, 1 << (size - 1))  # duplicate of row 0
        assert m.complete_elimination() == 2

    def test_zero_matrix_rank_zero(self, size):
        assert Matrix(size, size, 1).complete_elimination() == 0

    def test_result_in_rref(self, size):
        m = identity_matrix(size)
        m.complete_elimination()
        raw = np.asarray(m._mat, dtype=np.int64)
        assert np.array_equal(raw, np.eye(size, dtype=np.int64))


# ================================================================== #
#  Rank property                                                      #
# ================================================================== #

class TestRank:
    def test_identity_rank(self, size):
        assert identity_matrix(size).rank == size

    def test_zero_matrix_rank(self, size):
        assert Matrix(size, size, 1).rank == 0

    def test_rank_does_not_modify(self, size):
        m = identity_matrix(size)
        raw_before = np.copy(np.asarray(m._mat))
        _ = m.rank
        assert np.array_equal(np.asarray(m._mat), raw_before)

    def test_rank_with_duplicate_rows(self, size):
        m = Matrix(size, size, 1)
        for i in range(size):
            m[i, 0] = BitVect(size, 1 << (size - 1))  # all rows identical
        assert m.rank == 1


# ================================================================== #
#  Diag                                                               #
# ================================================================== #

# ================================================================== #
#  Transpose                                                          #
# ================================================================== #

class TestTranspose:
    def test_shape_square(self, size):
        m = identity_matrix(size)
        T = m.transpose()
        assert T.nblignes == size
        assert T.l == size
        assert T.t == 1

    def test_identity_is_self_transpose(self, size):
        m = identity_matrix(size)
        T = m.transpose()
        assert np.array_equal(
            np.asarray(m._mat, dtype=np.int64),
            np.asarray(T._mat, dtype=np.int64),
        )

    def test_double_transpose(self, size):
        m = Matrix(size, size, 1)
        m[0, 0] = BitVect(size, alt_val(size))
        m[size // 2, 0] = BitVect(size, all_ones(size))
        TT = m.transpose().transpose()
        assert np.array_equal(
            np.asarray(m._mat, dtype=np.int64),
            np.asarray(TT._mat, dtype=np.int64),
        )

    def test_values_non_square(self, size):
        """Transpose of a 2-row matrix has 2 columns."""
        half = max(1, size // 2)
        m = Matrix(2, half, 1)
        m[0, 0] = BitVect(half, all_ones(half))   # row 0: all ones
        m[1, 0] = BitVect(half, 0)                # row 1: all zeros
        T = m.transpose()
        # T has shape (half, 2, 1)
        assert T.nblignes == half and T.l == 2 and T.t == 1
        # Each column of m becomes a row of T.
        # Column 0 of m: [1, 0] → T row 0 = [1, 0] → _val with n=2: 0b10 = 2
        # Column k (for k < half): m[0,0][k]=1, m[1,0][k]=0 → T[k,0]._val = 0b10
        for k in range(half):
            assert T[k, 0]._val == 0b10  # [1, 0] in n=2 MSB-first

    def test_multi_group_shape(self, size):
        m = Matrix(size, size, 2)
        T = m.transpose()
        assert T.nblignes == size
        assert T.l == size
        assert T.t == 2


# ================================================================== #
#  Inverse                                                            #
# ================================================================== #

class TestInverse:
    def test_identity_inverse_is_identity(self, size):
        m = identity_matrix(size)
        ok, inv = m.inverse()
        assert ok is True
        assert np.array_equal(
            np.asarray(inv._mat, dtype=np.int64),
            np.eye(size, dtype=np.int64),
        )

    def test_invertible_product_is_identity(self, size):
        """M @ M_inv = I for a full-rank matrix."""
        m = identity_matrix(size)
        # Perturb: XOR row 0 with row 1 to get a non-trivial invertible matrix
        m.xor_rows(0, 1)
        ok, inv = m.inverse()
        assert ok is True
        raw_m = np.asarray(m._mat, dtype=np.int64)
        raw_inv = np.asarray(inv._mat, dtype=np.int64)
        prod = (raw_m @ raw_inv) % 2
        assert np.array_equal(prod, np.eye(size, dtype=np.int64))

    def test_singular_returns_false(self, size):
        m = Matrix(size, size, 1)
        # All rows identical → singular
        for i in range(size):
            m[i, 0] = BitVect(size, alt_val(size))
        ok, inv = m.inverse()
        assert ok is False and inv is None

    def test_does_not_modify_original(self, size):
        m = identity_matrix(size)
        raw_before = np.copy(np.asarray(m._mat))
        m.inverse()
        assert np.array_equal(np.asarray(m._mat), raw_before)

    def test_non_square_raises(self, size):
        m = Matrix(size, size, 2)
        with pytest.raises(ValueError, match="square"):
            m.inverse()

    def test_non_square_rows_raises(self, size):
        m = Matrix(size + 1, size, 1)
        with pytest.raises(ValueError, match="square"):
            m.inverse()


# ================================================================== #
#  Display                                                            #
# ================================================================== #

class TestDisplay:
    def test_display_row_count(self, size, capsys):
        m = identity_matrix(size)
        m.display()
        out = capsys.readouterr().out
        data_lines = re.findall(r"^\s*\d+\[.*\]", out, re.MULTILINE)
        assert len(data_lines) == size

    def test_display_partial_row_count(self, size, capsys):
        half = size // 2
        m = identity_matrix(size)
        m.display(kg=half)
        out = capsys.readouterr().out
        data_lines = re.findall(r"^\s*\d+\[.*\]", out, re.MULTILINE)
        assert len(data_lines) == half

    def test_display_col_width(self, size, capsys):
        m = identity_matrix(size)
        m.display()
        out = capsys.readouterr().out
        data_lines = re.findall(r"^\s*\d+\[(.*?)\]", out, re.MULTILINE)
        assert all(len(l) == size for l in data_lines)

    def test_repr(self, size):
        m = identity_matrix(size)
        r = repr(m)
        assert f"nblignes={size}" in r
        assert f"l={size}" in r
        assert "t=1" in r


# ================================================================== #
#  Internal helpers                                                   #
# ================================================================== #

class TestHelpers:
    def test_bv_to_np_roundtrip(self, size):
        bv = BitVect(size, alt_val(size))
        arr = Matrix._bv_to_np(bv, size)
        assert arr.shape == (size,)
        assert arr.dtype == np.int64
        bv2 = Matrix._np_to_bv(arr, size)
        assert bv2._val == bv._val and bv2.n == bv.n

    def test_bv_to_np_msb_first(self, size):
        bv = BitVect(size, 1 << (size - 1))  # only MSB
        arr = Matrix._bv_to_np(bv, size)
        assert arr[0] == 1
        assert all(arr[i] == 0 for i in range(1, size))

    def test_np_to_bv_msb_first(self, size):
        arr = np.zeros(size, dtype=np.int64)
        arr[0] = 1  # MSB
        bv = Matrix._np_to_bv(arr, size)
        assert bv._val == 1 << (size - 1)

    def test_raw_commit_roundtrip(self, size):
        m = identity_matrix(size)
        raw = m._raw()
        assert raw.dtype == np.int64
        raw[0, :] = 0  # zero out first row
        m._commit(raw)
        assert m[0, 0]._val == 0

    def test_col_slice(self, size):
        t = 3
        m = Matrix(2, size, t)
        for j in range(t):
            s = m._col_slice(j)
            assert s == slice(j * size, (j + 1) * size)


# ================================================================== #
#  Integration                                                        #
# ================================================================== #

class TestIntegration:
    def test_complete_elimination_rank_matches_rank_property(self, size):
        m = Matrix(size, size, 1)
        for i in range(size):
            m[i, 0] = BitVect(size, 1 << (size - 1 - i))
        m.xor_rows(0, size // 2)
        expected = m.rank
        m2 = m.copy()
        assert m2.complete_elimination() == expected

    def test_symmetric_matrix_transpose_equals_self(self, size):
        """Symmetric M: M == M^T."""
        m = identity_matrix(size)  # identity is symmetric
        T = m.transpose()
        assert np.array_equal(
            np.asarray(m._mat, dtype=np.int64),
            np.asarray(T._mat, dtype=np.int64),
        )

    def test_rank_after_row_operations(self, size):
        """XOR-ing a row with itself (via xor_rows) drops rank by 1."""
        m = identity_matrix(size)
        before = m.rank
        # XOR row 0 with row 1: row 0 becomes e_0 ^ e_1, still full rank
        m.xor_rows(0, 1)
        assert m.rank == before  # XOR of distinct standard basis vecs keeps rank

