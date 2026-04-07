"""
test_bitvect.py — Complete test battery for BitVect, parametrised by size.

Run with default size (8):
    pytest test_bitvect.py -v

Run with multiple sizes:
    pytest test_bitvect.py -v --sizes 8,16,32,64,94
"""

import pytest
from bitvect import BitVect


# ------------------------------------------------------------------ #
#  Helpers                                                            #
# ------------------------------------------------------------------ #

def alt_val(n):
    """n-bit alternating pattern, MSB = 1: 10101..."""
    v = 0
    for i in range(0, n, 2):
        v |= 1 << (n - 1 - i)
    return v


def alt_inv_val(n):
    """n-bit alternating pattern, MSB = 0: 01010..."""
    return alt_val(n) ^ ((1 << n) - 1)


def top_k(n, k):
    """Value with top k bits set (0..k-1 in public indexing)."""
    k = max(0, min(k, n))
    if k == 0:
        return 0
    return ((1 << k) - 1) << (n - k)


def all_ones(n):
    return (1 << n) - 1


# ================================================================== #
#  Construction / factory methods                                      #
# ================================================================== #

class TestConstruction:
    def test_basic_init(self, size):
        v = alt_val(size)
        bv = BitVect(size, v)
        assert bv.n == size
        assert bv._val == v

    def test_init_masks_excess_bits(self, size):
        bv = BitVect(size, all_ones(size + 8))
        assert bv._val == all_ones(size)

    def test_init_default_val(self, size):
        bv = BitVect(size)
        assert bv._val == 0

    def test_zeros(self, size):
        bv = BitVect.zeros(size)
        assert bv.n == size and bv._val == 0

    def test_canonic_msb(self, size):
        bv = BitVect.canonic(size, 0)
        assert bv._val == 1 << (size - 1)

    def test_canonic_lsb(self, size):
        bv = BitVect.canonic(size, size - 1)
        assert bv._val == 1

    def test_canonic_mid(self, size):
        mid = size // 2
        bv = BitVect.canonic(size, mid)
        assert bv._val == 1 << (size - 1 - mid)

    def test_canonic_out_of_range(self):
        with pytest.raises(IndexError):
            BitVect.canonic(8, 8)
        with pytest.raises(IndexError):
            BitVect.canonic(8, -1)

    def test_init_zero_n_raises(self):
        with pytest.raises(ValueError):
            BitVect(0, 0)

    def test_all_ones(self, size):
        bv = BitVect.all_ones(size)
        assert bv._val == all_ones(size)

    def test_mask(self, size):
        k = size // 2
        bv = BitVect.mask(size, k)
        assert bv._val == top_k(size, k)
        assert BitVect.mask(size, 0)._val == 0
        assert BitVect.mask(size, size)._val == all_ones(size)
        assert BitVect.mask(size, size + 4)._val == all_ones(size)

    def test_invmask(self, size):
        k = size // 2
        bv = BitVect.invmask(size, k)
        assert bv._val == all_ones(size) - top_k(size, k)
        assert BitVect.invmask(size, 0)._val == all_ones(size)
        assert BitVect.invmask(size, size)._val == 0
        assert BitVect.invmask(size, size + 4)._val == 0

    def test_mask_invmask_complement(self, size):
        """mask(n, l) XOR invmask(n, l) should be all_ones(n)."""
        for l in [0, size // 4, size // 2, 3 * size // 4, size]:
            m = BitVect.mask(size, l)
            im = BitVect.invmask(size, l)
            assert (m ^ im)._val == all_ones(size)

    def test_random_width(self, size):
        bv = BitVect.random(size)
        assert bv.n == size
        assert 0 <= bv._val <= all_ones(size)


# ================================================================== #
#  In-place reset / canonic / mask helpers                            #
# ================================================================== #

class TestInPlaceReset:
    def test_put_to_zero(self, size):
        bv = BitVect(size, all_ones(size))
        bv.put_to_zero()
        assert bv._val == 0

    def test_put_to_all_ones(self, size):
        bv = BitVect(size, 0)
        bv.put_to_all_ones()
        assert bv._val == all_ones(size)

    def test_set_canonic(self, size):
        bv = BitVect(size, all_ones(size))
        bv.set_canonic(size // 2)
        assert bv._val == 1 << (size - 1 - size // 2)

    def test_set_canonic_out_of_range(self):
        with pytest.raises(IndexError):
            BitVect(8, 0).set_canonic(8)

    def test_set_mask(self, size):
        k = size // 3
        bv = BitVect(size, all_ones(size))
        bv.set_mask(k)
        assert bv._val == top_k(size, k)
        bv.set_mask(0)
        assert bv._val == 0
        bv.set_mask(size + 4)
        assert bv._val == all_ones(size)

    def test_set_invmask(self, size):
        k = size // 3
        bv = BitVect(size, all_ones(size))
        bv.set_invmask(k)
        assert bv._val == all_ones(size) - top_k(size, k)
        bv.set_invmask(0)
        assert bv._val == all_ones(size)
        bv.set_invmask(size + 4)
        assert bv._val == 0


# ================================================================== #
#  Bit access                                                          #
# ================================================================== #

class TestBitAccess:
    def test_get_bit_msb(self, size):
        bv = BitVect(size, 1 << (size - 1))  # only MSB set
        assert bv.get_bit(0) == 1
        assert bv.get_bit(size - 1) == 0

    def test_get_bit_lsb(self, size):
        bv = BitVect(size, 1)  # only LSB set
        assert bv.get_bit(0) == 0
        assert bv.get_bit(size - 1) == 1

    def test_get_bit_alternating(self, size):
        bv = BitVect(size, alt_val(size))
        for i in range(size):
            assert bv.get_bit(i) == (1 if i % 2 == 0 else 0)

    def test_put_bit(self, size):
        bv = BitVect(size, 0)
        bv.put_bit(0, 1)
        assert bv._val == 1 << (size - 1)
        bv.put_bit(size - 1, 1)
        assert bv._val == (1 << (size - 1)) | 1
        bv.put_bit(0, 0)
        assert bv._val == 1

    def test_getitem_setitem(self, size):
        bv = BitVect(size, 0)
        bv[size // 2] = 1
        assert bv[size // 2] == 1
        assert bv[0] == 0

    def test_bit_out_of_range(self):
        bv = BitVect(8, 0)
        with pytest.raises(IndexError):
            bv.get_bit(8)
        with pytest.raises(IndexError):
            bv.get_bit(-1)
        with pytest.raises(IndexError):
            bv.put_bit(8, 1)

    def test_len(self, size):
        assert len(BitVect(size, 0)) == size


# ================================================================== #
#  Predicates                                                          #
# ================================================================== #

class TestPredicates:
    def test_have_common_bits_true(self, size):
        a = BitVect(size, 1 << (size - 1))   # MSB
        b = BitVect(size, all_ones(size))      # all set
        assert a.have_common_bits(b)

    def test_have_common_bits_false(self, size):
        a = BitVect(size, top_k(size, size // 2))      # top half
        b = BitVect(size, all_ones(size) - top_k(size, size // 2))  # bottom half
        assert not a.have_common_bits(b)

    def test_have_common_bits_size_mismatch(self):
        with pytest.raises(ValueError):
            BitVect(8, 0).have_common_bits(BitVect(4, 0))


# ================================================================== #
#  Copy                                                                #
# ================================================================== #

class TestCopy:
    def test_copy_independence(self, size):
        v = alt_val(size)
        a = BitVect(size, v)
        b = a.copy()
        assert a._val == b._val
        b._val = 0
        assert a._val == v

    def test_copy_from(self, size):
        a = BitVect(size, alt_val(size))
        b = BitVect(size, 0)
        b.copy_from(a)
        assert b._val == a._val
        a._val = 0
        assert b._val == alt_val(size)

    def test_copy_from_size_mismatch(self):
        with pytest.raises(ValueError):
            BitVect(8, 0).copy_from(BitVect(4, 0))

    def test_copy_part_from_half(self, size):
        k = size // 2
        src = BitVect(size, all_ones(size))
        dst = BitVect(size, 0)
        dst.copy_part_from(src, k)
        # Top k bits from src (all 1s) placed in top k of dst; rest zeroed.
        assert dst._val == top_k(size, k)

    def test_copy_part_from_full(self, size):
        src = BitVect(size, alt_val(size))
        dst = BitVect(size, 0)
        dst.copy_part_from(src, size)
        assert dst._val == src._val

    def test_copy_part_from_zero_l(self, size):
        dst = BitVect(size, all_ones(size))
        before = dst._val
        dst.copy_part_from(BitVect(size, 0), 0)
        assert dst._val == before  # unchanged

    def test_copy_part_from_dst_too_small(self):
        with pytest.raises(ValueError):
            BitVect(4, 0).copy_part_from(BitVect(8, 0xFF), 8)


# ================================================================== #
#  Boolean operations — new BitVect                                   #
# ================================================================== #

class TestBooleanNew:
    def test_xor(self, size):
        a = BitVect(size, alt_val(size))
        b = BitVect(size, alt_inv_val(size))
        assert (a ^ b)._val == all_ones(size)

    def test_xor_self(self, size):
        a = BitVect(size, alt_val(size))
        assert (a ^ a)._val == 0

    def test_xor_size_mismatch(self):
        with pytest.raises(ValueError):
            BitVect(8, 0) ^ BitVect(4, 0)

    def test_and(self, size):
        a = BitVect(size, alt_val(size))
        b = BitVect(size, all_ones(size))
        assert (a & b)._val == a._val

    def test_and_complement(self, size):
        a = BitVect(size, alt_val(size))
        b = BitVect(size, alt_inv_val(size))
        assert (a & b)._val == 0

    def test_and_size_mismatch(self):
        with pytest.raises(ValueError):
            BitVect(8, 0) & BitVect(4, 0)

    def test_invert(self, size):
        bv = BitVect(size, alt_val(size))
        assert (~bv)._val == alt_inv_val(size)
        assert (~BitVect(size, 0))._val == all_ones(size)
        assert (~BitVect.all_ones(size))._val == 0

    def test_double_invert(self, size):
        bv = BitVect(size, alt_val(size))
        assert (~~bv)._val == bv._val

    def test_xor3(self, size):
        a = BitVect(size, alt_val(size))
        b = BitVect(size, alt_inv_val(size))
        c = BitVect(size, top_k(size, size // 2))
        expected = alt_val(size) ^ alt_inv_val(size) ^ top_k(size, size // 2)
        assert a.xor3(b, c)._val == expected & all_ones(size)

    def test_xor3_size_mismatch(self):
        with pytest.raises(ValueError):
            BitVect(8, 0).xor3(BitVect(4, 0), BitVect(8, 0))

    def test_and_mask(self, size):
        bv = BitVect.all_ones(size)
        k = size // 2
        result = bv.and_mask(k)
        assert result._val == top_k(size, k)
        assert bv.and_mask(0)._val == 0
        assert bv.and_mask(size)._val == bv._val
        assert bv.and_mask(size + 4)._val == bv._val

    def test_and_invmask(self, size):
        bv = BitVect.all_ones(size)
        k = size // 2
        result = bv.and_invmask(k)
        assert result._val == all_ones(size) - top_k(size, k)
        assert bv.and_invmask(0)._val == bv._val
        assert bv.and_invmask(size)._val == 0
        assert bv.and_invmask(size + 4)._val == 0


# ================================================================== #
#  Boolean in-place operations                                        #
# ================================================================== #

class TestBooleanInPlace:
    def test_ixor(self, size):
        a = BitVect(size, alt_val(size))
        a ^= BitVect(size, alt_inv_val(size))
        assert a._val == all_ones(size)

    def test_ixor_size_mismatch(self):
        a = BitVect(8, 0)
        with pytest.raises(ValueError):
            a ^= BitVect(4, 0)

    def test_iand(self, size):
        a = BitVect(size, alt_val(size))
        a &= BitVect.all_ones(size)
        assert a._val == alt_val(size)

    def test_iand_zeroes_out(self, size):
        a = BitVect(size, alt_val(size))
        a &= BitVect(size, alt_inv_val(size))
        assert a._val == 0

    def test_iand_size_mismatch(self):
        a = BitVect(8, 0)
        with pytest.raises(ValueError):
            a &= BitVect(4, 0)

    def test_iand_mask(self, size):
        k = size // 2
        bv = BitVect.all_ones(size)
        bv.iand_mask(k)
        assert bv._val == top_k(size, k)
        bv2 = BitVect.all_ones(size)
        bv2.iand_mask(0)
        assert bv2._val == 0
        bv3 = BitVect.all_ones(size)
        bv3.iand_mask(size + 4)
        assert bv3._val == all_ones(size)

    def test_iand_invmask(self, size):
        k = size // 2
        bv = BitVect.all_ones(size)
        bv.iand_invmask(k)
        assert bv._val == all_ones(size) - top_k(size, k)
        bv2 = BitVect.all_ones(size)
        bv2.iand_invmask(0)
        assert bv2._val == all_ones(size)
        bv3 = BitVect.all_ones(size)
        bv3.iand_invmask(size + 4)
        assert bv3._val == 0

    def test_inverse_self(self, size):
        bv = BitVect(size, alt_val(size))
        bv.inverse_self()
        assert bv._val == alt_inv_val(size)
        bv.inverse_self()
        assert bv._val == alt_val(size)


# ================================================================== #
#  Shift operations — new BitVect                                     #
# ================================================================== #

class TestShiftNew:
    def test_lshift_moves_toward_msb(self, size):
        k = size // 4
        # LSB region shifted left should now be in MSB region
        bv = BitVect(size, top_k(size, size // 2))  # top half set
        result = bv << k
        assert result._val == top_k(size, size // 2 - k)

    def test_lshift_zero(self, size):
        bv = BitVect(size, alt_val(size))
        assert (bv << 0)._val == bv._val

    def test_lshift_by_n(self, size):
        bv = BitVect.all_ones(size)
        assert (bv << size)._val == 0

    def test_lshift_beyond_n(self, size):
        bv = BitVect.all_ones(size)
        assert (bv << (size + 5))._val == 0

    def test_lshift_msb_lost(self, size):
        bv = BitVect(size, 1 << (size - 1))  # only MSB
        assert (bv << 1)._val == 0

    def test_rshift_moves_toward_lsb(self, size):
        k = size // 4
        bv = BitVect(size, top_k(size, size // 2))
        result = bv >> k
        assert result._val == top_k(size, size // 2 + k) - top_k(size, k)

    def test_rshift_zero(self, size):
        bv = BitVect(size, alt_val(size))
        assert (bv >> 0)._val == bv._val

    def test_rshift_by_n(self, size):
        bv = BitVect.all_ones(size)
        assert (bv >> size)._val == 0

    def test_rshift_beyond_n(self, size):
        bv = BitVect.all_ones(size)
        assert (bv >> (size + 5))._val == 0

    def test_rshift_lsb_lost(self, size):
        bv = BitVect(size, 1)  # only LSB
        assert (bv >> 1)._val == 0

    def test_lshift_rshift_inner_bits_preserved(self, size):
        """Left then right by k preserves bits not touching edges."""
        k = size // 4
        # Start with bits only in the middle quartile
        inner = top_k(size, 3 * size // 4) - top_k(size, size // 4)
        bv = BitVect(size, inner)
        assert ((bv << k) >> k)._val == bv._val

    def test_width_never_exceeded(self, size):
        """No shift result should produce _val >= 2^size."""
        bv = BitVect.all_ones(size)
        for k in [0, 1, size // 2, size - 1, size, size + 5]:
            assert (bv << k)._val <= all_ones(size)
            assert (bv >> k)._val <= all_ones(size)


# ================================================================== #
#  Shift in-place operations                                          #
# ================================================================== #

class TestShiftInPlace:
    def test_ilshift(self, size):
        k = size // 4
        bv = BitVect(size, top_k(size, size // 2))
        bv <<= k
        assert bv._val == top_k(size, size // 2 - k)

    def test_ilshift_zero(self, size):
        v = alt_val(size)
        bv = BitVect(size, v)
        bv <<= 0
        assert bv._val == v

    def test_ilshift_by_n(self, size):
        bv = BitVect.all_ones(size)
        bv <<= size
        assert bv._val == 0

    def test_irshift(self, size):
        k = size // 4
        bv = BitVect(size, top_k(size, size // 2))
        bv >>= k
        assert bv._val == top_k(size, size // 2 + k) - top_k(size, k)

    def test_irshift_zero(self, size):
        v = alt_val(size)
        bv = BitVect(size, v)
        bv >>= 0
        assert bv._val == v

    def test_irshift_by_n(self, size):
        bv = BitVect.all_ones(size)
        bv >>= size
        assert bv._val == 0

    def test_lshift1_self(self, size):
        """Single-step left shift in place."""
        bv = BitVect(size, 1 << (size - 2))  # second-to-last bit (public index size-2)
        bv.lshift1_self()
        assert bv._val == 1 << (size - 1)  # now MSB

    def test_lshift1_self_msb_lost(self, size):
        bv = BitVect(size, 1 << (size - 1))  # MSB only
        bv.lshift1_self()
        assert bv._val == 0


# ================================================================== #
#  Rotative shifts                                                     #
# ================================================================== #

class TestRotativeShifts:
    def test_lrot_full_width_msb_wraps(self, size):
        """Rotating all bits left by 1: MSB wraps to LSB."""
        bv = BitVect(size, 1 << (size - 1))  # MSB set
        result = bv.lrot_shift(1, size)
        assert result._val == 1  # LSB

    def test_rot_zero_shift(self, size):
        bv = BitVect(size, alt_val(size))
        w = size // 2
        assert bv.lrot_shift(0, w)._val == bv._val

    def test_rot_full_rotation_is_identity(self, size):
        """Rotating by w positions (full window) returns original."""
        w = size // 2
        bv = BitVect(size, alt_val(size))
        assert bv.lrot_shift(w, w)._val == bv._val

    def test_rot_bottom_bits_unchanged(self, size):
        """Bits outside the rotation window are not modified."""
        w = size // 2
        bv = BitVect(size, alt_val(size))
        bottom_mask = all_ones(size) - top_k(size, w)
        bottom = bv._val & bottom_mask
        result = bv.lrot_shift(size // 4, w)
        assert (result._val & bottom_mask) == bottom

    def test_rot_partial_window(self, size):
        """Rotate only the top quarter of bits left, leave rest alone."""
        w = max(2, size // 4)
        k = 1
        bv = BitVect(size, all_ones(size))  # all ones — rotation of all-ones is still all-ones
        result = bv.lrot_shift(k, w)
        assert result._val == bv._val  # all-ones is invariant under any rotation

    def test_rot_invalid_w(self):
        bv = BitVect(8, 0)
        with pytest.raises(ValueError):
            bv.lrot_shift(1, 0)
        with pytest.raises(ValueError):
            bv.lrot_shift(1, 9)


# ================================================================== #
#  Display / repr / str                                               #
# ================================================================== #

class TestDisplay:
    def test_str_length(self, size):
        bv = BitVect(size, alt_val(size))
        assert len(str(bv)) == size

    def test_str_only_binary_chars(self, size):
        bv = BitVect(size, alt_val(size))
        assert all(c in "01" for c in str(bv))

    def test_str_all_zeros(self, size):
        assert str(BitVect(size, 0)) == "0" * size

    def test_str_all_ones(self, size):
        assert str(BitVect.all_ones(size)) == "1" * size

    def test_str_alternating(self, size):
        bv = BitVect(size, alt_val(size))
        for i, ch in enumerate(str(bv)):
            assert ch == ("1" if i % 2 == 0 else "0")

    def test_str_leading_zeros_lsb_only(self, size):
        bv = BitVect(size, 1)  # only LSB
        assert str(bv) == "0" * (size - 1) + "1"

    def test_repr_contains_n(self, size):
        bv = BitVect(size, alt_val(size))
        assert f"n={size}" in repr(bv)

    def test_display_stdout(self, size, capsys):
        bv = BitVect(size, alt_val(size))
        bv.display()
        out = capsys.readouterr().out
        assert len(out) == size
        assert all(c in "01" for c in out)

    def test_display_partial(self, size, capsys):
        k = size // 2
        bv = BitVect(size, all_ones(size))
        bv.display(l=k)
        out = capsys.readouterr().out
        assert len(out) == k
        assert out == "1" * k



# ================================================================== #
#  Edge cases (fixed, no size parametrization)                       #
# ================================================================== #

class TestEdgeCases:
    def test_width_1(self):
        bv0 = BitVect(1, 0)
        bv1 = BitVect(1, 1)
        assert str(bv0) == "0"
        assert str(bv1) == "1"
        assert (bv0 ^ bv1)._val == 1
        assert (~bv0)._val == 1
        assert (~bv1)._val == 0

    def test_all_operations_preserve_width(self, size):
        """No operation should ever produce _val >= 2^n."""
        a = BitVect.all_ones(size)
        b = BitVect(size, alt_val(size))
        k = size // 4
        w = max(2, size // 2)
        for result in [
            a ^ b, a & b, ~a,
            a << k, a >> k,
            a.and_mask(k), a.and_invmask(k),
            a.xor3(b, BitVect(size, alt_inv_val(size))),
            a.copy(),
            a.lrot_shift(k, w),
        ]:
            assert result._val <= all_ones(size), f"Width overflow for size={size}"
