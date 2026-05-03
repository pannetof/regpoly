"""
bitvect.py — BitVect class implementing all bit-vector operations from vectorsF2.c.

Bit-indexing convention (matches the C source):
    bit 0   = most significant bit (MSB, leftmost when printed)
    bit n-1 = least significant bit (LSB, rightmost when printed)

Internal storage:
    A single Python int `_val` where Python bit position (n-1-i) holds public bit i.
    This mapping makes Python's << correspond to C's BVLShift (shift toward bit 0)
    and Python's >> correspond to C's BVRShift (shift toward bit n-1).

    Example for n=8, value "10110010":
        public bits:  0 1 2 3 4 5 6 7
        bit values:   1 0 1 1 0 0 1 0
        _val = 0b10110010 = 178  (bit 7 of _val = public bit 0 = 1, etc.)

"""

import random as _random


class BitVect:
    """
    Fixed-width bit vector over GF(2) with big-endian bit indexing (bit 0 = MSB).

    All operations that modify the vector can be done either in place (methods
    ending in _self or using augmented assignment operators) or as pure functions
    that return a new BitVect.
    """

    # ------------------------------------------------------------------ #
    #  Construction / factory methods                                      #
    # ------------------------------------------------------------------ #

    def __init__(self, n: int, val: int = 0):
        """
        Create a bit vector of n bits initialised to val.

        val is interpreted as a Python int where Python bit 0 is public bit n-1
        (i.e., val=1 sets only the LSB, public bit n-1).
        Bits above position n-1 in val are silently masked out.
        """
        if n <= 0:
            raise ValueError(f"n must be positive, got {n}")
        self.n = n
        self._val = val & ((1 << n) - 1)

    @staticmethod
    def zeros(n: int) -> "BitVect":
        """AllocBV / PutBVToZero: create a zero-filled bit vector of n bits."""
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = n
        bv._val = 0
        return bv

    @staticmethod
    def canonic(n: int, l: int) -> "BitVect":
        """
        BVCanonic: create an n-bit vector with only bit l set to 1.
        Bit 0 is the MSB; bit l corresponds to Python bit position n-1-l.
        """
        if l < 0 or l >= n:
            raise IndexError(f"Bit index {l} out of range for n={n}")
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = n
        bv._val = 1 << (n - 1 - l)
        return bv

    @staticmethod
    def all_ones(n: int) -> "BitVect":
        """AllOnes: create an n-bit vector with every bit set to 1."""
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = n
        bv._val = (1 << n) - 1
        return bv

    @staticmethod
    def mask(n: int, l: int) -> "BitVect":
        """
        mask(): create an n-bit vector with the first l bits (bits 0..l-1) set to 1
        and the remaining bits set to 0.
        """
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = n
        if l <= 0:
            bv._val = 0
        elif l >= n:
            bv._val = (1 << n) - 1
        else:
            bv._val = ((1 << l) - 1) << (n - l)
        return bv

    @staticmethod
    def invmask(n: int, l: int) -> "BitVect":
        """
        invmask(): create an n-bit vector with bits l..n-1 set to 1 and
        bits 0..l-1 set to 0.  Complement of mask().
        """
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = n
        if l <= 0:
            bv._val = (1 << n) - 1
        elif l >= n:
            bv._val = 0
        else:
            bv._val = (1 << (n - l)) - 1
        return bv

    @staticmethod
    def random(n: int) -> "BitVect":
        """RandVect: create a random n-bit vector."""
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = n
        bv._val = _random.getrandbits(n)
        return bv

    # ------------------------------------------------------------------ #
    #  In-place reset / canonic / mask helpers                            #
    # ------------------------------------------------------------------ #

    def put_to_zero(self) -> None:
        """PutBVToZero: set all bits to 0 in place."""
        self._val = 0

    def put_to_all_ones(self) -> None:
        """AllOnes (in place): set all bits to 1."""
        self._val = (1 << self.n) - 1

    def set_canonic(self, l: int) -> None:
        """BVCanonic (in place): zero the vector then set only bit l."""
        if l < 0 or l >= self.n:
            raise IndexError(f"Bit index {l} out of range for n={self.n}")
        self._val = 1 << (self.n - 1 - l)

    def set_mask(self, l: int) -> None:
        """mask() (in place): set first l bits to 1, rest to 0."""
        if l <= 0:
            self._val = 0
        elif l >= self.n:
            self._val = (1 << self.n) - 1
        else:
            self._val = ((1 << l) - 1) << (self.n - l)

    def set_invmask(self, l: int) -> None:
        """invmask() (in place): set bits l..n-1 to 1, bits 0..l-1 to 0."""
        if l <= 0:
            self._val = (1 << self.n) - 1
        elif l >= self.n:
            self._val = 0
        else:
            self._val = (1 << (self.n - l)) - 1

    # ------------------------------------------------------------------ #
    #  Bit access                                                          #
    # ------------------------------------------------------------------ #

    def get_bit(self, i: int) -> int:
        """ValBitBV: return the value (0 or 1) of bit i (bit 0 = MSB)."""
        if i < 0 or i >= self.n:
            raise IndexError(f"Bit index {i} out of range for n={self.n}")
        return (self._val >> (self.n - 1 - i)) & 1

    def put_bit(self, i: int, v: int) -> None:
        """PutBitBV: set bit i to v (0 or 1). Bit 0 = MSB."""
        if i < 0 or i >= self.n:
            raise IndexError(f"Bit index {i} out of range for n={self.n}")
        pos = self.n - 1 - i
        if v:
            self._val |= 1 << pos
        else:
            self._val &= ~(1 << pos)

    def __getitem__(self, i: int) -> int:
        """Shorthand for get_bit(i)."""
        return self.get_bit(i)

    def __setitem__(self, i: int, v: int) -> None:
        """Shorthand for put_bit(i, v)."""
        self.put_bit(i, v)

    def __len__(self) -> int:
        """Return the number of bits n."""
        return self.n

    # ------------------------------------------------------------------ #
    #  Predicates                                                          #
    # ------------------------------------------------------------------ #

    def have_common_bits(self, other: "BitVect") -> bool:
        """VerifBitsCommuns: return True if self and other share at least one set bit."""
        if self.n != other.n:
            raise ValueError(
                f"VerifBitsCommuns: size mismatch ({self.n} vs {other.n})"
            )
        return bool(self._val & other._val)

    # ------------------------------------------------------------------ #
    #  Copy                                                                #
    # ------------------------------------------------------------------ #

    def copy(self) -> "BitVect":
        """CopyBV: return a new BitVect with the same width and value."""
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        bv._val = self._val
        return bv

    def copy_from(self, other: "BitVect") -> None:
        """CopyBV (in place): copy other into self. Widths must match."""
        if self.n != other.n:
            raise ValueError(
                f"CopyBV: size mismatch ({self.n} vs {other.n})"
            )
        self._val = other._val

    def copy_part_from(self, other: "BitVect", l: int) -> None:
        """
        CopyBVPart: copy the first l bits of other into self (bits 0..l-1).
        self.n must be >= l.  Bits l..self.n-1 are preserved (not zeroed),
        matching the C behaviour when l is word-aligned.
        When l is not a multiple of 32 the C source masks the last word, which
        here corresponds to zeroing bits l..self.n-1.  We always apply the
        masking for simplicity (safe and consistent with the spirit of the C code).
        """
        if self.n < l:
            raise ValueError(
                f"CopyBVPart: destination (n={self.n}) too small for l={l} bits"
            )
        if l <= 0:
            return
        # Extract the top l bits of other.
        src_top: int
        if l >= other.n:
            src_top = other._val << (l - other.n)
        else:
            src_top = other._val >> (other.n - l)
        src_top &= (1 << l) - 1
        # Place those l bits at the top of self; zero the remaining bits.
        self._val = src_top << (self.n - l)

    # ------------------------------------------------------------------ #
    #  Boolean operations — return a new BitVect                          #
    # ------------------------------------------------------------------ #

    def __xor__(self, other: "BitVect") -> "BitVect":
        """XORBV: return self XOR other (new BitVect)."""
        if self.n != other.n:
            raise ValueError(f"XORBV: size mismatch ({self.n} vs {other.n})")
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        bv._val = self._val ^ other._val
        return bv

    def __and__(self, other: "BitVect") -> "BitVect":
        """ANDBV: return self AND other (new BitVect)."""
        if self.n != other.n:
            raise ValueError(f"ANDBV: size mismatch ({self.n} vs {other.n})")
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        bv._val = self._val & other._val
        return bv

    def __invert__(self) -> "BitVect":
        """InverseBV: return bitwise NOT of self (new BitVect)."""
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        bv._val = self._val ^ ((1 << self.n) - 1)
        return bv

    def xor3(self, b: "BitVect", c: "BitVect") -> "BitVect":
        """
        XOR2BV: return self XOR b XOR c (new BitVect).
        Equivalent to the three-operand XOR A = B ^ C ^ D in C.
        """
        if self.n != b.n or b.n != c.n:
            raise ValueError("XOR2BV: size mismatch")
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        bv._val = self._val ^ b._val ^ c._val
        return bv

    def and_mask(self, t: int) -> "BitVect":
        """
        ANDBVMask: return a copy of self with only the first t bits kept.
        Bits t..n-1 are zeroed in the result.
        """
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        if t <= 0:
            bv._val = 0
        elif t >= self.n:
            bv._val = self._val
        else:
            bv._val = self._val & (((1 << t) - 1) << (self.n - t))
        return bv

    def and_invmask(self, t: int) -> "BitVect":
        """
        ANDBVInvMask: return a copy of self with bits 0..t-1 zeroed
        (only bits t..n-1 kept).
        """
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        if t <= 0:
            bv._val = self._val
        elif t >= self.n:
            bv._val = 0
        else:
            bv._val = self._val & ((1 << (self.n - t)) - 1)
        return bv

    # ------------------------------------------------------------------ #
    #  Boolean in-place operations                                        #
    # ------------------------------------------------------------------ #

    def __ixor__(self, other: "BitVect") -> "BitVect":
        """XORBVSelf: self ^= other (in place)."""
        if self.n != other.n:
            raise ValueError(f"XORBVSelf: size mismatch ({self.n} vs {other.n})")
        self._val ^= other._val
        return self

    def __iand__(self, other: "BitVect") -> "BitVect":
        """ANDBVSelf: self &= other (in place)."""
        if self.n != other.n:
            raise ValueError(f"ANDBVSelf: size mismatch ({self.n} vs {other.n})")
        self._val &= other._val
        return self

    def iand_mask(self, t: int) -> None:
        """ANDBVMask (in place): keep only the first t bits of self."""
        if t <= 0:
            self._val = 0
        elif t < self.n:
            self._val &= ((1 << t) - 1) << (self.n - t)

    def iand_invmask(self, t: int) -> None:
        """ANDBVInvMask (in place): zero bits 0..t-1, keep bits t..n-1."""
        if t <= 0:
            pass
        elif t >= self.n:
            self._val = 0
        else:
            self._val &= (1 << (self.n - t)) - 1

    def inverse_self(self) -> None:
        """InverseBV (in place): self = ~self (bitwise NOT within n bits)."""
        self._val ^= (1 << self.n) - 1

    # ------------------------------------------------------------------ #
    #  Shift operations — return a new BitVect                            #
    # ------------------------------------------------------------------ #

    def __lshift__(self, k: int) -> "BitVect":
        """
        BVLShift: logical left shift by k positions (toward bit 0 / MSB).
        Bits shifted beyond bit 0 are lost; zeros enter from bit n-1.
        """
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        if k <= 0:
            bv._val = self._val
        elif k >= self.n:
            bv._val = 0
        else:
            bv._val = (self._val << k) & ((1 << self.n) - 1)
        return bv

    def __rshift__(self, k: int) -> "BitVect":
        """
        BVRShift: logical right shift by k positions (toward bit n-1 / LSB).
        Bits shifted beyond bit n-1 are lost; zeros enter from bit 0.
        """
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        if k <= 0:
            bv._val = self._val
        elif k >= self.n:
            bv._val = 0
        else:
            bv._val = self._val >> k
        return bv

    # ------------------------------------------------------------------ #
    #  Shift in-place operations                                          #
    # ------------------------------------------------------------------ #

    def __ilshift__(self, k: int) -> "BitVect":
        """BVLShiftSelf: self <<= k (in place)."""
        if k <= 0:
            return self
        if k >= self.n:
            self._val = 0
        else:
            self._val = (self._val << k) & ((1 << self.n) - 1)
        return self

    def __irshift__(self, k: int) -> "BitVect":
        """BVRShiftSelf: self >>= k (in place)."""
        if k <= 0:
            return self
        if k >= self.n:
            self._val = 0
        else:
            self._val >>= k
        return self

    def lshift1_self(self) -> None:
        """
        BVLS1Self: self <<= 1 (in place).
        Specialised single-bit left shift; same result as <<= 1 but named
        explicitly to match the C function used in hot inner loops.
        """
        self._val = (self._val << 1) & ((1 << self.n) - 1)

    # ------------------------------------------------------------------ #
    #  Rotative shifts — return a new BitVect                             #
    # ------------------------------------------------------------------ #

    def lrot_shift(self, k: int, w: int) -> "BitVect":
        """
        BVLRotativeShift: rotate the first w bits left by k positions,
        leaving bits w..n-1 unchanged.

        'Left' means toward bit 0 (MSB): bit i moves to bit (i-k) % w for i < w.
        """
        if w <= 0 or w > self.n:
            raise ValueError(f"lrot_shift: w={w} out of range for n={self.n}")
        k = k % w
        if k == 0:
            return self.copy()
        # Split: top w public bits (Python bits n-1..n-w) and bottom n-w bits.
        top: int = self._val >> (self.n - w)
        bottom: int = self._val & ((1 << (self.n - w)) - 1)
        rotated: int = ((top << k) | (top >> (w - k))) & ((1 << w) - 1)
        bv: BitVect = BitVect.__new__(BitVect)
        bv.n = self.n
        bv._val = (rotated << (self.n - w)) | bottom
        return bv

    # ------------------------------------------------------------------ #
    #  Display                                                             #
    # ------------------------------------------------------------------ #

    def __repr__(self) -> str:
        hex_digits: int = (self.n + 3) // 4
        return f"BitVect(n={self.n}, val=0x{self._val:0{hex_digits}x})"

    def __str__(self) -> str:
        """Return the bit string from bit 0 (MSB) to bit n-1 (LSB)."""
        return format(self._val, f"0{self.n}b")

    def display(self, l: int = -1) -> None:
        """DispBitVect: print the first l bits (default: all n bits)."""
        if l < 0:
            l = self.n
        print("".join(str(self.get_bit(i)) for i in range(l)), end="")
