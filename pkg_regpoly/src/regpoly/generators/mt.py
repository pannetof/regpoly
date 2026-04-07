"""
mt.py — Mersenne Twister (MT) generator.

The state consists of r words of w bits (total w*r bits), with an offset
parameter p so that the degree is k = w*r - p.  The recurrence uses a
twist matrix defined by coefficient a.

Parameters: w (word size), r (register count), m (middle distance),
p (offset), a (twist matrix coefficient).

Implements the Generateur abstract base class.
"""

import cython
import sys

from regpoly.bitvect import BitVect
from regpoly.generateur import Generateur


class MersenneTwister(Generateur):
    """
    Mersenne Twister generator.

    Attributes
    ----------
    k    : int      — degree = w * r - p  (inherited)
    w    : int      — word size in bits (≤ 32)
    r    : int      — number of w-bit words in the state register
    m    : int      — middle distance parameter
    p    : int      — offset parameter (0 ≤ p < w)
    a    : int      — twist matrix coefficient (w bits)
    i    : int      — current index into the circular buffer
    uu   : int      — upper mask: maskw << p
    ll   : int      — lower mask: (1 << p) - 1, or 0 if p == 0
    maskw: int      — w-bit mask: (1 << w) - 1
    L    : int      — output resolution in bits
    """

    def __init__(
        self,
        w: cython.int,
        r: cython.int,
        m: cython.int,
        p: cython.int,
        a: object,
        L: cython.int,
    ) -> None:
        k: cython.int = w * r - p
        super().__init__(k)
        self.w: cython.int = w
        self.r: cython.int = r
        self.m: cython.int = m
        self.p: cython.int = p
        self.a: object = a
        self.i: cython.int = 0
        self.L = L
        self.smax = sys.maxsize
        self.maskw: object = (1 << w) - 1
        self.uu: object = 0
        self.ll: object = 0
        # State storage uses w*r bits (not k) to match C's circular buffer.
        # The extra p bits at the low end are always zero.
        self._state_bits: cython.int = w * r

    # -- Generateur interface ---------------------------------------------

    @classmethod
    def name(cls) -> str:
        return "Mersenne Twister"

    def display(self) -> None:
        """DispMT: print w, r, m, p and the twist matrix coefficient a."""
        print(f" w= {self.w:3d}  r= {self.r:3d}  m= {self.m:3d}"
              f"  p= {self.p:3d} a= {self.a:10x}  ")

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """
        InitMT: copy init into state, set up masks, return first L output bits.

        The state is allocated with w*r bits (not k = w*r-p) to match C's
        circular buffer layout.  The init vector has k bits; the extra p
        low-order bits of the state remain zero.

        The initial output is a straight copy (no circular rotation),
        matching C's InitMT.
        """
        self.gen_state = BitVect(self._state_bits, 0)
        self.gen_state.copy_part_from(init_bv, self.k)

        self.i = 0
        self.maskw = (1 << self.w) - 1
        self.uu = self.maskw << self.p

        if self.p != 0:
            self.ll = (1 << self.p) - 1
        else:
            self.ll = 0

        # Initial output: straight copy of first L bits (no rotation)
        sb: int = self._state_bits
        L: int = self.L
        if L <= sb:
            result: object = (self.gen_state._val >> (sb - L)) & ((1 << L) - 1)
        else:
            result = self.gen_state._val << (L - sb)
        return BitVect(L, result)

    def __next__(self) -> BitVect:
        """
        MT: one recurrence step.

        Algorithm (matching the C code):
          Y = (V(i % r) & uu) | (V((i+1) % r) & ll)
          if Y is odd:
              V_i = V((i+m) % r) ^ (Y >> 1) ^ a
          else:
              V_i = V((i+m) % r) ^ (Y >> 1)
          SetV(i, V_i)
          output L bits starting from current position
          i = (i + 1) % r
        """
        w: cython.int = self.w
        r: cython.int = self.r
        m: cython.int = self.m
        ii: cython.int = self.i

        Y: object = (self._V(ii % r) & self.uu) | (self._V((ii + 1) % r) & self.ll)

        if Y & 1:
            V_i: object = self._V((ii + m) % r) ^ (Y >> 1) ^ self.a
        else:
            V_i = self._V((ii + m) % r) ^ (Y >> 1)

        self._SetV(ii, V_i)

        retour: BitVect = self._extract_output()

        self.i = (ii + 1) % r
        return retour

    def copy(self) -> "MersenneTwister":
        new: MersenneTwister = MersenneTwister(
            self.w, self.r, self.m, self.p, self.a, self.L
        )
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        new.i = self.i
        new.uu = self.uu
        new.ll = self.ll
        new.maskw = self.maskw
        new._state_bits = self._state_bits
        return new

    # transition_matrix is inherited from Generateur — no override needed
    # since initialize_state handles the _state_bits allocation internally.

    # -- Word access (matching C's V / SetV) ------------------------------

    def _V(self, idx: cython.int) -> object:
        """
        V(Gen, i): return the i-th w-bit word from the state.

        The C code stores words in reverse order: word i is at position
        (r - i - 1) in the register, packed into uint32_t words with WL=32.
        In Python, the BitVect stores bits MSB-first: bit 0 is the MSB.
        Word (r-i-1) starts at public bit (r-i-1)*w.
        """
        pos: cython.int = (self.r - idx - 1) * self.w
        shift: cython.int = self._state_bits - pos - self.w
        return (self.gen_state._val >> shift) & self.maskw

    def _SetV(self, idx: cython.int, val: object) -> None:
        """SetV(Gen, i, val): set the i-th w-bit word of the state."""
        pos: cython.int = (self.r - idx - 1) * self.w
        shift: cython.int = self._state_bits - pos - self.w
        mask: object = self.maskw << shift
        val_masked: object = (val & self.maskw) << shift
        self.gen_state._val = (self.gen_state._val & ~mask) | val_masked

    def _extract_output(self) -> BitVect:
        """
        Extract L output bits from the circular buffer.

        C's MT iterate outputs a circularly rotated view of the state,
        starting from word TI (the word just updated).  This corresponds
        to a right circular rotation of the packed state integer by
        (TI + 1) * w bits.

        Called from __next__ BEFORE self.i is incremented.
        """
        total: int = self._state_bits
        w: int = self.w
        L: int = self.L
        sv: object = self.gen_state._val
        rotation: int = ((self.i + 1) * w) % total
        if rotation != 0:
            sv = ((sv >> rotation) | (sv << (total - rotation))) & ((1 << total) - 1)
        if L <= total:
            result: object = (sv >> (total - L)) & ((1 << L) - 1)
        else:
            result = sv << (L - total)
        return BitVect(L, result)

    def _transition_state(self) -> object:
        """
        Return the k-bit state for transition matrix computation.

        Same circular rotation as _extract_output, but always returns k bits
        and is called AFTER self.i has been incremented (so uses self.i * w
        as the rotation amount).
        """
        total: int = self._state_bits
        w: int = self.w
        sv: object = self.gen_state._val
        rotation: int = (self.i * w) % total
        if rotation != 0:
            sv = ((sv >> rotation) | (sv << (total - rotation))) & ((1 << total) - 1)
        K: int = self.k
        return (sv >> (total - K)) & ((1 << K) - 1)

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateListFromFile(cls, filename: str, L: cython.int) -> list:
        """
        ReadDataMT: read MT parameters from filename.

        File format
        -----------
        Line 1 : type tag  ("MT")
        Line 2 : nbgen
        Lines 3 .. nbgen+2:
            w  r  m  p  hex_a

        Returns a list of MersenneTwister instances.
        """
        with open(filename) as f:
            f.readline()                        # skip type tag
            nbgen: cython.int = int(f.readline().split()[0])

            generators: list = []
            for _ in range(nbgen):
                tokens: list = f.readline().split()
                w: cython.int = int(tokens[0])
                r: cython.int = int(tokens[1])
                m: cython.int = int(tokens[2])
                p: cython.int = int(tokens[3])
                a_val: object = int(tokens[4], 16)

                if w > 32:
                    raise ValueError(f"w must be <= 32 bits, got {w}")

                gen_L: cython.int = L
                generators.append(cls(w, r, m, p, a_val, gen_L))

        return generators
