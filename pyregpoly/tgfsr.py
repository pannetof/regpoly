"""
tgfsr.py — Twisted GFSR (TGFSR) generator.

The state consists of r words of w bits (k = w*r total).
Each step applies a twist matrix A to two state words at prescribed
distances, then shifts and XORs the result into the state register.

Characteristic polynomial: computed via Berlekamp-Massey (inherited).

Implements the Generateur abstract base class.
"""

from __future__ import annotations

import sys

from bitvect import BitVect
from generateur import Generateur


class TGFSR(Generateur):
    """
    Twisted Generalized Feedback Shift Register.

    Attributes
    ----------
    k    : int     — degree = w * r  (inherited)
    w    : int     — word size in bits
    r    : int     — number of w-bit words in the state register
    m    : int     — middle distance parameter
    a    : BitVect — twist matrix A (k bits, only first w bits used)
    L    : int     — output resolution in bits
    """

    def __init__(self, w: int, r: int, m: int, a: BitVect, L: int) -> None:
        k = w * r
        super().__init__(k)
        self.w: int = w
        self.r: int = r
        self.m: int = m
        self.a: BitVect = a.copy()
        self.L = L
        self.step = 1
        self.smax = sys.maxsize

    # -- Generateur interface ---------------------------------------------

    @classmethod
    def name(cls) -> str:
        return "TGFSR"

    def display(self) -> None:
        """DispTGFSR: print w, r, m and the twist matrix a."""
        a_word = self.a._val >> (self.a.n - 32) if self.a.n >= 32 else self.a._val
        print(f" w= {self.w:3d}  r= {self.r:3d}  m= {self.m:3d}  a= {a_word:08x}")

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """InitTGFSR: copy init into state, return first L bits."""
        self.gen_state = BitVect(self.k, 0)
        self.gen_state.copy_part_from(init_bv, self.k)
        retour = BitVect(self.L, 0)
        retour.copy_part_from(init_bv, self.L)
        return retour

    def __next__(self) -> BitVect:
        """
        TGFSR: one recurrence step.

        Algorithm (matching the C code):
          Temp2 = (state << w*(r-m-1)) masked to w bits   [word at distance m+1]
          Temp  = (state << w*(r-1))   masked to w bits   [last word]
          if bit w-1 of Temp == 0:
              Temp = Temp >> 1  XOR  Temp2
          else:
              Temp = Temp >> 1  XOR  Temp2  XOR  A
          Temp masked to w bits
          state = (state >> w)  XOR  Temp
          return first L bits of state
        """
        w = self.w
        r = self.r
        m = self.m
        k = self.k
        state = self.gen_state

        # Extract word at distance (r-m-1) and last word (r-1)
        temp2 = (state << (w * (r - m - 1))).and_mask(w)
        temp  = (state << (w * (r - 1))).and_mask(w)

        # Check LSB of the w-bit word (bit w-1 in public indexing)
        if temp.get_bit(w - 1) == 0:
            temp >>= 1
            temp = temp2 ^ temp
        else:
            temp >>= 1
            temp = temp2.xor3(temp, self.a)

        # Shift state right by w, XOR with masked temp
        state >>= w
        temp = temp.and_mask(w)
        self.gen_state = state ^ temp

        # Return first L bits
        retour = BitVect(self.L, 0)
        retour.copy_part_from(self.gen_state, self.L)
        return retour

    def copy(self) -> "TGFSR":
        new = TGFSR(self.w, self.r, self.m, self.a, self.L)
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        new.step = self.step
        return new

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateListFromFile(cls, filename: str, L: int) -> list["TGFSR"]:
        """
        ReadDataTGFSR: read TGFSR parameters from filename.

        File format
        -----------
        Line 1 : type tag  ("tgfsr")
        Line 2 : w  r
        Line 3 : nbgen
        Lines 4 .. nbgen+3:
            hex_a  m
            hex_a is the twist matrix value (stored in vect[0])
            m is the middle distance parameter

        Returns a list of TGFSR instances.
        """
        with open(filename) as f:
            f.readline()                        # skip type tag
            parts = f.readline().split()
            w = int(parts[0])
            r = int(parts[1])
            if w > 32:
                raise ValueError(f"Error: w must be <= 32 bits, got {w}")
            k = w * r
            nbgen = int(f.readline().split()[0])

            generators = []
            for _ in range(nbgen):
                tokens = f.readline().split()
                a_val = int(tokens[0], 16)
                m_val = int(tokens[1])

                # Build a k-bit BitVect with the hex value in the top 32 bits
                # (matching C's vect[0] storage)
                a_bv = BitVect(k, a_val << (k - 32) if k > 32 else a_val)

                gen_L = min(w, L)
                generators.append(cls(w, r, m_val, a_bv, gen_L))

        return generators
