"""
tgfsr.py — Twisted GFSR (TGFSR) generator.

The state consists of r words of w bits (k = w*r total).
Each step applies a twist matrix A to two state words at prescribed
distances, then shifts and XORs the result into the state register.

Characteristic polynomial: computed via Berlekamp-Massey (inherited).

Implements the Generateur abstract base class.
"""


import cython
import sys

from regpoly.bitvect import BitVect
from regpoly.generateur import Generateur


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

    def __init__(self, w: cython.int, r: cython.int, m: cython.int, a: BitVect, L: cython.int) -> None:
        k: cython.int = w * r
        super().__init__(k)
        self.w: cython.int = w
        self.r: cython.int = r
        self.m: cython.int = m
        self.a: BitVect = a.copy()
        self.L = L
        self.smax = sys.maxsize

    # -- Generateur interface ---------------------------------------------

    @classmethod
    def name(cls) -> str:
        return "TGFSR"

    def display(self) -> None:
        a_word: object = self.a._val >> (self.a.n - 32) if self.a.n >= 32 else self.a._val
        print(f" w= {self.w:3d}  r= {self.r:3d}  m= {self.m:3d}  a= {a_word:08x}")

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        self.gen_state = BitVect(self.k, 0)
        self.gen_state.copy_part_from(init_bv, self.k)
        retour: BitVect = BitVect(self.L, 0)
        retour.copy_part_from(init_bv, self.L)
        return retour

    def __next__(self) -> BitVect:
        w: cython.int = self.w
        r: cython.int = self.r
        m: cython.int = self.m
        state: BitVect = self.gen_state

        temp2: BitVect = (state << (w * (r - m - 1))).and_mask(w)
        temp: BitVect = (state << (w * (r - 1))).and_mask(w)

        if temp.get_bit(w - 1) == 0:
            temp >>= 1
            temp = temp2 ^ temp
        else:
            temp >>= 1
            temp = temp2.xor3(temp, self.a)

        state >>= w
        temp = temp.and_mask(w)
        self.gen_state = state ^ temp

        retour: BitVect = BitVect(self.L, 0)
        retour.copy_part_from(self.gen_state, self.L)
        return retour

    def copy(self) -> "TGFSR":
        new: TGFSR = TGFSR(self.w, self.r, self.m, self.a, self.L)
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        return new

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateListFromFile(cls, filename: str, L: cython.int) -> list:
        with open(filename) as f:
            f.readline()
            parts: list = f.readline().split()
            w: cython.int = int(parts[0])
            r: cython.int = int(parts[1])
            if w > 32:
                raise ValueError(f"Error: w must be <= 32 bits, got {w}")
            k: cython.int = w * r
            nbgen: cython.int = int(f.readline().split()[0])

            generators: list = []
            for _ in range(nbgen):
                tokens: list = f.readline().split()
                a_val: object = int(tokens[0], 16)
                m_val: cython.int = int(tokens[1])

                a_bv: BitVect = BitVect(k, a_val << (k - 32) if k > 32 else a_val)

                gen_L: cython.int = min(w, L)
                generators.append(cls(w, r, m_val, a_bv, gen_L))

        return generators
