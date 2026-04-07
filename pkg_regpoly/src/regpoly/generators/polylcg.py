"""
polylcg.py — Polynomial LCG generator over GF(2).

At each step the k-bit state is shifted left by one position; if the outgoing
MSB (public bit 0) was 1 the polynomial is XORed in, keeping the state within
the k-bit window.  The first L public bits of the updated state are returned as
the output word.

Implements the Generateur abstract base class.
"""


import cython

from regpoly.bitvect import BitVect
from regpoly.generateur import Generateur


class PolyLCG(Generateur):
    """
    Polynomial LCG generator (LFSR over GF(2)).

    Attributes
    ----------
    k    : int     — degree; also the state width in bits (inherited)
    L    : int     — output resolution; first L bits of state are returned
    poly : BitVect — k-bit characteristic polynomial, WITHOUT the leading x^k
                     term.  Public bit i set ↔ x^(k-1-i) term present.
    """

    def __init__(self, k: cython.int, poly: BitVect, L: cython.int) -> None:
        super().__init__(k)
        self.poly: BitVect = poly.copy()
        self.L = L

    # -- Generateur interface ---------------------------------------------

    def display(self) -> None:
        exponents = [self.k - i - 1 for i in range(self.k)
                     if self.poly.get_bit(i) == 1]
        print(f" {self.k} " + " ".join(str(e) for e in exponents) + " ")
        hex_digits: cython.int = (self.k + 3) // 4
        print(f"hexadecimal notation:\n {self.poly._val:0{hex_digits}x} ")

    @classmethod
    def name(cls) -> str:
        return "Polynomial LCG"

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        self.gen_state = init_bv.copy()
        return self.gen_state.and_mask(self.L)

    def __next__(self) -> BitVect:
        xor: cython.bint = self.gen_state.get_bit(0) == 1
        self.gen_state <<= 1
        if xor:
            self.gen_state ^= self.poly
        return self.gen_state.and_mask(self.L)

    def copy(self) -> "PolyLCG":
        new: PolyLCG = PolyLCG(self.k, self.poly, self.L)
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        return new

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateListFromFile(cls, filename: str, L: cython.int) -> list:
        with open(filename) as f:
            f.readline()
            n: cython.int = int(f.readline())
            tokens: list = []
            for line in f:
                tokens.extend(line.split())

        generators: list = []
        idx: cython.int = 0
        k: cython.int
        e: cython.int
        for _ in range(n):
            k = int(tokens[idx]); idx += 1
            poly: BitVect = BitVect.zeros(k)
            while True:
                e = int(tokens[idx]); idx += 1
                poly.put_bit(k - e - 1, 1)
                if e == 0:
                    break
            generators.append(cls(k, poly, L))
        return generators
