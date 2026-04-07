"""
polylcg.py — Polynomial LCG generator over GF(2).

At each step the k-bit state is shifted left by one position; if the outgoing
MSB (public bit 0) was 1 the polynomial is XORed in, keeping the state within
the k-bit window.  The first L public bits of the updated state are returned as
the output word.

Implements the Generateur abstract base class.
"""

from __future__ import annotations

from bitvect import BitVect
from generateur import Generateur


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

    def __init__(self, k: int, poly: BitVect, L: int) -> None:
        """
        InitParamPolyLCG: initialise the generator.

        Parameters
        ----------
        k    : degree / state width in bits
        poly : k-bit polynomial BitVect (leading x^k term implicit)
        L    : output resolution in bits (L ≤ k)
        """
        super().__init__(k)
        self.poly: BitVect = poly.copy()
        self.L: int = L

    # -- Generateur interface ---------------------------------------------

    def display(self) -> None:
        """
        DispPolyLCG: print the polynomial in exponent and hexadecimal notation.

        Exponents are listed in decreasing order (degree first, constant last).
        """
        exponents = [self.k - i - 1 for i in range(self.k)
                     if self.poly.get_bit(i) == 1]
        print(f" {self.k} " + " ".join(str(e) for e in exponents) + " ")
        hex_digits = (self.k + 3) // 4
        print(f"hexadecimal notation:\n {self.poly._val:0{hex_digits}x} ")

    @classmethod
    def name(cls) -> str:
        """Return the generator type name."""
        return "Polynomial LCG"

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """
        InitPolyLCG: set the state to init_bv and return the first L output bits.
        """
        self.gen_state = init_bv.copy()
        return self.gen_state.and_mask(self.L)

    def __next__(self) -> BitVect:
        """
        PolyLCG: advance by one step and return the first L bits of the new state.

        Algorithm:
            xor  ← MSB of state (public bit 0)
            state <<= 1          (shift left, LSB fills with 0)
            if xor: state ^= poly
        """
        xor = self.gen_state.get_bit(0) == 1
        self.gen_state <<= 1
        if xor:
            self.gen_state ^= self.poly
        return self.gen_state.and_mask(self.L)

    def copy(self) -> "PolyLCG":
        """CopyPolyLCG: return a deep copy of self."""
        new = PolyLCG(self.k, self.poly, self.L)
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        new.step = self.step
        return new

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateListFromFile(cls, filename: str, L: int) -> list["PolyLCG"]:
        """
        ReadDataPolyLCG: read generator parameters from filename and return
        a list of PolyLCG instances.

        File format
        -----------
        Line 1 : type tag  ("polylcg")  — skipped
        Line 2 : n  (number of generators)
        Lines 3..n+2 :
            k  e1  e2  …  0
            First token is the degree k.  Subsequent tokens are non-leading
            exponents in decreasing order, terminated by 0 (which also
            represents the constant term x^0 and is always present in a
            primitive polynomial).

        The loop mirrors the C fscanf logic exactly: each exponent e (including
        the terminating 0) sets public bit (k − e − 1) of the polynomial
        BitVect, then stops when e == 0.
        """
        with open(filename) as f:
            f.readline()                     # skip type tag
            n = int(f.readline())

            # Collect remaining content as a flat token stream so that
            # exponents split across lines are handled transparently.
            tokens = []
            for line in f:
                tokens.extend(line.split())

        generators = []
        idx = 0
        for _ in range(n):
            k = int(tokens[idx]); idx += 1
            poly = BitVect.zeros(k)
            while True:
                e = int(tokens[idx]); idx += 1
                poly.put_bit(k - e - 1, 1)  # set the bit for exponent e
                if e == 0:                   # 0 = constant term, also terminates
                    break
            generators.append(cls(k, poly, L))
        return generators
