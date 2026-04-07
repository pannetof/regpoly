"""
permutation.py — Bit permutation tempering transformation.

Two modes, selected by the value of p:
  p == 10  : cyclic left rotation of the first w bits by q positions
  otherwise: linear bit permutation  i -> (i*p + q) mod w
             requires gcd(p, w) == 1 so the permutation is invertible

Implements the Transformation abstract base class.
"""

from __future__ import annotations

import copy
import math
import random

from regpoly.bitvect import BitVect
from regpoly.transformation import Transformation


class Permutation(Transformation):
    """
    Bit permutation tempering transformation.

    Attributes
    ----------
    w, w_original : int  — inherited from Transformation
    p, p_original : int  — permutation multiplier (-1 = randomise on update)
    q, q_original : int  — permutation offset    (-1 = randomise on update)
    """

    def __init__(self, w: int, p: int, q: int) -> None:
        """
        InitPermut: initialise a permutation with width w, multiplier p,
        and offset q.

        Raises ValueError if p != -1 and gcd(p, w) != 1.
        """
        super().__init__()
        if p != -1 and p != 10 and math.gcd(p, w) != 1:
            raise ValueError(
                f"Invalid Permutation: gcd({p}, {w}) = {math.gcd(p, w)} != 1"
            )
        self.w: int = w
        self.w_original: int = w
        self.p: int = p
        self.p_original: int = p
        self.q: int = q
        self.q_original: int = q

    # -- Transformation interface -----------------------------------------

    @property
    def display_name(self) -> str:
        return "Permutation"

    def display(self) -> None:
        """DispPermut: print the (p, q) parameters."""
        print(f"Permutation({self.p},{self.q})")

    def __call__(self, A: BitVect) -> BitVect:
        """
        Permut: apply the permutation to the first w bits of A.

        
        remap bit i to position (i*p + q) % w for i in 0..w-1.
        Bits w..n-1 are left unchanged.
        """
       

        n = len(A)
        temp = A.and_invmask(self.w)          # preserve bits w..n-1, zero bits 0..w-1
        permuted = BitVect(n, 0)
        for i in range(self.w):
            u = (i * self.p + self.q) % self.w
            permuted[i] = A[u]
        return temp ^ permuted

    def inverse(self, A: BitVect) -> BitVect:
        """
        InversePermut: invert the permutation on the first w bits of A.

        
        uses the modular inverse of p to reverse the mapping.
        Bits w..n-1 are left unchanged.
        """

        n = len(A)
        inv_p = pow(self.p, -1, self.w)       # modular inverse of p mod w
        temp = A.and_invmask(self.w)
        permuted = BitVect(n, 0)
        for i in range(self.w):
            u = ((i + self.w - self.q) * inv_p) % self.w
            permuted[i] = A[u]
        return temp ^ permuted

    def update_params(self) -> None:
        """
        UpdatePermut: recompute p and q from their originals.

        -1 means choose a new random valid value; otherwise restore the
        stored original.
        """
        if self.p_original == -1:
            p = 0
            while math.gcd(p, self.w) != 1:
                p = random.randrange(1, self.w)
            self.p = p
        else:
            self.p = self.p_original

        if self.q_original == -1:
            self.q = random.randrange(self.w)
        else:
            self.q = self.q_original

    @property
    def type_id(self) -> int:
        """GivePermutID: return PERMID = 3."""
        return 3

    def copy(self) -> "Permutation":
        """CopyPermut: return a deep copy of self."""
        return copy.copy(self)

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def read_params(cls, tokens: list[str]) -> "Permutation":
        """
        Parse a permutation line and return a new Permutation instance.

        File format (all on one line):
            permut <w> <p> <q>

        w  — bit width
        p  — permutation multiplier (-1 = randomise on update)
        q  — permutation offset     (-1 = randomise on update)
        """
        w = int(tokens[1])
        p = int(tokens[2])
        q = int(tokens[3])
        return cls(w, p, q)
