"""
tausworthe.py — Tausworthe generator (LFSR over GF(2)).

Two variants:
  QuickTaus   — fast iteration; the characteristic polynomial has the
                special form where the step s ranges from 1 to kmq.
  GeneralTaus — general polynomial; slower bit-by-bit recurrence;
                s ranges from 1 to smax (read from the data file).

Characteristic polynomial: x^k + x^Q[NbCoeff-2] + … + x^Q[1] + 1
where Q = [0, Q[1], …, Q[NbCoeff-2], k] in increasing order.

Implements the Generateur abstract base class.
"""


import cython
import sys

from regpoly.bitvect import BitVect
from regpoly.generateur import Generateur

_NLIGNES: cython.int = 256


class Tausworthe(Generateur):
    """
    Tausworthe (LFSR) generator.

    Attributes
    ----------
    k         : int       — degree (inherited)
    Q         : list[int] — polynomial exponents in increasing order:
                            [0, q1, …, qm, k]
    NbCoeff   : int       — len(Q)
    s         : int       — step parameter
    quicktaus : bool      — True = QuickTaus fast variant
    L         : int       — output resolution = state size in bits
    gen_kms   : int       — k − s  (QuickTaus only)
    """

    def __init__(
        self,
        k: cython.int,
        Q: list,
        s: cython.int,
        quicktaus: cython.bint,
        L: cython.int,
    ) -> None:
        super().__init__(k)
        self.Q: list = list(Q)
        self.NbCoeff: cython.int = len(Q)
        self.s: cython.int = s
        self.quicktaus: cython.bint = quicktaus
        self.L = L
        self.smax = sys.maxsize
        self.gen_kms: cython.int = k - s
        self.gen_state = BitVect.zeros(L)

    # -- Generateur interface ---------------------------------------------

    @classmethod
    def name(cls) -> str:
        return "Tausworthe Generator"

    def display(self) -> None:
        parts: list = []
        i: cython.int
        for i in range(self.NbCoeff - 1, 0, -1):
            e: cython.int = self.Q[i]
            if e == 1:
                parts.append(" x +")
            else:
                parts.append(f" x^{e} +")
        output: str = "".join(parts) + " 1 "
        print(f"{output:<40}   s={self.s}")

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        return self._init_general_taus(init_bv)

    def __next__(self) -> BitVect:
        if self.quicktaus:
            return self._iter_quick_taus()
        else:
            return self._iter_general_taus()

    def copy(self) -> "Tausworthe":
        new: Tausworthe = Tausworthe(self.k, self.Q, self.s, self.quicktaus, self.L)
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        new.gen_kms = self.gen_kms
        return new

    def _transition_state(self) -> object:
        """Extract first k bits from the L-bit state (gen_state.n == L >= k)."""
        sv: object = self.gen_state._val
        K: object = self.k
        sn: object = self.gen_state.n
        return (sv >> (sn - K)) & ((1 << K) - 1)

    # -- QuickTaus --------------------------------------------------------

    def _iter_quick_taus(self) -> BitVect:
        """
        QuickTaus: one iteration step using raw int arithmetic.

        Matches C's QuickTaus exactly:
            gen_B = XOR of BVLShift(state, Q[j]) for j = 1..NbCoeff-2
            gen_B ^= state
            gen_B >>= k − s
            state = ANDBVMask(state, k) << s
            state ^= gen_B

        Each BVLShift in C clips to L bits, so we mask after each left shift.
        """
        state_int: object = self.gen_state._val
        L_mask: object = (1 << self.L) - 1
        gen_B_int: object = 0
        j: cython.int
        for j in range(1, self.NbCoeff - 1):
            gen_B_int ^= (state_int << self.Q[j]) & L_mask
        gen_B_int ^= state_int
        gen_B_int >>= self.gen_kms

        self.gen_state.iand_mask(self.k)
        self.gen_state <<= self.s
        self.gen_state._val ^= gen_B_int

        return self.gen_state.copy()

    # -- GeneralTaus ------------------------------------------------------

    def _init_general_taus(self, init: BitVect) -> BitVect:
        self.gen_state = BitVect.zeros(self.L)
        self.gen_state.copy_part_from(init, self.k)
        j: cython.int
        i: cython.int
        bit: cython.int
        for j in range(self.k, self.L):
            bit = 0
            for i in range(self.NbCoeff - 1):
                bit ^= self.gen_state.get_bit(j - (self.k - self.Q[i]))
            if bit:
                self.gen_state.put_bit(j, 1)
        return self.gen_state.copy()

    def _iter_general_taus(self) -> BitVect:
        ss: cython.int = min(self.s, self.L - self.k)
        m: cython.int = self.s
        j: cython.int
        i: cython.int
        bit: cython.int
        while m > 0:
            self.gen_state <<= ss
            for j in range(self.L - ss, self.L):
                bit = 0
                for i in range(self.NbCoeff - 1):
                    bit ^= self.gen_state.get_bit(j - (self.k - self.Q[i]))
                if bit:
                    self.gen_state.put_bit(j, 1)
            m -= ss
            if m < ss:
                ss = m
        return self.gen_state.copy()

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateListFromFile(cls, filename: str, L: cython.int) -> list:
        with open(filename) as f:
            f.readline()
            header: list = f.readline().split()
            n: cython.int = int(header[0])
            quicktaus: cython.bint = bool(int(header[1]))
            file_smax: cython.int = int(header[2]) if not quicktaus else 0

            tokens: list = []
            for line in f:
                tokens.extend(line.split())

        generators: list = []
        idx: cython.int = 0
        e: cython.int
        s: cython.int
        k_poly: cython.int
        kmq: cython.int
        s_max: cython.int
        for _ in range(n):
            R: list = []
            while True:
                e = int(tokens[idx]); idx += 1
                R.append(e)
                if e == 0:
                    break
            Q: list = list(reversed(R))
            nq: cython.int = len(Q)
            k_poly = Q[nq - 1]
            kmq = Q[nq - 1] - Q[nq - 2]
            s_max = kmq if quicktaus else file_smax

            for s in range(1, s_max + 1):
                if k_poly <= _NLIGNES and _shares_small_factor(s, k_poly):
                    continue
                generators.append(cls(k_poly, Q, s, quicktaus, L))

        return generators


def _shares_small_factor(s: int, k: int) -> bool:
    """
    VerifBitsCommuns: return True if s and 2^k − 1 share a common prime factor.

    Uses Python's 3-arg pow (modular exponentiation) — must not be typed
    as cython.int to avoid Cython compiling pow() as C's floating-point pow.
    """
    temp = s
    p = 2
    while p * p <= temp:
        if temp % p == 0:
            if pow(2, k, p) == 1:
                return True
            while temp % p == 0:
                temp //= p
        p += 1
    if temp > 1 and pow(2, k, temp) == 1:
        return True
    return False
