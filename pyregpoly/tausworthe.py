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

from __future__ import annotations

import sys

from bitvect import BitVect
from generateur import Generateur

_NLIGNES = 256   # max degree for which the gcd check is applied


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
        k: int,
        Q: list[int],
        s: int,
        quicktaus: bool,
        L: int,
    ) -> None:
        super().__init__(k)
        self.Q: list[int] = list(Q)
        self.NbCoeff: int = len(Q)
        self.s: int = s
        self.quicktaus: bool = quicktaus
        self.L: int = L
        self.smax: int = sys.maxsize
        self.step: int = s
        self.gen_kms: int = k - s       # used by QuickTaus iteration
        # Override base-class gen_state (k bits) with L-bit state
        self.gen_state = BitVect.zeros(L)

    # -- Generateur interface ---------------------------------------------

    @classmethod
    def name(cls) -> str:
        return "Tausworthe Generator"

    def display(self) -> None:
        """DispTaus: print the polynomial in symbolic form and the step s."""
        parts = []
        for i in range(self.NbCoeff - 1, 0, -1):
            e = self.Q[i]
            if e == 1:
                parts.append(" x +")
            else:
                parts.append(f" x^{e} +")
        output = "".join(parts) + " 1 "
        print(f"{output:<40}   s={self.s}")

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """InitGeneralTaus: initialise state from init_bv.

        The C always assigns InitGeneralTaus to Gen->InitGen, for both
        QuickTaus and GeneralTaus variants.
        """
        return self._init_general_taus(init_bv)

    def __next__(self) -> BitVect:
        """Advance by one step and return the L-bit output word."""
        if self.quicktaus:
            return self._iter_quick_taus()
        else:
            return self._iter_general_taus()

    def copy(self) -> "Tausworthe":
        new = Tausworthe(self.k, self.Q, self.s, self.quicktaus, self.L)
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        new.step = self.step
        new.gen_kms = self.gen_kms
        return new

    # -- QuickTaus --------------------------------------------------------

    def _iter_quick_taus(self) -> BitVect:
        """
        QuickTaus: one iteration step.

            gen_B = XOR of (state << Q[j]) for j = 1..NbCoeff-2, then XOR state
            gen_B >>= k − s
            state  = (state & mask(k)) << s  XOR  gen_B

        When L ≤ 32 (single 32-bit word in C), all arithmetic is done on raw
        Python ints to replicate C's 64-bit ulong overflow behaviour: left-shift
        overflow bits above bit 31 are preserved across iterations and influence
        subsequent mask/shift/XOR operations.  For L > 32, standard BitVect
        operations are used (no overflow issue with multi-word states).
        """
        if self.L <= 32:
            return self._iter_quick_taus_32()
        return self._iter_quick_taus_wide()

    def _iter_quick_taus_32(self) -> BitVect:
        """QuickTaus for L ≤ 32: raw int arithmetic replicating C's 64-bit ulong."""
        _M64 = 0xFFFFFFFFFFFFFFFF               # 64-bit mask (C's ulong width)
        state = self.gen_state._val & _M64       # may carry overflow from prior iter

        # In C, each BVLShift/XOR operates on ulong (64-bit) with wrapping.
        # Python ints are arbitrary precision, so mask to 64 bits at each step.
        gen_B = 0
        for j in range(1, self.NbCoeff - 1):    # j = 1 .. NbCoeff-2
            gen_B = (gen_B ^ ((state << self.Q[j]) & _M64)) & _M64
        gen_B = (gen_B ^ state) & _M64          # XOR state (j=0 term)
        gen_B >>= self.gen_kms                  # right shift by k − s

        # ANDBVMask in C: mask = (ulong)0xFFFFFFFF << (WL - k%WL)
        # where WL=32.  The cast to ulong gives 0x00000000FFFFFFFF, so
        # the shifted mask only covers 32+shift bits, clearing upper overflow.
        c_mask = (0xFFFFFFFF << (32 - (self.k % 32))) & _M64
        state = ((state & c_mask) << self.s) & _M64   # mask then left-shift
        state ^= gen_B                          # XOR gen_B

        self.gen_state._val = state & _M64      # C's ulong is 64-bit
        return self.gen_state.copy()            # copy() masks to L bits

    def _iter_quick_taus_wide(self) -> BitVect:
        """QuickTaus for L > 32: raw ints for gen_B to preserve shift overflow."""
        state_int = self.gen_state._val
        gen_B_int = 0
        for j in range(1, self.NbCoeff - 1):    # j = 1 .. NbCoeff-2
            gen_B_int ^= (state_int << self.Q[j])
        gen_B_int ^= state_int                  # XOR state (j=0 term)
        gen_B_int >>= self.gen_kms              # right shift by k − s
        gen_B_int &= (1 << self.L) - 1          # mask to L bits

        self.gen_state.iand_mask(self.k)        # ANDBVMask(state, state, k)
        self.gen_state <<= self.s               # BVLShiftSelf by s
        self.gen_state._val ^= gen_B_int

        return self.gen_state.copy()

    # -- GeneralTaus ------------------------------------------------------

    def _init_general_taus(self, init: BitVect) -> BitVect:
        """
        InitGeneralTaus: copy first k bits from init, then extend the state
        to L bits using the linear recurrence.
        """
        self.gen_state = BitVect.zeros(self.L)
        self.gen_state.copy_part_from(init, self.k)
        for j in range(self.k, self.L):
            bit = 0
            for i in range(self.NbCoeff - 1):
                bit ^= self.gen_state.get_bit(j - (self.k - self.Q[i]))
            if bit:
                self.gen_state.put_bit(j, 1)
        return self.gen_state.copy()

    def _iter_general_taus(self) -> BitVect:
        """
        GeneralTaus: one iteration step, producing s new bits via the
        linear recurrence, ss bits at a time (ss = min(s, L − k)).
        """
        ss = min(self.s, self.L - self.k)
        m = self.s
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
    def CreateListFromFile(cls, filename: str, L: int) -> list["Tausworthe"]:
        """
        ReadDataTaus: read generator parameters from filename and return
        a list of Tausworthe instances.

        File format
        -----------
        Line 1 : type tag  (e.g. "taus2") — skipped
        Line 2 : n  quicktaus  [smax]
                   n         = number of polynomials
                   quicktaus = 1 for QuickTaus, 0 for GeneralTaus
                   smax      = max step (GeneralTaus only, must be < 256)
        Lines 3 .. n+2:
            k  e_{m-1}  …  e_1  0
            Exponents in strictly decreasing order, terminated by 0.

        For QuickTaus:   s = 1 … kmq  (kmq = k − e_{m-1})
        For GeneralTaus: s = 1 … smax
        Both: skip s when s and 2^k − 1 share a common prime factor
              (only checked when k ≤ 256, matching the C table limit).
        """
        with open(filename) as f:
            f.readline()                    # skip type tag
            header = f.readline().split()
            n          = int(header[0])
            quicktaus  = bool(int(header[1]))
            file_smax  = int(header[2]) if not quicktaus else None

            tokens = []
            for line in f:
                tokens.extend(line.split())

        generators = []
        idx = 0
        for _ in range(n):
            # Read exponents in decreasing order, terminated by 0
            R = []
            while True:
                e = int(tokens[idx]); idx += 1
                R.append(e)
                if e == 0:
                    break
            # R = [k, e_{m-1}, …, e_1, 0]  →  Q = [0, e_1, …, e_{m-1}, k]
            Q = list(reversed(R))
            k_poly = Q[-1]
            kmq    = Q[-1] - Q[-2]          # k − Q[NbCoeff-2]
            s_max  = kmq if quicktaus else file_smax

            for s in range(1, s_max + 1):
                if k_poly <= _NLIGNES and _shares_small_factor(s, k_poly):
                    continue
                generators.append(cls(k_poly, Q, s, quicktaus, L))

        return generators


def _shares_small_factor(s: int, k: int) -> bool:
    """
    VerifBitsCommuns: return True if s and 2^k − 1 share a common prime factor.

    Factors s into primes; for each prime p, checks p | 2^k − 1  via
    pow(2, k, p) == 1.  Equivalent to the C table-based check for s < 256.
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
