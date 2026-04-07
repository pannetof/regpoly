"""
temper_mk.py — Matsumoto-Kurita / Matsumoto-Nishimura tempering transformation.

Two types, selected at construction:
  type 1  Matsumoto-Kurita (type I):
              y  = x ^ ((x << eta) & b)
              y ^= (y << mu)  & c
  type 2  Matsumoto-Nishimura (type II):
              y  = x ^ (x >> u)
              y ^= (y << eta) & b
              y ^= (y << mu)  & c
              y ^= (y >> l)

All operations act on the first w public bits of the state BitVect;
bits w..n-1 are left unchanged.

Inverse is only implemented for type I.

Implements the Transformation abstract base class.
"""

from __future__ import annotations

import copy
import random as _random

from bitvect import BitVect
from transformation import Transformation

TEMPMKID = 1


class TemperMK(Transformation):
    """
    Matsumoto-Kurita / Matsumoto-Nishimura tempering transformation.

    Attributes (in addition to w, w_original from Transformation)
    ------
    type         : int     — 1 = MK type I, 2 = MN type II
    eta          : int     — left-shift amount for the b mask step
    mu           : int     — left-shift amount for the c mask step
    u            : int     — right-shift amount (type II only)
    l            : int     — final right-shift amount (type II only)
    random       : int     — -1 = randomise b/c on update, else restore originals
    optimize     : bool    — True = allow OptimizeTemper to modify b/c
    disp_progress: bool    — True = print progress during optimisation
    limit_v      : int     — optimisation resolution limit
    b, c         : BitVect — current mask vectors (w bits)
    b_original   : BitVect — initial b (restored when random != -1)
    c_original   : BitVect — initial c (restored when random != -1)
    best_b       : BitVect — best b found during optimisation
    best_c       : BitVect — best c found during optimisation
    """

    def __init__(
        self,
        w: int,
        eta: int,
        mu: int,
        u: int,
        l: int,
        b: BitVect,
        c: BitVect,
        *,
        random: int = 0,
        optimize: bool = False,
        type: int = 1,
        disp_progress: bool = False,
        limit_v: int = 0,
    ) -> None:
        """
        InitTemperMK: initialise the transformation.

        Parameters
        ----------
        w            : output bit width
        eta, mu      : shift amounts for b and c mask steps
        u, l         : shift amounts for the outer steps (type II only)
        b, c         : initial mask bit vectors (must be w bits wide)
        random       : -1 → randomise b/c on each update_params() call
        optimize     : enable the OptimizeTemper search
        type         : 1 = MK (type I), 2 = MN (type II)
        disp_progress: print progress during optimisation
        limit_v      : resolution limit for optimisation
        """
        super().__init__()
        if len(b) != w or len(c) != w:
            raise ValueError(f"b and c must be {w}-bit vectors")
        self.w: int = w
        self.w_original: int = w
        self.type: int = type
        self.eta: int = eta
        self.mu: int = mu
        self.u: int = u
        self.l: int = l
        self.random: int = random
        self.optimize: bool = optimize
        self.disp_progress: bool = disp_progress
        self.limit_v: int = limit_v
        self.b: BitVect = b.copy()
        self.b_original: BitVect = b.copy()
        self.best_b: BitVect = b.copy()
        self.c: BitVect = c.copy()
        self.c_original: BitVect = c.copy()
        self.best_c: BitVect = c.copy()

    # -- Helpers ----------------------------------------------------------

    def _pad(self, bv: BitVect, n: int) -> BitVect:
        """Zero-pad a w-bit vector to n bits, keeping bit 0 at position 0."""
        return BitVect(n, bv._val << (n - self.w))

    # -- Transformation interface -----------------------------------------

    @property
    def display_name(self) -> str:
        return ("Matsumoto-Kurita Tempering(II)"
                if self.type == 2
                else "Matsumoto-Kurita(I) Tempering")

    def display(self) -> None:
        """DispTemperMK: print type, shift parameters and best mask vectors."""
        hex_digits = (self.w + 3) // 4
        if self.type == 1:
            print(
                f"Matsumoto-Kurita Tempering (I)"
                f"({self.eta},{self.mu})(w={self.w})",
                end="",
            )
        else:
            print(
                f"Matsumoto-Nishimura Tempering (II)"
                f"({self.u},{self.eta},{self.mu},{self.l})(w={self.w})",
                end="",
            )
        print(
            f"  b = {self.best_b._val:0{hex_digits}x}"
            f"  c = {self.best_c._val:0{hex_digits}x}"
        )

    def __call__(self, A: BitVect) -> BitVect:
        """
        TemperingMK: apply the tempering to the first w bits of A.

        Returns a new BitVect with bits 0..w-1 tempered and w..n-1 unchanged.
        """
        n = len(A)
        state_w = A.and_mask(self.w)
        rest = A.and_invmask(self.w)

        b_n = self._pad(self.b, n)
        c_n = self._pad(self.c, n)

        if self.type == 2:
            state_w = (state_w ^ (state_w >> self.u)).and_mask(self.w)

        state_w = state_w ^ ((state_w << self.eta) & b_n)
        state_w = state_w ^ ((state_w << self.mu) & c_n)

        if self.type == 2:
            state_w = (state_w ^ (state_w >> self.l)).and_mask(self.w)

        return rest ^ state_w

    def inverse(self, A: BitVect) -> BitVect:
        """
        InverseTemperMK: invert the type I tempering on the first w bits of A.

        Only implemented for type I; raises NotImplementedError for type II.
        """
        if self.type == 2:
            raise NotImplementedError("TemperMK inverse is only implemented for type I")

        n = len(A)
        rest = A.and_invmask(self.w)
        state_w = A.and_mask(self.w)

        c_n = self._pad(self.c, n)
        b_n = self._pad(self.b, n)

        # Invert  y = x ^ ((x << mu) & c)
        mask = c_n.copy()
        inv1 = state_w.copy()
        i = 1
        while i * self.mu <= self.w:
            temp = (state_w << (i * self.mu)) & mask
            inv1 = inv1 ^ temp
            mask = mask & (c_n << (i * self.mu))
            i += 1

        # Invert  y = x ^ ((x << eta) & b)
        mask = b_n.copy()
        inv2 = inv1.copy()
        i = 1
        while i * self.eta < self.w:
            temp = (inv1 << (i * self.eta)) & mask
            inv2 = inv2 ^ temp
            mask = mask & (b_n << (i * self.eta))
            i += 1

        return rest ^ inv2

    def update_params(self) -> None:
        """
        UpdateTemperMK: randomise b and c (if random == -1) or restore originals.

        Masks both to w bits and updates best_b / best_c to match.
        """
        if self.random == -1:
            self.b = BitVect.random(self.w).and_mask(self.w)
            self.c = BitVect.random(self.w).and_mask(self.w)
        else:
            self.b = self.b_original.copy().and_mask(self.w)
            self.c = self.c_original.copy().and_mask(self.w)
        self.best_b = self.b.copy()
        self.best_c = self.c.copy()

    @property
    def type_id(self) -> int:
        """GiveTemperMKID: return TEMPMKID = 1."""
        return TEMPMKID

    def copy(self) -> "TemperMK":
        """CopyTemperMK: return a deep copy of self."""
        return copy.deepcopy(self)

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def _extra_tokens(cls, header_tokens: list[str], f) -> list[str]:
        """
        Read the b and c lines from the file when nb_words != -1.

        Matches C's ReadLn(f) + fscanf loop pattern: b is on the next line,
        c is on the line after that.
        """
        type_str = header_tokens[0]
        idx = 4
        if "2" in type_str:
            idx += 2   # skip u and l
        if int(header_tokens[idx]) == -1:
            return []
        return f.readline().split() + f.readline().split()

    @classmethod
    def read_params(cls, tokens: list[str]) -> "TemperMK":
        """
        Parse a TemperMK token list and return a new instance.

        tokens[0] is the type string; all data (including b/c hex words
        appended by _extra_tokens) is already present in the flat list.

        Token layout:
            [0]  type_str   {tempMK|tempMKopt|tempMK2|tempMK2opt}
            [1]  w
            [2]  eta
            [3]  mu
            [4]  u          (type II only)
            [5]  l          (type II only)
         [4|6]  nb_words   (-1 = randomise)
         [5|7]  disp        (opt variants only)
         [6|8]  limit_v     (opt variants only)
          then  nb_words hex tokens for b
          then  nb_words hex tokens for c   (both absent when nb_words == -1)
        """
        import sys as _sys

        type_str  = tokens[0]
        w         = int(tokens[1])
        eta       = int(tokens[2])
        mu        = int(tokens[3])
        idx       = 4

        if "2" in type_str:
            type_ = 2
            u = int(tokens[idx]); idx += 1
            l = int(tokens[idx]); idx += 1
        else:
            type_ = 1
            u = l = 0

        nb_words = int(tokens[idx]); idx += 1

        optimize      = False
        disp_progress = False
        limit_v       = _sys.maxsize
        if "opt" in type_str:
            optimize      = True
            disp_progress = bool(int(tokens[idx])); idx += 1
            limit_v       = int(tokens[idx]);        idx += 1

        if nb_words != -1:
            b = cls._read_bitvect(tokens[idx:idx + nb_words], nb_words, w)
            idx += nb_words
            c = cls._read_bitvect(tokens[idx:idx + nb_words], nb_words, w)
        else:
            b = BitVect.zeros(w)
            c = BitVect.zeros(w)

        return cls(
            w, eta, mu, u, l, b, c,
            random=nb_words,
            optimize=optimize,
            type=type_,
            disp_progress=disp_progress,
            limit_v=limit_v,
        )

    @staticmethod
    def _read_bitvect(hex_tokens: list[str], nb_words: int, w: int) -> BitVect:
        """
        Reconstruct a w-bit BitVect from nb_words 32-bit big-endian hex tokens.

        In the C storage, vect[0] holds the most significant 32 bits (public
        bits 0..31), vect[1] the next 32 bits, and so on.  When w is not a
        multiple of 32, the last word's padding bits are in its low-order
        positions and must be discarded.
        """
        val = 0
        for token in hex_tokens[:nb_words]:
            val = (val << 32) | int(token, 16)
        # Discard padding bits added by the last (possibly partial) word.
        shift = 32 * nb_words - w
        val >>= shift
        return BitVect(w, val)

    # -- Optimisation (depends on mecf / polynomials / myzen) -------------

    def optimize_temper(self, combinaison, params) -> None:
        """
        TestMETemperMK: optimise b and c to minimise the equidistribution gap.

        Not yet implemented — requires mecf, polynomials, and myzen modules.
        """
        raise NotImplementedError(
            "optimize_temper requires mecf, polynomials, and myzen — not yet translated"
        )
