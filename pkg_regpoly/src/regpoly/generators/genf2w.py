"""
genf2w.py — Generators over the extension field GF(2^w).

Two variants:
  GenF2wLFSR    — LFSR in GF(2^w)
  GenF2wPolyLCG — Polynomial LCG in GF(2^w)[z]/P(z)

Both use configurable coefficients in GF(2^w) with either polynomial-basis
or normal-basis arithmetic, and a shared file format.

GF(2^w) element arithmetic uses cython.ulonglong (C unsigned long long,
64-bit) for native register-level speed.  Supports w up to 64.

Implements the Generateur abstract base class.
"""

import cython
import sys

from regpoly.bitvect import BitVect
from regpoly.generateur import Generateur

_ALL_ONES_64: cython.ulonglong = ~cython.cast(cython.ulonglong, 0)


# ── GF(2^w) arithmetic (native 64-bit) ─────────────────────────────────

def _multiply_z(a: cython.ulonglong, k: cython.int,
                modM: cython.ulonglong) -> cython.ulonglong:
    """Multiply element a by z^k in GF(2^w) modulo modM (polynomial basis)."""
    i: cython.int
    for i in range(k):
        if a & cython.cast(cython.ulonglong, 1):
            a = (a >> 1) ^ modM
        else:
            a >>= 1
    return a


def _multiply_poly_basis(a: cython.ulonglong, b: cython.ulonglong,
                         w: cython.int,
                         modM: cython.ulonglong) -> cython.ulonglong:
    """Multiply a * b in GF(2^w) using polynomial basis."""
    res: cython.ulonglong = 0
    verif: cython.ulonglong = 1
    i: cython.int
    for i in range(w):
        if b & verif:
            res ^= _multiply_z(a, w - i - 1, modM)
        verif <<= 1
    return res


def _multiply_normal_basis(a: cython.ulonglong, b: cython.ulonglong,
                           w: cython.int,
                           table: list) -> cython.ulonglong:
    """Multiply a * b in GF(2^w) using normal basis with precomputed table."""
    if a == 0 or b == 0:
        return 0
    res: cython.ulonglong = 0
    verifa: cython.ulonglong = cython.cast(cython.ulonglong, 1) << (w - 1)
    j: cython.int
    i: cython.int
    for j in range(w):
        verifb: cython.ulonglong = cython.cast(cython.ulonglong, 1) << (w - 1)
        for i in range(w):
            if (verifa & a) and (verifb & b):
                res ^= cython.cast(cython.ulonglong, table[i * w + j])
            verifb >>= 1
        verifa >>= 1
    return res


def _make_table(w: cython.int, modM_val: cython.ulonglong) -> list:
    """
    Build the multiplication table for normal basis in GF(2^w).

    Returns a list of w*w entries where table[i*w + j] holds the product
    of the i-th and j-th normal basis elements, expressed in the normal
    basis representation.
    """
    one: cython.ulonglong = 1
    # x[j] = x^j in polynomial basis
    x: list = [cython.cast(cython.ulonglong, 0)] * (2 * w)
    x[0] = one << (w - 1)  # x^0 = 1
    j: cython.int
    xj: cython.ulonglong
    for j in range(2 * w - 1):
        xj = cython.cast(cython.ulonglong, x[j])
        if xj & one:
            x[j + 1] = (xj >> 1) ^ modM_val
        else:
            x[j + 1] = xj >> 1

    # xx[j] = x^(2^j) (Frobenius elements)
    xx: list = [cython.cast(cython.ulonglong, 0)] * w
    xx[0] = x[1]  # x^1
    i: cython.int
    xxj: cython.ulonglong
    for j in range(w - 1):
        verif: cython.ulonglong = one << (w - 1)
        xxj = cython.cast(cython.ulonglong, xx[j])
        val: cython.ulonglong = 0
        for i in range(w):
            if xxj & verif:
                val ^= cython.cast(cython.ulonglong, x[2 * i])
            verif >>= 1
        xx[j + 1] = val

    # Build matrix C: column j = polynomial repr of xx[j]
    # C is a w×w bit matrix stored as list of w ulonglong (rows)
    C: list = [cython.cast(cython.ulonglong, 0)] * w
    for j in range(w):
        verif = one << (w - 1)
        xxj = cython.cast(cython.ulonglong, xx[j])
        for i in range(w):
            if xxj & verif:
                C[i] = cython.cast(cython.ulonglong, C[i]) | (one << (w - 1 - j))
            verif >>= 1

    # Invert C over GF(2) via Gauss-Jordan
    InvC: list = [cython.cast(cython.ulonglong, 0)] * w
    for i in range(w):
        InvC[i] = one << (w - 1 - i)  # identity matrix

    work: list = list(C)
    col: cython.int
    for col in range(w):
        bit: cython.ulonglong = one << (w - 1 - col)
        # Find pivot
        pivot: cython.int = -1
        for i in range(col, w):
            if cython.cast(cython.ulonglong, work[i]) & bit:
                pivot = i
                break
        if pivot < 0:
            raise ValueError("No normal basis for this modM")
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            InvC[col], InvC[pivot] = InvC[pivot], InvC[col]
        wc: cython.ulonglong = cython.cast(cython.ulonglong, work[col])
        ic: cython.ulonglong = cython.cast(cython.ulonglong, InvC[col])
        for i in range(w):
            if i != col and (cython.cast(cython.ulonglong, work[i]) & bit):
                work[i] = cython.cast(cython.ulonglong, work[i]) ^ wc
                InvC[i] = cython.cast(cython.ulonglong, InvC[i]) ^ ic

    # Transpose InvC → TransInvC
    TransInvC: list = [cython.cast(cython.ulonglong, 0)] * w
    for i in range(w):
        inv_i: cython.ulonglong = cython.cast(cython.ulonglong, InvC[i])
        for j in range(w):
            if inv_i & (one << (w - 1 - j)):
                TransInvC[j] = cython.cast(cython.ulonglong, TransInvC[j]) | (one << (w - 1 - i))

    # Build multiplication table
    table: list = [cython.cast(cython.ulonglong, 0)] * (w * w)
    jj: cython.int
    for i in range(w):
        xxi: cython.ulonglong = cython.cast(cython.ulonglong, xx[i])
        for j in range(w):
            temp: cython.ulonglong = _multiply_poly_basis(
                xxi, cython.cast(cython.ulonglong, xx[j]), w, modM_val)
            entry: cython.ulonglong = 0
            for jj in range(w):
                if (temp >> (w - 1 - jj)) & one:
                    entry ^= cython.cast(cython.ulonglong, TransInvC[jj])
            table[i * w + j] = entry

    return table


# ── Base class ──────────────────────────────────────────────────────────

class GenF2w(Generateur):
    """
    Common base for LFSR and PolyLCG generators over GF(2^w).

    Attributes
    ----------
    k          : int      — degree = w * r  (inherited)
    w          : int      — extension degree (bits per element)
    r          : int      — number of register stages
    nbcoeff    : int      — number of recurrence coefficients
    nocoeff    : list[int] — coefficient positions
    coeff      : list[int] — coefficient values in GF(2^w)
    modM       : int      — irreducible polynomial (polynomial basis modulus)
    normal_basis : bool   — True = use normal basis, False = polynomial basis
    table      : list     — multiplication table (normal basis only)
    maskw      : int      — w-bit mask
    _step_count: int      — number of iterations per call
    L          : int      — output resolution in bits
    """

    def __init__(
        self,
        w: cython.int,
        r: cython.int,
        nbcoeff: cython.int,
        nocoeff: list,
        coeff: list,
        modM: object,
        normal_basis: cython.bint,
        step_count: cython.int,
        L: cython.int,
    ) -> None:
        k: cython.int = w * r
        super().__init__(k)
        self.w: cython.int = w
        self.r: cython.int = r
        self.nbcoeff: cython.int = nbcoeff
        self.nocoeff: list = list(nocoeff)
        self.coeff: list = list(coeff)
        self.modM: cython.ulonglong = cython.cast(cython.ulonglong, modM)
        self.normal_basis: cython.bint = normal_basis
        self.maskw: cython.ulonglong = _ALL_ONES_64 >> (64 - w)
        self._step_count: cython.int = step_count
        self.smax = sys.maxsize
        self.L = min(L, k)

        if normal_basis:
            self.table: list = _make_table(w, self.modM)
        else:
            self.table = []

    def _multiply(self, a: cython.ulonglong,
                  b: cython.ulonglong) -> cython.ulonglong:
        """Multiply a * b in GF(2^w)."""
        if self.normal_basis:
            return _multiply_normal_basis(a, b, self.w, self.table)
        else:
            return _multiply_poly_basis(a, b, self.w, self.modM)

    # -- Word access (V / SetV from left) ---------------------------------

    def _V(self, idx: cython.int) -> cython.ulonglong:
        """Return the idx-th w-bit element from the left of gen_state."""
        shift: object = self.gen_state.n - (idx + 1) * self.w
        # Mask to w bits first (as Python object), then cast to ulonglong
        return cython.cast(cython.ulonglong,
                           (self.gen_state._val >> shift) & cython.cast(object, self.maskw))

    def _SetV(self, idx: cython.int, val: cython.ulonglong) -> None:
        """Set the idx-th w-bit element from the left of gen_state."""
        shift: object = self.gen_state.n - (idx + 1) * self.w
        maskw_obj: object = self.maskw
        mask: object = maskw_obj << shift
        self.gen_state._val = (self.gen_state._val & ~mask) | (
            (cython.cast(object, val) & maskw_obj) << shift)

    # -- Display ----------------------------------------------------------

    def display(self) -> None:
        """DispGenF2w: print generator type, recurrence polynomial, and step."""
        print(self._type_name())
        parts: list = []
        j: cython.int
        for j in range(self.nbcoeff - 1):
            parts.append(f"({self.coeff[j]:08x})z^{self.nocoeff[j]}")
        last: cython.int = self.nbcoeff - 1
        if self.nocoeff[last] == 0:
            parts.append(f"({self.coeff[last]:08x})")
        else:
            parts.append(f"({self.coeff[last]:08x})z^{self.nocoeff[last]}")
        print(f"z^{self.r} + " + " +".join(parts)
              + f" -- coeffs in F_2/P(z) where P(z)={self.modM:08x}")
        print(f"Step = {self._step_count}")

    @classmethod
    def name(cls) -> str:
        return "Generator in F_{2^w}"

    @classmethod
    def display_name(cls) -> str:
        return "Generator that uses F_{2^w}"

    @classmethod
    def _type_name(cls) -> str:
        raise NotImplementedError

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateListFromFile(cls, filename: str, L: cython.int) -> list:
        """
        ReadDataGenF2w: read genf2w parameters from filename.

        File format
        -----------
        Line 1 : type tag ("genf2w")
        Line 2 : type (1=LFSR, 0=PolyLCG)
        Line 3 : normalbasis (0 or 1)
        Line 4 : nbgen
        Per generator:
            w r modM(hex) step nbcoeff
            coeff1(hex) nocoeff1  coeff2(hex) nocoeff2 ...
        """
        import re
        _token_re = re.compile(r'\S+')

        with open(filename) as f:
            f.readline()                                    # skip tag
            gen_type: cython.int = int(f.readline().split()[0])
            normal_basis: cython.bint = bool(int(f.readline().split()[0]))
            nbgen: cython.int = int(f.readline().split()[0])

            # Token iterator: reads lines on demand, not the whole file
            def _tokens():
                for line in f:
                    for m in _token_re.finditer(line):
                        yield m.group()

            it = _tokens()
            generators: list = []
            for _ in range(nbgen):
                w: cython.int = int(next(it))
                r: cython.int = int(next(it))
                modM: object = int(next(it), 16)
                step_count: cython.int = int(next(it))
                nbcoeff: cython.int = int(next(it))

                coeff: list = []
                nocoeff: list = []
                j: cython.int
                for j in range(nbcoeff):
                    c: object = int(next(it), 16)
                    n: cython.int = int(next(it))
                    coeff.append(c)
                    nocoeff.append(n)

                if gen_type == 1:  # LFSR
                    gen = GenF2wLFSR(w, r, nbcoeff, nocoeff, coeff,
                                     modM, normal_basis, step_count, L)
                else:              # PolyLCG
                    gen = GenF2wPolyLCG(w, r, nbcoeff, nocoeff, coeff,
                                        modM, normal_basis, step_count, L)
                generators.append(gen)

        return generators


# ── LFSR in GF(2^w) ────────────────────────────────────────────────────

class GenF2wLFSR(GenF2w):
    """
    LFSR in GF(2^w).

    Recurrence: X_n = sum_{j} coeff[j] * X_{n-r+nocoeff[j]}
    The state stores r elements of w bits (leftmost = oldest output).
    Each step: left-shift by w, compute feedback, place at position r-1.
    """

    def __init__(self, w, r, nbcoeff, nocoeff, coeff, modM,
                 normal_basis, step_count, L) -> None:
        super().__init__(w, r, nbcoeff, nocoeff, coeff, modM,
                         normal_basis, step_count, L)
        self._p: cython.int = r
        k_bits: cython.int = w * r
        if k_bits < 32:
            self.L = ((L - k_bits) // w) * w + k_bits
        self._state_bits: cython.int = max(self.L, k_bits)

    @classmethod
    def _type_name(cls) -> str:
        return "LFSR in F_{2^w}"

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """InitGenF2wLFSR: copy init, pre-compute extended output, return L bits."""
        self.gen_state = BitVect(self._state_bits, 0)
        self.gen_state.copy_part_from(init_bv, self.k)
        self._p = self.r

        p: cython.int = self.L // self.w - self.r
        m: cython.int
        j: cython.int
        for m in range(p):
            res: cython.ulonglong = 0
            for j in range(self.nbcoeff):
                res ^= self._multiply(
                    self._V(self.nocoeff[j] + m),
                    cython.cast(cython.ulonglong, self.coeff[j]))
            self._SetV(self.r + m, res)

        sb: int = self._state_bits
        L: int = self.L
        if L <= sb:
            result: object = (self.gen_state._val >> (sb - L)) & ((1 << L) - 1)
        else:
            result = self.gen_state._val << (L - sb)
        return BitVect(L, result)

    def __next__(self) -> BitVect:
        """GenF2wLFSR: STEP iterations of the LFSR recurrence."""
        w: cython.int = self.w
        r: cython.int = self.r
        i: cython.int
        j: cython.int

        for i in range(self._step_count):
            res: cython.ulonglong = 0
            for j in range(self.nbcoeff):
                res ^= self._multiply(
                    self._V(self.nocoeff[j]),
                    cython.cast(cython.ulonglong, self.coeff[j]))
            self.gen_state <<= w
            self._SetV(r - 1, res)

        p: cython.int = self.L // w - r
        m: cython.int
        for m in range(p):
            res = 0
            for j in range(self.nbcoeff):
                res ^= self._multiply(
                    self._V(self.nocoeff[j] + m),
                    cython.cast(cython.ulonglong, self.coeff[j]))
            self._SetV(r + m, res)

        sb: int = self._state_bits
        L: int = self.L
        if L <= sb:
            result: object = (self.gen_state._val >> (sb - L)) & ((1 << L) - 1)
        else:
            result = self.gen_state._val << (L - sb)
        return BitVect(L, result)

    def _transition_state(self) -> object:
        """Extract first k bits from the state (may be larger than k)."""
        sv: object = self.gen_state._val
        sn: object = self.gen_state.n
        K: object = self.k
        if sn == K:
            return sv
        return (sv >> (sn - K)) & ((1 << K) - 1)

    def copy(self) -> "GenF2wLFSR":
        new: GenF2wLFSR = GenF2wLFSR(
            self.w, self.r, self.nbcoeff, self.nocoeff, self.coeff,
            self.modM, self.normal_basis, self._step_count, self.L,
        )
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        new._p = self._p
        new._state_bits = self._state_bits
        new.table = self.table
        return new


# ── Polynomial LCG in GF(2^w)[z]/P(z) ──────────────────────────────────

class GenF2wPolyLCG(GenF2w):
    """
    Polynomial LCG in GF(2^w)[z]/P(z).

    Each step: save last element VP, right-shift state by w bits,
    then for each coefficient, XOR multiply(VP, coeff) into the
    corresponding position.
    """

    @classmethod
    def _type_name(cls) -> str:
        return "Polynomial LCG in F_{2^w}[z]/P(z)"

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """InitGenF2w: copy init into state, return first L bits."""
        self.gen_state = BitVect(self.k, 0)
        self.gen_state.copy_part_from(init_bv, self.k)

        L: int = self.L
        result: object = (self.gen_state._val >> (self.k - L)) & ((1 << L) - 1)
        return BitVect(L, result)

    def __next__(self) -> BitVect:
        """GenF2wPolyLCG: STEP iterations of the polynomial LCG."""
        w: cython.int = self.w
        r: cython.int = self.r
        i: cython.int
        j: cython.int

        for i in range(self._step_count):
            VP: cython.ulonglong = self._V(r - 1)
            self.gen_state >>= w
            if VP:
                for j in range(self.nbcoeff):
                    res: cython.ulonglong = self._multiply(
                        VP, cython.cast(cython.ulonglong, self.coeff[j]))
                    res ^= self._V(self.nocoeff[j])
                    self._SetV(self.nocoeff[j], res)

        L: int = self.L
        result: object = (self.gen_state._val >> (self.k - L)) & ((1 << L) - 1)
        return BitVect(L, result)

    def copy(self) -> "GenF2wPolyLCG":
        new: GenF2wPolyLCG = GenF2wPolyLCG(
            self.w, self.r, self.nbcoeff, self.nocoeff, self.coeff,
            self.modM, self.normal_basis, self._step_count, self.L,
        )
        new.gen_state = self.gen_state.copy()
        new.smax = self.smax
        new.table = self.table
        return new
