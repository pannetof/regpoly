"""Smoke tests to verify the package loads and basic operations work."""

from regpoly import BitVect, Matrix
from regpoly.generators import PolyLCG, Tausworthe, TGFSR
from regpoly.transformations import TemperMK, Permutation
from regpoly.tests import EquidistributionTest, CollisionFreeTest, TupletsTest


def test_bitvect_basics():
    bv = BitVect.zeros(32)
    bv.put_bit(0, 1)
    assert bv.get_bit(0) == 1
    assert bv.get_bit(1) == 0


def test_bitvect_xor():
    a = BitVect(8, 0b10110010)
    b = BitVect(8, 0b01101001)
    c = a ^ b
    assert c._val == 0b11011011


def test_bitvect_shift():
    bv = BitVect(8, 0b10000000)
    shifted = bv >> 3
    assert shifted._val == 0b00010000


def test_matrix_create():
    m = Matrix(4, 8, 1)
    assert m.nblignes == 4
    assert m.l == 8
    assert m.t == 1


def test_generators_importable():
    assert PolyLCG.name() == "Polynomial LCG"
    assert Tausworthe.name() == "Tausworthe Generator"
    assert TGFSR.name() == "TGFSR"


def test_polylcg_iteration():
    poly = BitVect(3, 0b011)  # x^3 + x + 1 → poly bits = 011
    gen = PolyLCG(3, poly, 3)
    init = BitVect.zeros(3)
    init.put_bit(0, 1)
    gen.initialize_state(init)
    out = next(gen)
    assert out.n == 3
