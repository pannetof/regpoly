"""Smoke tests to verify the package loads and basic operations work."""

from regpoly import BitVect, BitMatrix, Generator, Transformation
from regpoly.core.generator import resolve_family


# ── BitVect ──────────────────────────────────────────────────────────────

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


# ── BitMatrix ────────────────────────────────────────────────────────────

def test_bitmatrix_create():
    m = BitMatrix(4, 8)
    assert m.nblignes == 4
    assert m.l == 8
    m.set_bit(0, 0, 1)
    assert m.get_bit(0, 0) == 1
    assert m.get_bit(0, 1) == 0


def test_bitmatrix_display():
    m = BitMatrix(2, 4)
    m.set_bit(0, 0, 1)
    m.set_bit(1, 3, 1)
    text = m.display("raw")
    assert "1" in text


# ── Generator.create ────────────────────────────────────────────────────

def test_create_polylcg():
    gen = Generator.create("PolyLCG", 3, k=3, poly=[1])
    assert gen.k == 3
    assert gen.L == 3
    assert gen.name() != ""


def test_polylcg_iteration():
    gen = Generator.create("PolyLCG", 3, k=3, poly=[1])
    init_bv = BitVect.zeros(3)
    init_bv.put_bit(0, 1)
    gen.initialize_state(init_bv)
    out = next(gen)
    assert out.n == 3


def test_create_tausworthe():
    # Canonical poly form is [0, q_1, ..., q_{t-2}, k]; nb_terms must be odd >= 3.
    gen = Generator.create("Tausworthe", 7, poly=[0, 3, 7], s=3, quicktaus=True)
    assert gen.k == 7


def test_create_tausworthe_auto_s():
    gen = Generator.create("Tausworthe", 7, poly=[0, 3, 7])
    assert gen.k == 7


def test_create_tgfsr():
    gen = Generator.create("TGFSR", 8, w=8, r=3, m=1, a=0b10110011)
    assert gen.k == 24  # w * r


def test_create_with_legacy_name():
    gen = Generator.create("polylcg", 3, k=3, poly=[1])
    assert gen.k == 3


def test_create_with_random_params():
    gen = Generator.create("TGFSR", 32, w=32, r=3)
    assert gen.k == 96  # w * r = 32 * 3


def test_generateur_display_returns_string():
    gen = Generator.create("PolyLCG", 3, k=3, poly=[1])
    result = gen.display()
    assert isinstance(result, str)


def test_char_poly():
    gen = Generator.create("PolyLCG", 3, k=3, poly=[1])
    init_bv = BitVect.zeros(3)
    init_bv.put_bit(0, 1)
    gen.initialize_state(init_bv)
    poly = gen.char_poly()
    assert poly.n == 3


def test_transition_matrix():
    gen = Generator.create("PolyLCG", 3, k=3, poly=[1])
    mat = gen.transition_matrix()
    assert isinstance(mat, BitMatrix)
    assert mat.nblignes == 3
    assert mat.l == 3


# ── Generator.parameters ────────────────────────────────────────────────

def test_parameters_tgfsr():
    specs = Generator.parameters("TGFSR")
    names = [s["name"] for s in specs]
    assert "w" in names
    assert "r" in names
    assert "m" in names
    assert "a" in names
    w_spec = next(s for s in specs if s["name"] == "w")
    assert w_spec["structural"] is True
    m_spec = next(s for s in specs if s["name"] == "m")
    assert m_spec["structural"] is False
    assert m_spec["rand_type"] == "range"


def test_parameters_legacy_alias():
    specs = Generator.parameters("tgfsr")
    assert len(specs) == len(Generator.parameters("TGFSR"))


# ── resolve_family ───────────────────────────────────────────────────────

def test_resolve_family_aliases():
    assert resolve_family("polylcg") == "PolyLCGGen"
    assert resolve_family("taus") == "TauswortheGen"
    assert resolve_family("taus2") == "TauswortheGen"
    assert resolve_family("tgfsr") == "TGFSRGen"
    assert resolve_family("MT") == "MTGen"
    assert resolve_family("matsumoto") == "MatsumotoGen"
    assert resolve_family("marsaxorshift") == "MarsaXorshiftGen"
    assert resolve_family("AC1D") == "AC1DGen"
    assert resolve_family("carry") == "WELLGen"


def test_resolve_family_passthrough():
    # Old canonical names resolve forward to the new -Gen names.
    assert resolve_family("PolyLCG") == "PolyLCGGen"
    assert resolve_family("MersenneTwister") == "MTGen"
    # New canonical names pass through unchanged.
    assert resolve_family("PolyLCGGen") == "PolyLCGGen"
    assert resolve_family("MTGen") == "MTGen"


def test_resolve_family_genf2w():
    assert resolve_family("genf2w", {"type": "lfsr"}) == "F2wLFSRGen"
    assert resolve_family("genf2w", {"type": "polylcg"}) == "F2wPolyLCGGen"
    assert resolve_family("genf2w") == "F2wPolyLCGGen"


# ── Transformation ───────────────────────────────────────────────────────

def test_create_permutation():
    t = Transformation.create("permut", w=8, p=3, q=1)
    assert t.w == 8
    assert isinstance(t.display(), str)


def test_transformation_copy():
    t = Transformation.create("permut", w=8, p=3, q=1)
    t2 = t.copy()
    assert t2.w == t.w
    assert t2 is not t


# ── Analyses imports ─────────────────────────────────────────────────────

def test_analyses_importable():
    from regpoly.analyses import (
        EquidistributionTest,
        CollisionFreeTest,
        TupletsTest,
    )
    assert EquidistributionTest is not None
    assert CollisionFreeTest is not None
    assert TupletsTest is not None
