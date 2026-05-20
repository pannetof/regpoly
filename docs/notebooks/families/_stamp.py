"""Stamp 14 per-family equidistribution-demo notebooks plus the
CombinedGenerator demo from a single skeleton.

Run from the repo root:

    uv run python docs/notebooks/families/_stamp.py

Each generated notebook follows the same explicit recipe:

  1. A paragraph of family-specific intro + links to theory / paper.
  2. Imports.
  3. Construct the generator with literal parameter values from a
     canonical published instance (no `run_for_family` wrapper, no
     catalog walk — the call is the demo).
  4. Wrap the generator (plus any companion tempering) in a single-
     component `Combination`.
  5. Build an `EquidistributionTest` with the family-appropriate
     method (matricial / harase / notprimitive / simd_notprimitive)
     and run it.
  6. Print SE and the per-resolution gap profile, with a brief
     interpretation cell.

Re-running the script overwrites only notebooks that still carry the
`# stamp:auto-generated` marker on their first code cell; manually-
extended notebooks (where the maintainer removed that marker) are
left alone.
"""

from __future__ import annotations

import json
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent

AUTO_MARKER = "# stamp:auto-generated"


def _md(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True),
    }


def _code(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.splitlines(keepends=True),
    }


def _is_auto(path: Path) -> bool:
    """An existing notebook counts as auto-stubbed iff its first code
    cell starts with the AUTO_MARKER. Once the marker is removed, the
    stamper leaves the file alone."""
    try:
        nb = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code":
            src = "".join(cell.get("source", []))
            return src.lstrip().startswith(AUTO_MARKER)
    return False


# ─── Per-family configurations ─────────────────────────────────────────
#
# Each entry contributes a full equidistribution demo for one family.
# Fields:
#   title         — display name in the page H1
#   paper         — short citation
#   library_id    — catalog id (link target only; the demo doesn't use
#                   the catalog)
#   intro         — markdown intro paragraph(s)
#   construct     — Python source for the "build the generator" cell
#                   (must define a variable named `gen`)
#   method_const  — name of the METHOD_* constant to import + use
#   method_blurb  — short markdown explaining why this method
#   tempering     — optional Python source appended to `construct`
#                   that attaches a Transformation; None = skip
#   ce_se         — expected SE value (informational; included as a
#                   one-line comment in the print cell)
#   followup      — optional trailing markdown cell

FAMILY_CONFIGS: dict[str, dict] = {
    # ── MT family ─────────────────────────────────────────────────────
    "MTGen": {
        "title": "MTGen — Mersenne Twister (MT19937) equidistribution",
        "paper": "Matsumoto & Nishimura (1998)",
        "library_id": "mt19937",
        "intro": """
**Mersenne Twister** (Matsumoto & Nishimura, 1998) is the de-facto
standard PRNG for simulation work. The canonical MT19937 parameter
set has state size $k = wr - p = 32 \\cdot 624 - 31 = 19937$ and
period $2^{19937} - 1$.

The default `tempMK2` tempering brings the gap profile down to SE = 0
(maximally equidistributed) at $L = 32$. This notebook builds the
bare generator first, runs the equidistribution test, then adds the
tempering chain and re-runs to show the improvement.

See [generators/MTGen.md](../../generators/MTGen.md) for the family
page and [theory/equidistribution-spec.md](../../theory/equidistribution-spec.md)
for the matricial method.
""",
        "construct": """\
# Canonical MT19937 parameters (Matsumoto & Nishimura, 1998 — Table 2).
gen = Generator.create("MTGen", L=32,
    w=32, r=624, m=397, p=31,
    a=0x9908B0DF)
print(gen.display())
""",
        "tempering": """\
# MT19937 tempering (Matsumoto-Kurita "tempMK2" form).
tempering = Transformation.create("tempMK2",
    w=32, eta=7, mu=15, u=11, l=18,
    b=0x9D2C5680, c=0xEFC60000)
""",
        "method_const": "METHOD_HARASE",
        "method_blurb": "MT19937's characteristic polynomial is primitive (full period), so the **Harase** primal-lattice method is the right choice — it is the fastest method for full-period generators.",
    },

    # ── WELL ──────────────────────────────────────────────────────────
    "WELLGen": {
        "title": "WELLGen — WELL512a equidistribution",
        "paper": "Panneton, L'Ecuyer & Matsumoto (2006)",
        "library_id": "well512a",
        "intro": """
**WELL** (Well Equidistributed Long-period Linear, Panneton-L'Ecuyer-
Matsumoto 2006) addresses MT's slow recovery from zero-heavy states
and its weaker low-resolution equidistribution. The canonical WELL512a
parameter set has $k = wr - p = 32 \\cdot 16 = 512$ and is reported
maximally equidistributed (ME) at $L = 32$.

The family is parameterised by eight matrix recipes $T_0, \\ldots, T_7$
(each picked from a small algebra $\\{M_0, M_1, M_2, M_3, M_4, M_5\\}$)
and three offsets $m_1, m_2, m_3$. See
[generators/WELLGen.md](../../generators/WELLGen.md).
""",
        "construct": """\
# Canonical WELL512a parameters (Panneton-L'Ecuyer-Matsumoto 2006).
gen = Generator.create("WELLGen", L=32,
    w=32, r=16, p=0, m1=13, m2=9, m3=5,
    matrices={
        "T0": {"M": 3, "t": -16},
        "T1": {"M": 3, "t": -15},
        "T2": {"M": 3, "t":  11},
        "T3": {"M": 0},
        "T4": {"M": 3, "t":  -2},
        "T5": {"M": 3, "t": -18},
        "T6": {"M": 2, "t": -28},
        "T7": {"M": 5, "b": 0xDA442D24, "t": -5},
    })
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "WELL is full-period, so the **Harase** method applies.",
    },

    # ── MELG ──────────────────────────────────────────────────────────
    "MELGGen": {
        "title": "MELGGen — MELG-607-64 equidistribution",
        "paper": "Harase & Kimoto (2018)",
        "library_id": "melg607-64",
        "intro": """
**MELG** (Maximally Equidistributed F$_2$-Linear, Harase-Kimoto 2018)
is a 64-bit family designed to be ME on every resolution by
construction. Uses a lagged-feedback recurrence with two shift
parameters $\\sigma_1, \\sigma_2$ and a companion `laggedTempering`
output map that XORs in a delayed state word.

This notebook uses MELG-607-64 ($k = 607$, period $2^{607} - 1$).
See [generators/MELGGen.md](../../generators/MELGGen.md).
""",
        "construct": """\
# Canonical MELG-607-64 parameters (Harase-Kimoto 2018 — Table 1).
gen = Generator.create("MELGGen", L=64,
    w=64, N=10, M=5, r=33,
    sigma1=13, sigma2=35,
    a=0x81F1FD68012348BC)
print(gen.display())
""",
        "tempering": """\
# MELG-style two-reference tempering: y' = (y ^ (y << sigma)) ^
# (state_word_at_lag_L & b).
tempering = Transformation.create("laggedTempering",
    w=64, sigma=30, L=3,
    b=0x66EDC62A6BF8C826)
""",
        "method_const": "METHOD_HARASE",
        "method_blurb": "MELG is full-period; the **Harase** method is the fastest applicable choice.",
    },

    # ── SFMT ──────────────────────────────────────────────────────────
    "SFMTGen": {
        "title": "SFMTGen — SFMT607 equidistribution",
        "paper": "Saito & Matsumoto (2008)",
        "library_id": "sfmt607",
        "intro": """
**SFMT** (SIMD-oriented Fast Mersenne Twister, Saito-Matsumoto 2008)
batches state updates into 128-bit super-words and exploits SIMD
shifts. Its characteristic polynomial is **not primitive** — the
state is $128N$ bits but the period exponent is the prime
$\\mathrm{MEXP}$ — so the equidistribution test must use the
**`simd_notprimitive`** method (Saito-Matsumoto 2008, §3.2),
which projects onto the largest primitive factor and accounts for
the 4 lanes per 128-bit word.

This notebook uses SFMT607 (the smallest variant). See
[generators/SFMTGen.md](../../generators/SFMTGen.md).
""",
        "construct": """\
# SFMT607 parameters (Saito-Matsumoto 2008 — supplementary table).
gen = Generator.create("SFMTGen", L=32,
    mexp=607, pos1=2, sl1=15, sl2=3, sr1=13, sr2=3,
    msk=[
        0xFDFF37FF, 0xEF7F3F7D, 0xFF777B7D, 0x7FF7FB2F,
    ])
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_SIMD_NOTPRIMITIVE",
        "method_blurb": "SFMT's $\\chi_f$ is reducible (state has $128N \\ne $ period exponent bits). The **`simd_notprimitive`** method (Saito-Matsumoto 2008, §3.2) projects onto the largest primitive factor and accounts for the 4 lanes per 128-bit word.",
    },

    # ── dSFMT ─────────────────────────────────────────────────────────
    "DSFMTGen": {
        "title": "DSFMTGen — dSFMT521 equidistribution",
        "paper": "Saito (2009)",
        "library_id": "dsfmt521",
        "intro": """
**dSFMT** (double-precision SFMT, Saito 2009) emits 52-bit mantissas
directly, skipping the 64→52-bit conversion the application would
otherwise do. Structurally similar to SFMT (SIMD-aware, non-primitive
$\\chi_f$), but $L = 52$ to match the IEEE-754 mantissa width.

This notebook uses dSFMT521 (the smallest variant). See
[generators/DSFMTGen.md](../../generators/DSFMTGen.md).
""",
        "construct": """\
# dSFMT521 parameters (Saito 2009 — supplementary table).
gen = Generator.create("DSFMTGen", L=52,
    mexp=521, pos1=3, sl1=25,
    msk1=0x000FBFEFFF77EFFF, msk2=0x000FFEEBFBDFBFDF,
    fix1=0xCFB393D661638469, fix2=0xC166867883AE2ADB,
    pcv1=0xCCAA588000000000, pcv2=0x0000000000000001)
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_SIMD_NOTPRIMITIVE",
        "method_blurb": "Same SIMD-aware path as SFMT.",
    },

    # ── MTGP (no inline demo) ─────────────────────────────────────────
    "MTGPGen": {
        "title": "MTGPGen — Mersenne Twister for Graphic Processors",
        "paper": "Saito & Matsumoto (2013)",
        "library_id": None,
        "intro": """
**MTGP** (Mersenne Twister for Graphic Processors, Saito-Matsumoto
2013) is a SIMD-aware MT variant tuned for GPU thread groups. Its
parameter set includes two 16-entry lookup tables (`tbl` and
`tmp_tbl`) populated from the upstream MTGPDC code generator.

There is no inline parameter set in REGPOLY's catalog yet, so this
demo intentionally stops at the family overview — see the
[generators/MTGPGen.md](../../generators/MTGPGen.md) family page for
the construction recipe and, when a catalog entry lands, this
notebook will be re-stamped to include a full equidistribution run.

The equidistribution method for MTGP is the same as SFMT/dSFMT —
**`simd_notprimitive`** (Saito-Matsumoto 2008, §3.2).
""",
        "construct": None,
        "method_const": "METHOD_SIMD_NOTPRIMITIVE",
    },

    # ── TinyMT32 ──────────────────────────────────────────────────────
    "TinyMT32Gen": {
        "title": "TinyMT32Gen — TinyMT-32 equidistribution",
        "paper": "Saito & Matsumoto (2011)",
        "library_id": "tinymt32-default",
        "intro": """
**TinyMT** (Saito-Matsumoto 2011) is a 127-bit-state MT variant
designed for parameter parallelism — each parameter set
$(\\mathrm{mat}_1, \\mathrm{mat}_2, \\mathrm{tmat})$ defines an
independent generator, useful for spawning many uncorrelated streams.

This notebook uses the upstream default parameter set. Because the
state is small ($k = 127$, output $L = 32$) the matricial method
runs in milliseconds. See
[generators/TinyMT32Gen.md](../../generators/TinyMT32Gen.md).
""",
        "construct": """\
# TinyMT-32 default parameter set.
gen = Generator.create("TinyMT32Gen", L=32,
    mat1=0x8F7011EE, mat2=0xFC78FF1F, tmat=0x3793FDFF)
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_NOTPRIMITIVE",
        "method_blurb": "TinyMT-32's $\\chi_f$ is reducible (state $= 127$ bits but only $127 - 1 = 126$ are reachable from the typical seed); the **`notprimitive`** matricial-on-invariant-subspace method applies.",
    },

    # ── RMT64 (no inline demo) ────────────────────────────────────────
    "RMT64Gen": {
        "title": "RMT64Gen — RMT64 family",
        "paper": "Saito (MTToolBox)",
        "library_id": None,
        "intro": """
**RMT64** is a 64-bit MT-style family with parameterisation drawn
from the upstream MTToolBox `RMT64DC` code generator. Each instance
requires five integers (`mexp`, `pos`, `mata`, `maskb`, `maskc`)
from a precomputed table.

There is no inline parameter set in REGPOLY's catalog yet, so this
demo intentionally stops at the family overview — see the
[generators/RMT64Gen.md](../../generators/RMT64Gen.md) family page
for the construction recipe and, when a catalog entry lands, this
notebook will be re-stamped to include a full equidistribution run.
""",
        "construct": None,
        "method_const": "METHOD_HARASE",
    },

    # ── Tausworthe ────────────────────────────────────────────────────
    "TauswortheGen": {
        "title": "TauswortheGen — taus88 equidistribution",
        "paper": "L'Ecuyer (1996)",
        "library_id": "taus88",
        "intro": """
**Tausworthe** generators (L'Ecuyer 1996) are LFSRs over $\\mathbb{F}_2$
defined by a primitive polynomial and a decimation step $s$. Single
components have modest equidistribution; the textbook combined
construction taus88 XORs three Tausworthe streams to bring the
combined period to $\\approx 2^{88}$ and the equidistribution gap
profile to a small constant.

This notebook builds the *first* component of taus88 ($z^{31} +
z^{13} + 1$, $s = 12$). See
[generators/TauswortheGen.md](../../generators/TauswortheGen.md) for
the full three-component recipe.
""",
        "construct": """\
# taus88 component 1 (L'Ecuyer 1996 — Table III). The full taus88
# combines this with two more Tausworthe streams (k=29, k=28).
gen = Generator.create("TauswortheGen", L=31,
    k=31, nb_terms=3, quicktaus=True,
    poly=[0, 13, 31], s=12)
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "Tausworthe generators are full-period when the polynomial is primitive; **Harase** applies.",
    },

    # ── TGFSR ─────────────────────────────────────────────────────────
    "TGFSRGen": {
        "title": "TGFSRGen — TT800 equidistribution",
        "paper": "Matsumoto & Kurita (1992)",
        "library_id": None,
        "intro": """
**TGFSR** (Twisted Generalised Feedback Shift Register, Matsumoto-Kurita
1992) is the immediate predecessor of MT. The famous TT800 instance
has $k = wr = 800$ and period $2^{800} - 1$.

This notebook uses the canonical TT800 parameter set. See
[generators/TGFSRGen.md](../../generators/TGFSRGen.md).
""",
        "construct": """\
# TT800 parameters (Matsumoto-Kurita 1992).
gen = Generator.create("TGFSRGen", L=32,
    w=32, r=25, m=7, a=0x8EBFD028)
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "TT800 is full-period; **Harase** applies.",
    },

    # ── PolyLCG ───────────────────────────────────────────────────────
    "PolyLCGGen": {
        "title": "PolyLCGGen — polynomial-LCG-over-GF(2) equidistribution",
        "paper": "L'Ecuyer (1999) — combined LFSR family",
        "library_id": None,
        "intro": """
**PolyLCGGen** implements a linear congruential generator over the
polynomial ring $\\mathbb{F}_2[t] / p(t)$ — a unified abstraction
underlying several practical families (combined LFSRs, F2w
generators, etc.). The parameter is the modulus polynomial $p(t)$,
specified as the exponent list with the leading 1 in ascending order.

This notebook uses the primitive trinomial $z^{31} + z^3 + 1$. See
[generators/PolyLCGGen.md](../../generators/PolyLCGGen.md).
""",
        "construct": """\
# PolyLCG modulus = primitive trinomial z^31 + z^3 + 1 over F_2.
gen = Generator.create("PolyLCGGen", L=31,
    k=31, poly=[0, 3, 31])
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "PolyLCG over a primitive polynomial is full-period; **Harase** applies.",
    },

    # ── F2w-LFSR ──────────────────────────────────────────────────────
    "F2wLFSRGen": {
        "title": "F2wLFSRGen — F2w-LFSR equidistribution",
        "paper": "Panneton & L'Ecuyer (2004)",
        "library_id": "f2w-lfsr-t1-r01",
        "intro": """
**F2wLFSR** (Panneton-L'Ecuyer 2004) lifts a Tausworthe-style LFSR
from $\\mathbb{F}_2$ to $\\mathrm{GF}(2^w)$, processing $w$ bits per
recurrence step instead of one. Parameters are the field defining
polynomial `modM`, the recurrence coefficients `coeff` (one element
of $\\mathrm{GF}(2^w)$ per non-zero recurrence term), and the
matching position list `nocoeff`.

This notebook uses the first entry from Panneton-L'Ecuyer 2004
Table 1. See [generators/F2wLFSRGen.md](../../generators/F2wLFSRGen.md).
""",
        "construct": """\
# F2w-LFSR, Panneton-L'Ecuyer 2004 — Table 1 row 1 (r=3, w=32).
gen = Generator.create("F2wLFSRGen", L=32,
    w=32, r=3, nb_terms=2, step=1, normal_basis=False,
    modM=0xCCB06F34,
    coeff=[0x30A2DEE7, 0x53782E5F],
    nocoeff=[1, 0])
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "F2w-LFSR with primitive `modM` and good coefficients is full-period; **Harase** applies.",
    },

    # ── F2w-PolyLCG ───────────────────────────────────────────────────
    "F2wPolyLCGGen": {
        "title": "F2wPolyLCGGen — F2w-PolyLCG equidistribution",
        "paper": "Panneton & L'Ecuyer (2004)",
        "library_id": "f2w-polylcg-t1-r01",
        "intro": """
**F2wPolyLCG** (Panneton-L'Ecuyer 2004) is the PolyLCG counterpart
to `F2wLFSRGen`: a linear congruential recurrence over the ring
$\\mathrm{GF}(2^w)[t] / p(t)$ where the field defining polynomial
`modM` over $\\mathbb{F}_2$ gives the underlying $\\mathrm{GF}(2^w)$.

This notebook uses the first entry from Panneton-L'Ecuyer 2004
Table 1 — identical parameter shape to the F2w-LFSR variant but a
different recurrence. See
[generators/F2wPolyLCGGen.md](../../generators/F2wPolyLCGGen.md).
""",
        "construct": """\
# F2w-PolyLCG, Panneton-L'Ecuyer 2004 — Table 1 row 1 (r=3, w=32).
gen = Generator.create("F2wPolyLCGGen", L=32,
    w=32, r=3, nb_terms=2, step=1, normal_basis=False,
    modM=0xCCB06F34,
    coeff=[0x30A2DEE7, 0x53782E5F],
    nocoeff=[1, 0])
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "Same as F2w-LFSR.",
    },

    # ── MarsaXorshift ─────────────────────────────────────────────────
    "MarsaXorshiftGen": {
        "title": "MarsaXorshiftGen — Marsaglia xorshift triple equidistribution",
        "paper": "Marsaglia (2003)",
        "library_id": "marsaglia-2003-xor",
        "intro": """
**Marsaglia xorshift** (Marsaglia 2003) is a family of three-shift
LFSRs of the form $y \\leftarrow y \\oplus (y \\ll a) \\oplus
(y \\gg b) \\oplus (y \\ll c)$, parameterised by the three shift
amounts. The 32-bit base instance Marsaglia gave in the paper has
shifts $(13, 17, 5)$ and produces $k = 32$ bits of state.

This notebook uses Marsaglia's original 32-bit triple. See
[generators/MarsaXorshiftGen.md](../../generators/MarsaXorshiftGen.md).
""",
        "construct": """\
# Marsaglia (2003) — original 32-bit xor() triple.
# Signed shifts: negative = left shift, positive = right shift.
gen = Generator.create("MarsaXorshiftGen", L=32,
    w=32, r=1, type=1, shifts=[-13, 17, -5])
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "32-bit xorshift over a primitive polynomial is full-period; **Harase** applies.",
    },

    # ── Xoroshiro ─────────────────────────────────────────────────────
    "XoroshiroGen": {
        "title": "XoroshiroGen — xoroshiro128 equidistribution",
        "paper": "Blackman & Vigna (2018, 2022)",
        "library_id": "blackman-vigna-2022-xoroshiro128",
        "intro": """
**xoroshiro** (Blackman-Vigna 2018) is a 128-bit-state lagged-feedback
generator using a rotation in the recurrence. The canonical
xoroshiro128 parameter set $(A, B, C) = (24, 16, 37)$ is the
"good triple" published in the 2022 erratum.

See [generators/XoroshiroGen.md](../../generators/XoroshiroGen.md).
""",
        "construct": """\
# xoroshiro128 — Blackman-Vigna 2022 erratum (good triple).
gen = Generator.create("XoroshiroGen", L=64,
    w=64, r=2, A=24, B=16, C=37)
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "xoroshiro128 is full-period; **Harase** applies.",
    },

    # ── Xoshiro ───────────────────────────────────────────────────────
    "XoshiroGen": {
        "title": "XoshiroGen — xoshiro256 equidistribution",
        "paper": "Blackman & Vigna (2018, 2022)",
        "library_id": "blackman-vigna-2022-xoshiro256",
        "intro": """
**xoshiro** (Blackman-Vigna 2018) is the 256-bit-state cousin of
xoroshiro. xoshiro256 uses parameters $(A, B) = (17, 45)$ from the
2022 erratum.

See [generators/XoshiroGen.md](../../generators/XoshiroGen.md).
""",
        "construct": """\
# xoshiro256 — Blackman-Vigna 2022 erratum.
gen = Generator.create("XoshiroGen", L=64,
    w=64, r=4, A=17, B=45)
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "xoshiro256 is full-period; **Harase** applies.",
    },

    # ── CellularAutomata ──────────────────────────────────────────────
    "CellularAutomataGen": {
        "title": "CellularAutomataGen — Rule-150 CA equidistribution",
        "paper": "Bhuvaneswari & Bhattacharjee (2026)",
        "library_id": "bhuv2026-r1",
        "intro": """
**Cellular Automata** (Bhuvaneswari-Bhattacharjee 2026) — a hybrid
CA over $k$ cells where each cell at position $i$ applies either
Rule 90 or Rule 150 according to a position list, producing a state-
update map with low-cost shift-and-XOR per cell.

This notebook uses the paper's R1 ($k = 32$) configuration. See
[generators/CellularAutomataGen.md](../../generators/CellularAutomataGen.md).
""",
        "construct": """\
# R1 hybrid CA, k=32 (Bhuvaneswari-Bhattacharjee 2026 — Table 1).
# `rule150_positions` lists the cells that apply Rule 150 (the
# remaining cells apply Rule 90).
gen = Generator.create("CellularAutomataGen", L=32,
    k=32, s=1,
    rule150_positions=[1, 5, 6, 12, 15, 16, 18, 19, 20,
                       22, 23, 24, 25, 27, 29, 31])
print(gen.display())
""",
        "tempering": None,
        "method_const": "METHOD_HARASE",
        "method_blurb": "Hybrid Rule-150/Rule-90 CAs have a maximal-period parameter regime; **Harase** applies there.",
    },
}


# ─── Common cell builders ──────────────────────────────────────────────

IMPORTS_CELL = """\
{auto}
from regpoly.core.generator import Generator
from regpoly.core.combination import Combination
from regpoly.core.transformation import Transformation
from regpoly.analyses.equidistribution_test import (
    EquidistributionTest,
    METHOD_MATRICIAL, METHOD_HARASE,
    METHOD_NOTPRIMITIVE, METHOD_SIMD_NOTPRIMITIVE,
)

INT_MAX = 2**31 - 1
"""


def _construct_cell(cfg: dict) -> dict:
    """Generator-construction code cell."""
    src = cfg["construct"]
    if cfg.get("tempering"):
        src += "\n" + cfg["tempering"]
    return _code(src)


COMB_CELL_NO_TEMPERING = """\
# Wrap the generator in a single-component Combination. The
# Combination is the live object the search loop iterates and the
# equidistribution test consumes.
comb = Combination(J=1, Lmax=gen.L)
comb.components[0].add_gen(gen)
comb.reset()
print(f"k_g = {comb.k_g}, L = {comb.L}")
"""

COMB_CELL_WITH_TEMPERING = """\
# Wrap the generator + tempering chain in a single-component Combination.
comb = Combination(J=1, Lmax=gen.L)
comb.components[0].add_gen(gen)
comb.components[0].add_trans(tempering)
comb.reset()
print(f"k_g = {comb.k_g}, L = {comb.L}")
"""


def _test_cell(cfg: dict) -> dict:
    method = cfg["method_const"]
    src = f"""\
# Build the equidistribution test and run it. We cap `delta` at
# INT_MAX so nothing is rejected — we just want the score.
test = EquidistributionTest(
    L=gen.L,
    delta=[INT_MAX] * (gen.L + 1),
    mse=INT_MAX,
    method={method},
)
result = test.run(comb)

print(f"SE (Σ gaps)   = {{result.se}}")
print(f"verified      = {{result.verified}}")
print(f"first 10 gaps = {{[result.ecart[i] for i in range(1, min(11, gen.L + 1))]}}")
"""
    return _code(src)


def _build_family_notebook(family: str, cfg: dict) -> dict:
    """Build a 14-family equidistribution-demo notebook."""
    title_md = f"# {cfg['title']}\n\n{cfg['intro'].strip()}\n"

    if cfg["construct"] is None:
        # Families with no inline parameter set yet — stop after the intro.
        return {
            "cells": [
                _md(title_md),
                _code(IMPORTS_CELL.format(auto=AUTO_MARKER)),
                _md(
                    "_No inline demo: this family needs upstream lookup tables. "
                    "When a catalog entry lands, re-running "
                    "`uv run python docs/notebooks/families/_stamp.py` will fill "
                    "in the rest of the cells._\n"
                ),
            ],
            "metadata": _notebook_metadata(),
            "nbformat": 4,
            "nbformat_minor": 5,
        }

    method_blurb = (
        "## Equidistribution test\n\n"
        + cfg.get("method_blurb", "")
        + "\n"
    )

    cells = [
        _md(title_md),
        _md("## Imports\n"),
        _code(IMPORTS_CELL.format(auto=AUTO_MARKER)),
        _md(f"## Construct the generator — _{cfg['paper']}_\n"),
        _construct_cell(cfg),
        _md("## Wrap in a `Combination`\n"),
        _code(
            COMB_CELL_WITH_TEMPERING
            if cfg.get("tempering")
            else COMB_CELL_NO_TEMPERING
        ),
        _md(method_blurb),
        _test_cell(cfg),
    ]

    if cfg.get("library_id"):
        cells.append(
            _md(
                "## Catalog entry\n\n"
                f"The published version of this parameter set lives in the "
                f"REGPOLY catalog under `library_id = "
                f'"{cfg["library_id"]}"`. To load it programmatically without '
                "hard-coding parameters:\n\n"
                "```python\n"
                "from regpoly.library import Catalog\n"
                "cat = Catalog('docs/library')\n"
                "cat.load()\n"
                f"_, entry = cat.generator({cfg['library_id']!r})\n"
                "# entry.components[0] carries the same params as constructed above\n"
                "```\n"
            )
        )

    return {
        "cells": cells,
        "metadata": _notebook_metadata(),
        "nbformat": 4,
        "nbformat_minor": 5,
    }


# ─── CombinedGenerator (separate demo shape) ───────────────────────────

COMBINED_TITLE_MD = """\
# CombinedGenerator — XORing two MT19937s

`CombinedGenerator` XORs $J$ component generators bitwise to widen
the state space without paying the per-step cost of a single huge
generator. The combined characteristic polynomial is the product of
the components' polynomials (when pairwise coprime), so a
2-component combination of MT19937 with two distinct twist
parameters has state size $k_g = 2 \\cdot 19937 = 39874$.

This notebook stamps two MT19937 instances with different twist
parameters $a$, combines them, and runs the equidistribution test on
the resulting 2-component combination.

See [theory/f2-linear-generators.md#combined-generators](../../theory/f2-linear-generators.md#combined-generators).
"""

COMBINED_BODY = """\
{auto}
from regpoly.core.generator import Generator
from regpoly.core.combination import Combination
from regpoly.analyses.equidistribution_test import (
    EquidistributionTest, METHOD_HARASE,
)
INT_MAX = 2**31 - 1

# Two MT19937 components with different `a` twist constants.
gen_a = Generator.create("MTGen", L=32,
    w=32, r=624, m=397, p=31, a=0x9908B0DF)  # canonical MT19937
gen_b = Generator.create("MTGen", L=32,
    w=32, r=624, m=397, p=31, a=0xDFFFFFFF)  # alternative twist

print(f"component A: k={{gen_a.k}}, L={{gen_a.L}}")
print(f"component B: k={{gen_b.k}}, L={{gen_b.L}}")

# 2-component Combination (J = 2).
comb = Combination(J=2, Lmax=32)
comb.components[0].add_gen(gen_a)
comb.components[1].add_gen(gen_b)
comb.reset()
print(f"combined: k_g = {{comb.k_g}}, L = {{comb.L}}")

# Equidistribution test on the combination — MT19937 has primitive
# χ_f, so the product is primitive too and Harase applies.
test = EquidistributionTest(
    L=comb.L,
    delta=[INT_MAX] * (comb.L + 1),
    mse=INT_MAX,
    method=METHOD_HARASE,
)
result = test.run(comb)
print(f"SE (Σ gaps) = {{result.se}}")
print(f"verified    = {{result.verified}}")
""".format(auto=AUTO_MARKER)


def _build_combined_notebook() -> dict:
    return {
        "cells": [_md(COMBINED_TITLE_MD), _code(COMBINED_BODY)],
        "metadata": _notebook_metadata(),
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def _notebook_metadata() -> dict:
    return {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python"},
    }


# ─── Driver ────────────────────────────────────────────────────────────


def main() -> int:
    written = 0
    skipped = 0

    for family, cfg in FAMILY_CONFIGS.items():
        out = OUTPUT_DIR / f"{family}.ipynb"
        if out.exists() and not _is_auto(out):
            skipped += 1
            continue
        nb = _build_family_notebook(family, cfg)
        out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
        written += 1

    out = OUTPUT_DIR / "CombinedGenerator.ipynb"
    if out.exists() and not _is_auto(out):
        skipped += 1
    else:
        out.write_text(json.dumps(_build_combined_notebook(), indent=1) + "\n", encoding="utf-8")
        written += 1

    print(f"stamped {written} notebook(s); skipped {skipped} (manually edited)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
