// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "marsaxorshift.h"

#include "gen_enumerator.h"

#include <NTL/ZZ.h>

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <stdexcept>

using namespace regpoly::core;


namespace regpoly::core {

namespace {

// Validate (type, w, r, payload) consistency.  Throws
// std::invalid_argument with a descriptive message when the payload
// doesn't match the runtime-type tag.  Returns `type` so it can sit
// in the initializer list (forcing validation to happen before any
// `const` field is bound).
int validate_ctor_args_(
    int type, int w, int r,
    const MarsaXorshiftGen::Type2xParams& t2x,
    const std::vector<MarsaXorshiftGen::Tap>& taps,
    const MarsaXorshiftGen::Type4Params& t4,
    const std::vector<MarsaXorshiftGen::MiEntry>& mi)
{
    if (type != 1 && type != 2 && type != 3 && type != 4 && type != 100)
        throw std::invalid_argument(
            "MarsaXorshiftGen: type must be one of {1, 2, 3, 4, 100}");
    if (w < 1)
        throw std::invalid_argument("MarsaXorshiftGen: w must be >= 1");
    // No upper bound: w > 64 routes through the BitVect-based wide
    // path in next() (see next_wide_()), which has no fixed register
    // width.  For w <= 64 the fast uint64_t path runs instead.
    if (r < 1)
        throw std::invalid_argument("MarsaXorshiftGen: r must be >= 1");
    if (type == 1 && r != 1)
        throw std::invalid_argument(
            "MarsaXorshiftGen: type=1 requires r=1");
    if (type == 2 && (t2x.p.size() != 3 || t2x.q.size() != 3))
        throw std::invalid_argument(
            "MarsaXorshiftGen: type=2 requires p and q vectors of length 3");
    if (type == 3 && taps.empty())
        throw std::invalid_argument(
            "MarsaXorshiftGen: type=3 requires a non-empty taps list");
    if (type == 4 && (t4.p.size() != 2 || t4.q.size() != 2))
        throw std::invalid_argument(
            "MarsaXorshiftGen: type=4 requires p and q vectors of length 2");
    if (type == 100) {
        if (mi.empty())
            throw std::invalid_argument(
                "MarsaXorshiftGen: type=100 requires a non-empty mi list");
        for (const auto& e : mi) {
            if (e.shifts.empty())
                throw std::invalid_argument(
                    "MarsaXorshiftGen: type=100 requires every MiEntry to "
                    "carry at least one shift");
        }
    }
    return type;
}

uint64_t compute_wmask_(int w) {
    return (w >= 64) ? ~uint64_t{0} : ((uint64_t{1} << w) - 1);
}

// Render a vector of ints into the stream as "[a, b, c, ...]".  Used
// by display_str() to avoid the manual comma-join boilerplate that
// previously appeared in every type branch.
template <class V>
void render_int_vec_(std::ostringstream& oss, const V& vec) {
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << vec[i];
    }
    oss << "]";
}

}  // anonymous namespace

MarsaXorshiftGen::MarsaXorshiftGen(
    int type, int w, int r, int m,
    const Type1Params& t1,
    const Type2xParams& t2x,
    const std::vector<Tap>& taps,
    const Type4Params& t4,
    const std::vector<MiEntry>& mi,
    int L)
    : Generator(w * r, L),
      // validate_ctor_args_ runs first (in the initializer list) so the
      // const fields below can never be bound to an inconsistent payload.
      type_(validate_ctor_args_(type, w, r, t2x, taps, t4, mi)),
      w_(w), r_(r), m_(m),
      t1_(t1), t2x_(t2x), taps_(taps), t4_(t4), mi_(mi),
      wmask_(compute_wmask_(w))
{
    state_ = BitVect(k_);
}

std::string MarsaXorshiftGen::name() const { return "Marsaglia Xor-shift"; }

std::string MarsaXorshiftGen::display_str() const {
    std::ostringstream oss;
    oss << "type = " << type_ << ", w = " << w_ << ", r = " << r_;

    if (type_ == 1) {
        oss << ", a = " << t1_.a << ", b = " << t1_.b << ", c = " << t1_.c;
    } else if (type_ == 2) {
        oss << ", m = " << m_ << ", p = ";
        render_int_vec_(oss, t2x_.p);
        oss << ", q = ";
        render_int_vec_(oss, t2x_.q);
    } else if (type_ == 3) {
        oss << ", taps = [";
        for (size_t i = 0; i < taps_.size(); i++) {
            if (i > 0) oss << ", ";
            oss << "(" << taps_[i].position << ", " << taps_[i].shift << ")";
        }
        oss << "]";
    } else if (type_ == 4) {
        oss << ", m = " << m_ << ", p = ";
        render_int_vec_(oss, t4_.p);
        oss << ", q = ";
        render_int_vec_(oss, t4_.q);
    } else if (type_ == 100) {
        oss << ", mi = [";
        for (size_t i = 0; i < mi_.size(); i++) {
            if (i > 0) oss << ", ";
            oss << "{pos=" << mi_[i].position << ", shifts=";
            render_int_vec_(oss, mi_[i].shifts);
            oss << "}";
        }
        oss << "]";
    }
    return oss.str();
}

void MarsaXorshiftGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

uint64_t MarsaXorshiftGen::ShiftR(uint64_t v, int s) {
    if (s > 0) return v >> s;
    return v << (-s);
}

// V / SetV index w-bit blocks of the state via the fast uint64_t
// kernel.  Only safe for w <= 64; the wide path below uses
// BitVect::get_slice / set_slice instead and has no width cap.
uint64_t MarsaXorshiftGen::V(int idx) const {
    return state_.get_word(idx, w_);
}

void MarsaXorshiftGen::SetV(int idx, uint64_t val) {
    state_.set_word(idx, w_, val & wmask_);
}

namespace {

// Wide-path xorshift: in-place v = v XOR ShiftR(v, s) on a w-bit
// BitVect.  Sign convention matches the uint64_t ShiftR helper —
// positive s shifts right, negative s shifts left.  Zero is a no-op
// (matches the Type-2 "skip this factor" semantic).
void wide_xorshift_inplace_(BitVect& v, int s) {
    if (s == 0) return;
    BitVect tmp = v.copy();
    if (s > 0) tmp.rshift(s);
    else       tmp.lshift(-s);
    v.xor_with(tmp);
}

}  // anonymous namespace

void MarsaXorshiftGen::next() {
    // Dispatch: native uint64_t arithmetic for w <= 64, BitVect for
    // anything wider.  Both paths preserve the same recurrence
    // semantics (w-bit masking after every shift, etc.).
    if (w_ > 64) { next_wide_(); return; }

    if (type_ == 1) {
        // Single word xorshift
        uint64_t v = V(0);
        v ^= ShiftR(v, t1_.a); v &= wmask_;
        v ^= ShiftR(v, t1_.b); v &= wmask_;
        v ^= ShiftR(v, t1_.c); v &= wmask_;
        SetV(0, v);

    } else if (type_ == 2) {
        // Two-component, r words: v_n = G v_{n-m} XOR H v_{n-r}, where
        // p[3] are the shifts that build H and q[3] are the shifts that
        // build G (zero entries mean "skip this factor").  Intermediate
        // values are masked to w bits after every shift so the
        // recurrence matches the paper's `unsigned int` arithmetic
        // (cf. Panneton & L'Ecuyer 2005, Tables III–IV): when a tap
        // composes a left shift with a right shift, the upper bits
        // produced by the left shift must be discarded before the
        // right shift can pull them down.
        uint64_t y = V(m_ - 1);
        uint64_t x = V(r_ - 1);
        state_.rshift(w_);
        for (int j = 0; j < 3; j++) {
            if (t2x_.p[j] != 0) { x ^= ShiftR(x, t2x_.p[j]); x &= wmask_; }
        }
        for (int j = 0; j < 3; j++) {
            if (t2x_.q[j] != 0) { y ^= ShiftR(y, t2x_.q[j]); y &= wmask_; }
        }
        SetV(0, x ^ y);

    } else if (type_ == 3) {
        // Multi-tap.  At w=32 this preserves Marsaglia's deliberate
        // cross-word read (his C code does `*(ulong*)&vect[i]` which
        // on little-endian gives vect[i] | (vect[i+1] << 32) — i.e. a
        // 64-bit operand built from two adjacent 32-bit words).  For
        // any other w the cross-word trick is meaningless (V already
        // returns a w-bit block), so the tap collapses to a regular
        // one-shift xorshift.  Only one shift per tap, so no
        // inter-shift masking is needed; SetV applies wmask_.
        uint64_t t = 0;
        const bool legacy_cross_word = (w_ == 32);
        for (const auto& tap : taps_) {
            int idx = tap.position - 1;
            uint64_t v_lo = V(idx);
            uint64_t shifted;
            if (legacy_cross_word) {
                uint64_t v64 = v_lo;
                if (idx + 1 < r_)
                    v64 |= V(idx + 1) << 32;
                shifted = ShiftR(v64, tap.shift);
            } else {
                shifted = ShiftR(v_lo, tap.shift);
            }
            t ^= v_lo ^ shifted;
        }
        state_.rshift(w_);
        SetV(0, t);

    } else if (type_ == 4) {
        // Two-component, 2 shifts each.  Mask between shifts (see Type 2x note).
        uint64_t y = V(m_ - 1);
        uint64_t x = V(r_ - 1);
        state_.rshift(w_);
        x ^= ShiftR(x, t4_.p[0]); x &= wmask_;
        x ^= ShiftR(x, t4_.p[1]); x &= wmask_;
        y ^= ShiftR(y, t4_.q[0]); y &= wmask_;
        y ^= ShiftR(y, t4_.q[1]); y &= wmask_;
        SetV(0, x ^ y);

    } else if (type_ == 100) {
        // General multi-tap with arbitrary shift counts per tap.  Mask
        // between shifts (see Type 2x note).
        uint64_t t = 0;
        for (const auto& entry : mi_) {
            uint64_t temp = V(entry.position - 1);
            for (int shift : entry.shifts) {
                temp ^= ShiftR(temp, shift);
                temp &= wmask_;
            }
            t ^= temp;
        }
        state_.rshift(w_);
        SetV(0, t);
    }
}

void MarsaXorshiftGen::next_wide_() {
    // Mirrors next()'s per-type branches, but every word is a w-bit
    // BitVect (no uint64_t register cap).  Block layout in state_
    // matches the fast path: bit range [idx*w, (idx+1)*w) is block idx,
    // with block 0 being the most-recent word and block r-1 the oldest.

    if (type_ == 1) {
        BitVect v = state_.get_slice(0, w_);
        wide_xorshift_inplace_(v, t1_.a);
        wide_xorshift_inplace_(v, t1_.b);
        wide_xorshift_inplace_(v, t1_.c);
        state_.set_slice(0, v);

    } else if (type_ == 2) {
        BitVect y = state_.get_slice((m_ - 1) * w_, w_);
        BitVect x = state_.get_slice((r_ - 1) * w_, w_);
        state_.rshift(w_);
        for (int j = 0; j < 3; j++) wide_xorshift_inplace_(x, t2x_.p[j]);
        for (int j = 0; j < 3; j++) wide_xorshift_inplace_(y, t2x_.q[j]);
        x.xor_with(y);
        state_.set_slice(0, x);

    } else if (type_ == 3) {
        // Type 3's deliberate cross-word read is a w=32-only artifact;
        // the wide path can never hit it (w > 64 here).  Each tap is a
        // straightforward one-shift xorshift on the corresponding
        // w-bit block.
        BitVect t(w_);
        for (const auto& tap : taps_) {
            BitVect v = state_.get_slice((tap.position - 1) * w_, w_);
            BitVect shifted = v.copy();
            if (tap.shift > 0)      shifted.rshift(tap.shift);
            else if (tap.shift < 0) shifted.lshift(-tap.shift);
            v.xor_with(shifted);
            t.xor_with(v);
        }
        state_.rshift(w_);
        state_.set_slice(0, t);

    } else if (type_ == 4) {
        BitVect y = state_.get_slice((m_ - 1) * w_, w_);
        BitVect x = state_.get_slice((r_ - 1) * w_, w_);
        state_.rshift(w_);
        wide_xorshift_inplace_(x, t4_.p[0]);
        wide_xorshift_inplace_(x, t4_.p[1]);
        wide_xorshift_inplace_(y, t4_.q[0]);
        wide_xorshift_inplace_(y, t4_.q[1]);
        x.xor_with(y);
        state_.set_slice(0, x);

    } else if (type_ == 100) {
        BitVect t(w_);
        for (const auto& entry : mi_) {
            BitVect temp = state_.get_slice((entry.position - 1) * w_, w_);
            for (int s : entry.shifts) wide_xorshift_inplace_(temp, s);
            t.xor_with(temp);
        }
        state_.rshift(w_);
        state_.set_slice(0, t);
    }
}

std::unique_ptr<Generator> MarsaXorshiftGen::copy() const {
    auto g = std::make_unique<MarsaXorshiftGen>(
        type_, w_, r_, m_, t1_, t2x_, taps_, t4_, mi_, L_);
    g->state_ = state_.copy();
    return g;
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generator> MarsaXorshiftGen::from_params(const Params& params, int L) {
    int type = (int)params.get_int("type");
    int w = (int)params.get_int("w", 32);
    int r = (int)params.get_int("r", 1);
    int m = (int)params.get_int("m", 0);

    Type1Params t1{};
    Type2xParams t2x;
    std::vector<Tap> taps;
    Type4Params t4;
    std::vector<MiEntry> mi;

    if (type == 1) {
        auto abc = params.get_int_vec("shifts");
        t1.a = abc.size() > 0 ? abc[0] : 0;
        t1.b = abc.size() > 1 ? abc[1] : 0;
        t1.c = abc.size() > 2 ? abc[2] : 0;
    } else if (type == 2) {
        t2x.p = params.get_int_vec("p");
        t2x.q = params.get_int_vec("q");
        t2x.p.resize(3, 0);
        t2x.q.resize(3, 0);
    } else if (type == 3) {
        auto tap_pos = params.get_int_vec("tap_positions");
        auto tap_shifts = params.get_int_vec("tap_shifts");
        for (size_t i = 0; i < tap_pos.size(); i++) {
            taps.push_back({tap_pos[i],
                (i < tap_shifts.size()) ? tap_shifts[i] : 0});
        }
    } else if (type == 4) {
        t4.p = params.get_int_vec("p");
        t4.q = params.get_int_vec("q");
        t4.p.resize(2, 0);
        t4.q.resize(2, 0);
    } else if (type == 100) {
        auto mi_pos = params.get_int_vec("mi_positions");
        auto mi_shifts = params.get_int_vec("mi_shifts");
        auto mi_counts = params.get_int_vec("mi_counts");
        int offset = 0;
        for (size_t i = 0; i < mi_pos.size(); i++) {
            MiEntry entry;
            entry.position = mi_pos[i];
            int count = (i < mi_counts.size()) ? mi_counts[i] : 1;
            for (int j = 0; j < count && offset < (int)mi_shifts.size(); j++) {
                entry.shifts.push_back(mi_shifts[offset++]);
            }
            mi.push_back(entry);
        }
    }

    return std::make_unique<MarsaXorshiftGen>(
        type, w, r, m, t1, t2x, taps, t4, mi, L);
}

std::vector<ParamSpec> MarsaXorshiftGen::param_specs() {
    // The vector params below are mutually exclusive across `type`
    // values; only one set is consumed by `from_params` per instance.
    // They are declared with `has_default=true` so the Python fill-
    // params layer leaves the unused ones empty rather than rejecting
    // the config for "no random generation available".
    return {
        {"type",           "int",     true,  false, 0,  "",     "", false},
        {"w",              "int",     true,  true,  32, "",     "", false},
        {"r",              "int",     true,  true,  1,  "",     "", false},
        {"m",              "int",     true,  true,  0,  "",     "", false},
        {"shifts",         "int_vec", false, true,  0,  "none", "", false},
        {"p",              "int_vec", false, true,  0,  "none", "", false},
        {"q",              "int_vec", false, true,  0,  "none", "", false},
        {"tap_positions",  "int_vec", false, true,  0,  "none", "", false},
        {"tap_shifts",     "int_vec", false, true,  0,  "none", "", false},
        {"mi_positions",   "int_vec", false, true,  0,  "none", "", false},
        {"mi_shifts",      "int_vec", false, true,  0,  "none", "", false},
        {"mi_counts",      "int_vec", false, true,  0,  "none", "", false},
        // Enumerator-time pins for type=100 (ignored by from_params; the
        // primitive-search driver strips them before instantiation):
        //   nb_taps:         number of distinct tap positions (default 3)
        //   shifts_per_tap:  per-tap shift count vector, length == nb_taps
        //                    (default a vector of 1s).
        {"nb_taps",        "int",     false, true,  3,  "none", "", false},
        {"shifts_per_tap", "int_vec", false, true,  0,  "none", "", false},
    };
}

// ── Exhaustive-search enumerator ───────────────────────────────────────
//
// Five per-type subclasses share a small base that holds the
// (axes, sizes, total) plumbing — same role-shape as
// TausworthePolyEnumerator (one enumerator per recurrence variant),
// but here type 1..4 and 100 each have distinct axis sets, so they
// each own a subclass with only the fields it needs.
//
// All counts are NTL::ZZ so the largest configurations (Type 2 with
// w=32 ≈ 6.3·10^10, Type 100 with non-trivial shifts_per_tap ≈
// 2·10^11) fit without truncation.

namespace {

// ── Helpers ────────────────────────────────────────────────────────────

// Signed-shift index decoders.  These are the inverse of the
// mixed-radix digit layout each enumerator uses; keeping them as free
// inline functions avoids member-overload noise.

// Alphabet {-(w-1), ..., 0, ..., w-1} of size 2w-1.  Index 0 → -(w-1).
inline int decode_signed_with_zero(int w, long idx) {
    return (int)idx - (w - 1);
}

// Alphabet {-(w-1), ..., -1, 1, ..., w-1} of size 2(w-1).  Index 0 →
// -(w-1); the "zero" slot is removed.
inline int decode_signed_nonzero(int w, long idx) {
    return (idx < w - 1) ? (int)idx - (w - 1) : (int)idx - (w - 2);
}

// Common GenEnumerator scaffolding shared by every MarsaXorshift
// subclass.  Subclasses populate axes_/sizes_/total_ from their ctor
// and only need to implement at().
class MarsaXorshiftEnumeratorBase : public GenEnumerator {
public:
    std::string size_dec() const final { return zz_to_dec(total_); }
    std::vector<Axis> axes()  const final { return axes_; }

protected:
    void push_axis_(const std::string& name, const NTL::ZZ& size,
                    const std::string& describe) {
        sizes_.push_back(size);
        Axis a;
        a.name = name;
        a.size_dec = zz_to_dec(size);
        a.describe = describe;
        axes_.push_back(std::move(a));
    }

    NTL::ZZ check_idx_(const std::string& idx_dec) const {
        NTL::ZZ idx = parse_zz(idx_dec);
        if (idx < 0 || idx >= total_)
            throw std::out_of_range(
                "MarsaXorshiftEnumerator: idx out of range");
        return idx;
    }

    std::vector<NTL::ZZ> sizes_;
    std::vector<Axis>    axes_;
    NTL::ZZ              total_;
};

// ── Type 1 ─────────────────────────────────────────────────────────────
//
// Axes (outermost first): pattern (4) × a (w-1) × b (w-1) × c (w-1).
// Pattern ∈ {X1, X2, X3, X5} per Proposition 4.5 of Panneton &
// L'Ecuyer (2005) — the four equivalence-class representatives whose
// equidistribution profiles are distinct.  Magnitudes a, b, c are
// in [1, w-1].

class MarsaXorshiftType1Enumerator final : public MarsaXorshiftEnumeratorBase {
public:
    explicit MarsaXorshiftType1Enumerator(int w) : w_(w) {
        const long mag = w_ - 1;
        push_axis_("pattern", NTL::ZZ(4),
                   "Type-I equivalence-class rep X_i ∈ {X1, X2, X3, X5}");
        push_axis_("a", NTL::ZZ(mag), "shift magnitude a ∈ [1, w-1]");
        push_axis_("b", NTL::ZZ(mag), "shift magnitude b ∈ [1, w-1]");
        push_axis_("c", NTL::ZZ(mag), "shift magnitude c ∈ [1, w-1]");
        total_ = NTL::ZZ(4) * NTL::ZZ(mag) * NTL::ZZ(mag) * NTL::ZZ(mag);
    }

    Params at(const std::string& idx_dec) const override {
        const NTL::ZZ idx = check_idx_(idx_dec);
        auto digits = mixed_radix_decode(sizes_, idx);
        const long pat = NTL::to_long(digits[0]);
        const int a = (int)NTL::to_long(digits[1]) + 1;
        const int b = (int)NTL::to_long(digits[2]) + 1;
        const int c = (int)NTL::to_long(digits[3]) + 1;
        // Mapping from paper's X_i to (t1.a, t1.b, t1.c) — the
        // application-order encoding used at runtime, where positive =
        // right shift, negative = left shift, and t1.a is applied
        // first (= rightmost factor in the matrix product).
        //   X1 = (I+L^c)(I+R^b)(I+L^a)  →  (-a, +b, -c)
        //   X2 = (I+L^a)(I+R^b)(I+L^c)  →  (-c, +b, -a)
        //   X3 = (I+R^c)(I+L^b)(I+R^a)  →  (+a, -b, +c)
        //   X5 = (I+R^b)(I+L^c)(I+L^a)  →  (-a, -c, +b)
        int t1a = 0, t1b = 0, t1c = 0;
        switch (pat) {
            case 0:  t1a = -a; t1b =  b; t1c = -c; break;  // X1
            case 1:  t1a = -c; t1b =  b; t1c = -a; break;  // X2
            case 2:  t1a =  a; t1b = -b; t1c =  c; break;  // X3
            case 3:  t1a = -a; t1b = -c; t1c =  b; break;  // X5
        }
        Params p;
        p.set_int("type", 1);
        p.set_int("w", w_);
        p.set_int("r", 1);
        p.set_int_vec("shifts", {t1a, t1b, t1c});
        return p;
    }

private:
    const int w_;
};

// ── Type 2 ─────────────────────────────────────────────────────────────
//
// Axes: m (r-1 if free, 1 if pinned) × p[0..2] (2w-1) × q[0..2] (2w-1).
// Each p/q slot is in ±[0, w-1] (size 2w-1, where 0 = "skip this
// factor" — the runtime branch in next() skips zero entries).

class MarsaXorshiftType2Enumerator final : public MarsaXorshiftEnumeratorBase {
public:
    MarsaXorshiftType2Enumerator(int w, int r, int m_user)
        : w_(w), r_(r), m_user_(m_user)
    {
        const long m_size = (m_user_ != 0) ? 1 : (r_ - 1);
        const long shift_size = 2 * w_ - 1;
        push_axis_("m", NTL::ZZ(m_size),
                   m_user_ != 0 ? std::string("m pinned by user")
                                : std::string("tap offset m ∈ [1, r-1]"));
        for (int i = 0; i < 3; ++i) {
            std::ostringstream nm; nm << "p[" << i << "]";
            push_axis_(nm.str(), NTL::ZZ(shift_size),
                       "signed shift in [-(w-1), w-1] (0 = skip)");
        }
        for (int i = 0; i < 3; ++i) {
            std::ostringstream nm; nm << "q[" << i << "]";
            push_axis_(nm.str(), NTL::ZZ(shift_size),
                       "signed shift in [-(w-1), w-1] (0 = skip)");
        }
        total_ = NTL::ZZ(m_size);
        for (int i = 0; i < 6; ++i) total_ *= NTL::ZZ(shift_size);
    }

    Params at(const std::string& idx_dec) const override {
        const NTL::ZZ idx = check_idx_(idx_dec);
        auto digits = mixed_radix_decode(sizes_, idx);
        const int m = (m_user_ != 0)
                          ? m_user_
                          : (int)NTL::to_long(digits[0]) + 1;
        std::vector<int> p(3), q(3);
        for (int i = 0; i < 3; ++i)
            p[i] = decode_signed_with_zero(w_, NTL::to_long(digits[1 + i]));
        for (int i = 0; i < 3; ++i)
            q[i] = decode_signed_with_zero(w_, NTL::to_long(digits[4 + i]));
        Params out;
        out.set_int("type", 2);
        out.set_int("w", w_);
        out.set_int("r", r_);
        out.set_int("m", m);
        out.set_int_vec("p", p);
        out.set_int_vec("q", q);
        return out;
    }

private:
    const int w_, r_, m_user_;
};

// ── Type 3 ─────────────────────────────────────────────────────────────
//
// Axes: tap-subset (C(r, 3)) × shift[0..2] (2(w-1)).
// Choose 3 distinct tap positions from {1, ..., r}; for each, pick a
// nonzero signed shift in ±[1, w-1].
//
// Behaviour at the boundary: when the enumerator picks the last tap
// at position == r, the runtime's cross-word read in next() (type 3)
// skips the high-32 OR (guarded by `idx + 1 < r`) so that tap reads
// only the last w-bit block — a regular xorshift, no cross-word
// extension.  See MarsaXorshiftGen::next() for the read pattern.

class MarsaXorshiftType3Enumerator final : public MarsaXorshiftEnumeratorBase {
public:
    MarsaXorshiftType3Enumerator(int w, int r) : w_(w), r_(r) {
        const NTL::ZZ tap_count = binomial_zz(r_, 3);
        const long shift_size = 2 * (w_ - 1);
        std::ostringstream tap_describe;
        tap_describe << "3-subset of {1, ..., r=" << r_
                     << "} (C(r, 3) = " << tap_count
                     << ").  Note: a tap at position r falls back to a "
                        "non-cross-word read in the runtime; see "
                        "next() guard `idx + 1 < r`.";
        push_axis_("taps", tap_count, tap_describe.str());
        for (int i = 0; i < 3; ++i) {
            std::ostringstream nm; nm << "shift[" << i << "]";
            push_axis_(nm.str(), NTL::ZZ(shift_size),
                       "signed nonzero shift in ±[1, w-1]");
        }
        total_ = tap_count
               * NTL::ZZ(shift_size) * NTL::ZZ(shift_size)
               * NTL::ZZ(shift_size);
    }

    Params at(const std::string& idx_dec) const override {
        const NTL::ZZ idx = check_idx_(idx_dec);
        auto digits = mixed_radix_decode(sizes_, idx);
        // unrank_combination returns ascending 0-based positions; the
        // runtime uses 1-based block indices, so add 1.
        std::vector<int> tap_idx = unrank_combination(r_, 3, digits[0]);
        std::vector<int> tap_positions(3), tap_shifts(3);
        for (int i = 0; i < 3; ++i) {
            tap_positions[i] = tap_idx[i] + 1;
            tap_shifts[i] = decode_signed_nonzero(w_,
                NTL::to_long(digits[1 + i]));
        }
        Params out;
        out.set_int("type", 3);
        out.set_int("w", w_);
        out.set_int("r", r_);
        out.set_int_vec("tap_positions", tap_positions);
        out.set_int_vec("tap_shifts", tap_shifts);
        return out;
    }

private:
    const int w_, r_;
};

// ── Type 4 ─────────────────────────────────────────────────────────────
//
// Axes: m (r-1 if free, 1 if pinned) × p[0..1] (2w-1) × q[0..1] (2w-1).

class MarsaXorshiftType4Enumerator final : public MarsaXorshiftEnumeratorBase {
public:
    MarsaXorshiftType4Enumerator(int w, int r, int m_user)
        : w_(w), r_(r), m_user_(m_user)
    {
        const long m_size = (m_user_ != 0) ? 1 : (r_ - 1);
        const long shift_size = 2 * w_ - 1;
        push_axis_("m", NTL::ZZ(m_size),
                   m_user_ != 0 ? std::string("m pinned by user")
                                : std::string("tap offset m ∈ [1, r-1]"));
        for (int i = 0; i < 2; ++i) {
            std::ostringstream nm; nm << "p[" << i << "]";
            push_axis_(nm.str(), NTL::ZZ(shift_size),
                       "signed shift in [-(w-1), w-1]");
        }
        for (int i = 0; i < 2; ++i) {
            std::ostringstream nm; nm << "q[" << i << "]";
            push_axis_(nm.str(), NTL::ZZ(shift_size),
                       "signed shift in [-(w-1), w-1]");
        }
        total_ = NTL::ZZ(m_size);
        for (int i = 0; i < 4; ++i) total_ *= NTL::ZZ(shift_size);
    }

    Params at(const std::string& idx_dec) const override {
        const NTL::ZZ idx = check_idx_(idx_dec);
        auto digits = mixed_radix_decode(sizes_, idx);
        const int m = (m_user_ != 0)
                          ? m_user_
                          : (int)NTL::to_long(digits[0]) + 1;
        std::vector<int> p(2), q(2);
        for (int i = 0; i < 2; ++i)
            p[i] = decode_signed_with_zero(w_, NTL::to_long(digits[1 + i]));
        for (int i = 0; i < 2; ++i)
            q[i] = decode_signed_with_zero(w_, NTL::to_long(digits[3 + i]));
        Params out;
        out.set_int("type", 4);
        out.set_int("w", w_);
        out.set_int("r", r_);
        out.set_int("m", m);
        out.set_int_vec("p", p);
        out.set_int_vec("q", q);
        return out;
    }

private:
    const int w_, r_, m_user_;
};

// ── Type 100 ───────────────────────────────────────────────────────────
//
// General multi-tap: vₙ = Σⱼ (I + S^{aⱼ,1})...(I + S^{aⱼ,kⱼ}) v_{n−mⱼ}
// where the j-th tap is at position mⱼ ∈ {1, ..., r} and carries
// shifts_per_tap_[j] signed nonzero shifts in ±[1, w-1].
//
// Axes (outermost first):
//   tap-subset (size C(r, nb_taps_))
//   then for each tap j ∈ [0, nb_taps_):
//     shift slots j,0 .. j,(shifts_per_tap_[j]-1)  — each size 2(w-1)
//
// Total = C(r, nb_taps_) × Π_j (2(w-1))^{shifts_per_tap_[j]}.

class MarsaXorshiftType100Enumerator final : public MarsaXorshiftEnumeratorBase {
public:
    MarsaXorshiftType100Enumerator(int w, int r, int nb_taps,
                                   std::vector<int> shifts_per_tap)
        : w_(w), r_(r), nb_taps_(nb_taps),
          shifts_per_tap_(std::move(shifts_per_tap))
    {
        const NTL::ZZ tap_count = binomial_zz(r_, nb_taps_);
        const long shift_size = 2 * (w_ - 1);
        std::ostringstream tap_describe;
        tap_describe << nb_taps_ << "-subset of {1, ..., r=" << r_
                     << "} (C(r, " << nb_taps_ << ") = " << tap_count << ")";
        push_axis_("taps", tap_count, tap_describe.str());

        total_ = tap_count;
        for (int j = 0; j < nb_taps_; ++j) {
            for (int s = 0; s < shifts_per_tap_[j]; ++s) {
                std::ostringstream nm;
                nm << "tap[" << j << "].shift[" << s << "]";
                push_axis_(nm.str(), NTL::ZZ(shift_size),
                           "signed nonzero shift in ±[1, w-1]");
                total_ *= NTL::ZZ(shift_size);
            }
        }
    }

    Params at(const std::string& idx_dec) const override {
        const NTL::ZZ idx = check_idx_(idx_dec);
        auto digits = mixed_radix_decode(sizes_, idx);
        // Outermost digit selects the position subset.
        std::vector<int> tap_idx = unrank_combination(r_, nb_taps_, digits[0]);
        std::vector<int> mi_positions(nb_taps_);
        std::vector<int> mi_counts(nb_taps_);
        std::vector<int> mi_shifts;
        size_t total_shifts = 0;
        for (int s : shifts_per_tap_) total_shifts += (size_t)s;
        mi_shifts.reserve(total_shifts);

        size_t digit_cursor = 1;   // skip the tap-subset digit at index 0
        for (int j = 0; j < nb_taps_; ++j) {
            mi_positions[j] = tap_idx[j] + 1;     // 0-based → 1-based
            mi_counts[j] = shifts_per_tap_[j];
            for (int s = 0; s < shifts_per_tap_[j]; ++s) {
                mi_shifts.push_back(decode_signed_nonzero(
                    w_, NTL::to_long(digits[digit_cursor++])));
            }
        }
        Params out;
        out.set_int("type", 100);
        out.set_int("w", w_);
        out.set_int("r", r_);
        out.set_int_vec("mi_positions", mi_positions);
        out.set_int_vec("mi_counts", mi_counts);
        out.set_int_vec("mi_shifts", mi_shifts);
        return out;
    }

private:
    const int w_, r_, nb_taps_;
    const std::vector<int> shifts_per_tap_;
};

}  // anonymous namespace

// ── Factory: validate up front, dispatch to the right subclass ─────────
//
// All structural validation happens here — once we hand off to a
// subclass ctor, the subclass trusts its inputs (matches the
// TauswortheGen::make_enumerator pattern).  The L parameter is
// unused: the recurrence's full period 2^(wr) - 1 is fixed by (w, r)
// alone and has no quicktaus-style L coupling.
std::unique_ptr<GenEnumerator> MarsaXorshiftGen::make_enumerator(
    const Params& resolved, int /*L*/)
{
    if (!resolved.has("type"))
        throw std::invalid_argument("needs_type");
    if (!resolved.has("w"))
        throw std::invalid_argument("needs_w");
    const int type = (int)resolved.get_int("type");
    const int w    = (int)resolved.get_int("w");
    if (w < 2)
        throw std::invalid_argument("needs_w_ge_2");
    // No upper bound: the runtime's wide path (BitVect-based) handles
    // w > 64 transparently; the enumerator's mixed-radix axes and
    // unrank_combination already use NTL::ZZ so size scales fine.

    // ── Type 1: r is implicitly 1 ────────────────────────────────────
    if (type == 1) {
        if (resolved.has("r")) {
            const int r = (int)resolved.get_int("r");
            if (r != 1)
                throw std::invalid_argument("needs_r_eq_1_for_type1");
        }
        return std::make_unique<MarsaXorshiftType1Enumerator>(w);
    }

    // ── Types 2, 3, 4, 100: r is required ────────────────────────────
    if (!resolved.has("r"))
        throw std::invalid_argument("needs_r");
    const int r = (int)resolved.get_int("r");

    auto read_optional_m = [&]() -> int {
        // 0 = "enumerate m"; >= 1 = "user pinned m to this value".
        // We reject a present-but-invalid m (e.g. m=0) instead of
        // silently falling back to the enumerated form — the user's
        // intent is unclear and silent fallbacks are surprising.
        if (!resolved.has("m")) return 0;
        const long m_val = resolved.get_int("m");
        if (m_val < 1 || m_val >= r)
            throw std::invalid_argument("needs_admissible_fixed_m");
        return (int)m_val;
    };

    if (type == 2) {
        if (r < 2) throw std::invalid_argument("needs_r_ge_2_for_type2");
        return std::make_unique<MarsaXorshiftType2Enumerator>(
            w, r, read_optional_m());
    }

    if (type == 3) {
        if (r < 3) throw std::invalid_argument("needs_r_ge_3_for_type3");
        return std::make_unique<MarsaXorshiftType3Enumerator>(w, r);
    }

    if (type == 4) {
        if (r < 2) throw std::invalid_argument("needs_r_ge_2_for_type4");
        return std::make_unique<MarsaXorshiftType4Enumerator>(
            w, r, read_optional_m());
    }

    if (type == 100) {
        // nb_taps defaults to 3; if user provides it, it must be ≥ 1.
        int nb_taps = 3;
        if (resolved.has("nb_taps")) {
            const long v = resolved.get_int("nb_taps");
            if (v < 1)
                throw std::invalid_argument("needs_nb_taps_ge_1");
            nb_taps = (int)v;
        }
        if (r < nb_taps)
            throw std::invalid_argument("needs_r_ge_nb_taps");

        // shifts_per_tap defaults to a length-nb_taps vector of 1s
        // (Table IV form).  If user provides one, it must match
        // nb_taps in length and have every entry ≥ 1.
        std::vector<int> shifts_per_tap;
        if (resolved.has("shifts_per_tap"))
            shifts_per_tap = resolved.get_int_vec("shifts_per_tap");
        if (shifts_per_tap.empty()) {
            shifts_per_tap.assign((size_t)nb_taps, 1);
        } else {
            if ((int)shifts_per_tap.size() != nb_taps)
                throw std::invalid_argument(
                    "needs_shifts_per_tap_len_eq_nb_taps");
            for (int s : shifts_per_tap)
                if (s < 1)
                    throw std::invalid_argument("needs_shifts_per_tap_ge_1");
        }
        return std::make_unique<MarsaXorshiftType100Enumerator>(
            w, r, nb_taps, std::move(shifts_per_tap));
    }

    throw std::invalid_argument(
        "MarsaXorshiftGen::make_enumerator: unsupported type (must be "
        "one of 1, 2, 3, 4, 100)");
}

}  // namespace regpoly::core
