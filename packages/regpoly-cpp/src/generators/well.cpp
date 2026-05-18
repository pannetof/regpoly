// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "well.h"
#include <algorithm>
#include <array>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <stdexcept>

using namespace regpoly::core;


namespace regpoly::core {

WELLGen::WELLGen(int w, int r, int p, int m1, int m2, int m3,
                     const std::vector<MatrixEntry>& matrices, int L)
    : Generator(w * r - p, L),
      w_(w), r_(r), p_(p), m1_(m1), m2_(m2), m3_(m3),
      matrices_(matrices), i_(0), state_bits_(w * r), maskp_(0), umaskp_(0)
{
    state_ = BitVect(state_bits_);
}

std::string WELLGen::name() const { return "Carry Generator"; }

// ── Display matching POL output ─────────────────────────────────────────

int WELLGen::type_cost(int Mi) {
    return static_cost_for_Mi(Mi);
}

int WELLGen::static_cost_for_Mi(int Mi) {
    // Costs indexed by paper M-class (M0..M6).
    // Intentionally non-monotonic: M4 (conditional XOR) is dearer than
    // M5 (masked shift) on most architectures.
    static const int costs[7] = {0, 1, 2, 3, 5, 4, 8};
    return (Mi >= 0 && Mi < 7) ? costs[Mi] : 0;
}

int WELLGen::total_cost() const {
    int sum = 0;
    for (const auto& m : matrices_) sum += static_cost_for_Mi(m.Mi);
    return sum;
}

std::string WELLGen::type_display(const MatrixEntry& m) {
    std::ostringstream oss;
    auto hex8 = [&](uint32_t v) {
        oss << "0x" << std::hex << std::setfill('0') << std::setw(8)
            << v << std::dec;
    };
    switch (m.Mi) {
        case 0: oss << "M0"; break;
        case 1: oss << "M1"; break;
        case 2: oss << "M2(" << m.t << ")"; break;
        case 3: oss << "M3(" << m.t << ")"; break;
        case 4: oss << "M4("; hex8(m.a); oss << ")"; break;
        case 5: oss << "M5(" << m.t << ", "; hex8(m.b); oss << ")"; break;
        case 6: oss << "M6(" << m.q << ", " << m.t << ", " << m.s << ", ";
                hex8(m.a); oss << ")"; break;
        default: oss << "Unknown(" << m.Mi << ")"; break;
    }
    return oss.str();
}

std::string WELLGen::display_str() const {
    std::ostringstream oss;
    oss << k_ << "\n";
    oss << " w= " << std::setw(3) << w_
        << "  r=" << std::setw(3) << r_
        << "  p= " << std::setw(3) << p_
        << "  m1=" << std::setw(3) << m1_
        << "  m2=" << std::setw(3) << m2_
        << "  m3=" << std::setw(3) << m3_
        << "  wordno= " << std::setw(3) << 0;
    int cost = 0;
    // 8 algorithm slots in the WELL recurrence (T0..T7); each slot's
    // class is one of M0..M6 from paper Table I.
    for (int j = 0; j < 8 && j < (int)matrices_.size(); j++) {
        oss << "\nA_" << j << " = " << type_display(matrices_[j]);
        cost += type_cost(matrices_[j].Mi);
    }
    oss << "\nCost = " << cost;
    return oss.str();
}

// ── Core operations (32-bit semantics via & M32) ────────────────────────

void WELLGen::init(const BitVect& init_bv) {
    state_ = BitVect(state_bits_);
    state_.copy_part_from(init_bv, k_);
    i_ = 0;

    if (p_ > 0 && p_ < 32)
        maskp_ = (1ULL << p_) - 1;
    else if (p_ >= 32)
        maskp_ = M32;
    else
        maskp_ = 0;
    umaskp_ = (~maskp_) & M32;
}

// ShiftR: 32-bit shift matching C's uint32_t behavior.
// Result is always masked to 32 bits.
uint64_t WELLGen::ShiftR(uint64_t v, int s) {
    v &= M32;
    if (s > 0)
        return v >> s;
    else
        return (v << (-s)) & M32;
}

// apply_matrix: paper Table I. All arithmetic in 32-bit semantics.
//
// `m.Mi` is a paper Mi class index (0..6). M6's d_s and test masks are
// synthesised from the small integers `s`, `t` at runtime — matching
// the paper exactly. The default branch preserves identity-fallback so
// search loops with stale candidates don't crash inside next();
// from_params() is the fail-fast gate at construction time.
uint64_t WELLGen::apply_matrix(const MatrixEntry& m, uint64_t v) {
    v &= M32;
    switch (m.Mi) {
        case 0:  // M0 — zero
            return 0;
        case 1:  // M1 — identity
            return v;
        case 2:  // M2(t) — shift
            return ShiftR(v, m.t);
        case 3:  // M3(t) — x ⊕ shift(x, t)
            return (v ^ ShiftR(v, m.t)) & M32;
        case 4:  // M4(a) — MT-style: (v>>1) ⊕ a if LSB(v) else (v>>1)
            return (v & 1ULL) ? (((v >> 1) ^ m.a) & M32) : (v >> 1);
        case 5:  // M5(t, b) — x ⊕ (shift(x, t) & b)
            return (v ^ (ShiftR(v, m.t) & m.b)) & M32;
        case 6: {  // M6(q, t, s, a) — rotate, mask, conditional XOR
            // Paper d_s has the (s+1)-th bit zero (1-indexed from MSB),
            // i.e. ~(1 << (31 - s)). Test bit at MSB-position t is
            // (1 << (31 - t)).
            uint32_t d_s   = ~(1u << (31 - m.s));
            uint32_t test  =  (1u << (31 - m.t));
            uint32_t cond  =  static_cast<uint32_t>(v & test);
            uint32_t rot   =  static_cast<uint32_t>(
                                  ((v << m.q) | (v >> (32 - m.q))) & M32);
            return ((rot & d_s) ^ (cond ? m.a : 0u)) & M32;
        }
        default:
            return v;
    }
}

uint64_t WELLGen::TMAT(int j, uint64_t val) const {
    return apply_matrix(matrices_[j], val);
}

uint64_t WELLGen::V(int idx) const {
    return (uint64_t)state_.get_word(idx, 32);
}

void WELLGen::SetV(int idx, uint64_t val) {
    state_.set_word(idx, 32, val & M32);
}

void WELLGen::next() {
    uint64_t z0 = (V((i_ + r_ - 1) % r_) & umaskp_)
                | (V((i_ + r_ - 2) % r_) & maskp_);
    uint64_t z1 = TMAT(0, V(i_)) ^ TMAT(1, V((i_ + m1_) % r_));
    uint64_t z2 = TMAT(2, V((i_ + m2_) % r_)) ^ TMAT(3, V((i_ + m3_) % r_));
    uint64_t z3 = z1 ^ z2;
    uint64_t z4 = TMAT(4, z0) ^ TMAT(5, z1) ^ TMAT(6, z2) ^ TMAT(7, z3);

    SetV((i_ + r_ - 1) % r_, z4);
    SetV(i_, z3);

    i_ = (i_ + r_ - 1) % r_;
}

BitVect WELLGen::get_output() const {
    BitVect out(L_);
    int nw = (L_ + 31) / 32;
    for (int j = 0; j < nw; j++) {
        uint64_t word = V((i_ + j) % r_);
        out.set_word(j, 32, word);
    }
    return out;
}

void WELLGen::get_transition_state(uint64_t* out_words, int out_nwords) const {
    BitVect tmp(k_);
    int nw = (k_ + 31) / 32;
    for (int j = 0; j < nw; j++) {
        uint64_t word = V((i_ + j) % r_);
        tmp.set_word(j, 32, word);
    }
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

std::unique_ptr<Generator> WELLGen::copy() const {
    auto g = std::make_unique<WELLGen>(w_, r_, p_, m1_, m2_, m3_, matrices_, L_);
    g->state_ = state_.copy();
    g->i_ = i_;
    g->state_bits_ = state_bits_;
    g->maskp_ = maskp_;
    g->umaskp_ = umaskp_;
    return g;
}

// ── Factory methods ────────────────────────────────────────────────────

namespace {

// Decode an `M:` value (string "M0".."M6" or integer 0..6) into 0..6.
int decode_mi(const ParamScalar& v, const std::string& slot_key) {
    if (auto pi = std::get_if<int64_t>(&v)) {
        int64_t mi = *pi;
        if (mi < 0 || mi > 6)
            throw std::out_of_range(
                "WELLGen: matrices['" + slot_key + "'].M = "
                + std::to_string(mi) + " is out of range (0..6)");
        return static_cast<int>(mi);
    }
    if (auto ps = std::get_if<std::string>(&v)) {
        const std::string& s = *ps;
        if (s.size() == 2 && s[0] == 'M' && s[1] >= '0' && s[1] <= '6')
            return s[1] - '0';
        throw std::runtime_error(
            "WELLGen: matrices['" + slot_key + "'].M = '" + s
            + "' is not a paper Mi class (use 0..6 or 'M0'..'M6')");
    }
    throw std::runtime_error(
        "WELLGen: matrices['" + slot_key + "'].M must be int 0..6 or "
        "string 'M0'..'M6'");
}

// Pull a signed integer arg by name; throw if missing or wrong type.
int get_int_arg(const StructEntry& e, const std::string& arg,
                const std::string& slot_key, const std::string& mi_label) {
    auto it = e.find(arg);
    if (it == e.end())
        throw std::runtime_error(
            "WELLGen: matrices['" + slot_key + "'] (" + mi_label
            + ") is missing required arg '" + arg + "'");
    if (auto pi = std::get_if<int64_t>(&it->second))
        return static_cast<int>(*pi);
    if (auto pu = std::get_if<uint64_t>(&it->second))
        return static_cast<int>(*pu);
    throw std::runtime_error(
        "WELLGen: matrices['" + slot_key + "']." + arg
        + " must be an integer");
}

// Pull a 32-bit unsigned mask arg by name.
uint32_t get_u32_arg(const StructEntry& e, const std::string& arg,
                     const std::string& slot_key, const std::string& mi_label) {
    auto it = e.find(arg);
    if (it == e.end())
        throw std::runtime_error(
            "WELLGen: matrices['" + slot_key + "'] (" + mi_label
            + ") is missing required arg '" + arg + "'");
    uint64_t raw = 0;
    if (auto pu = std::get_if<uint64_t>(&it->second)) raw = *pu;
    else if (auto pi = std::get_if<int64_t>(&it->second)) raw = static_cast<uint64_t>(*pi);
    else
        throw std::runtime_error(
            "WELLGen: matrices['" + slot_key + "']." + arg
            + " must be a 32-bit unsigned integer");
    if (raw > 0xFFFFFFFFull)
        throw std::out_of_range(
            "WELLGen: matrices['" + slot_key + "']." + arg
            + " does not fit in 32 bits");
    return static_cast<uint32_t>(raw);
}

// Reject any arg key that isn't part of this Mi's signature.
void reject_extra_args(const StructEntry& e,
                       const std::vector<std::string>& allowed,
                       const std::string& slot_key,
                       const std::string& mi_label) {
    for (const auto& kv : e) {
        if (kv.first == "M") continue;
        if (std::find(allowed.begin(), allowed.end(), kv.first) == allowed.end())
            throw std::runtime_error(
                "WELLGen: matrices['" + slot_key + "'] (" + mi_label
                + ") has unexpected arg '" + kv.first + "' (allowed: "
                + (allowed.empty() ? std::string("none") : [&](){
                      std::string out;
                      for (size_t i = 0; i < allowed.size(); ++i) {
                          if (i) out += ", ";
                          out += allowed[i];
                      }
                      return out;
                  }()) + ")");
    }
}

WELLGen::MatrixEntry decode_matrix_entry(const std::string& slot_key,
                                          const StructEntry& e) {
    auto m_it = e.find("M");
    if (m_it == e.end())
        throw std::runtime_error(
            "WELLGen: matrices['" + slot_key + "'] is missing 'M' key");
    int Mi = decode_mi(m_it->second, slot_key);
    std::string label = "M" + std::to_string(Mi);

    WELLGen::MatrixEntry out;
    out.Mi = Mi;
    switch (Mi) {
        case 0:
        case 1:
            reject_extra_args(e, {}, slot_key, label);
            break;
        case 2:
        case 3:
            out.t = get_int_arg(e, "t", slot_key, label);
            if (out.t < -32 || out.t > 32)
                throw std::out_of_range(
                    "WELLGen: matrices['" + slot_key + "'].t = "
                    + std::to_string(out.t) + " is out of range [-32, 32]");
            reject_extra_args(e, {"t"}, slot_key, label);
            break;
        case 4:
            out.a = get_u32_arg(e, "a", slot_key, label);
            reject_extra_args(e, {"a"}, slot_key, label);
            break;
        case 5:
            out.t = get_int_arg(e, "t", slot_key, label);
            out.b = get_u32_arg(e, "b", slot_key, label);
            if (out.t < -32 || out.t > 32)
                throw std::out_of_range(
                    "WELLGen: matrices['" + slot_key + "'].t = "
                    + std::to_string(out.t) + " is out of range [-32, 32]");
            reject_extra_args(e, {"t", "b"}, slot_key, label);
            break;
        case 6:
            out.q = get_int_arg(e, "q", slot_key, label);
            out.t = get_int_arg(e, "t", slot_key, label);
            out.s = get_int_arg(e, "s", slot_key, label);
            out.a = get_u32_arg(e, "a", slot_key, label);
            if (out.q < 0 || out.q > 31 || out.t < 0 || out.t > 31
                || out.s < 0 || out.s > 31)
                throw std::out_of_range(
                    "WELLGen: matrices['" + slot_key + "'] M6 args q/t/s "
                    "must each be in [0, 31]");
            reject_extra_args(e, {"q", "t", "s", "a"}, slot_key, label);
            break;
    }
    return out;
}

// Reject the legacy flat-triple keys with a clear migration pointer.
void reject_legacy_triple(const Params& params) {
    static const char* keys[] = {"mat_types", "mat_pi", "mat_pu"};
    for (const char* k : keys) {
        if (params.has(k))
            throw std::runtime_error(
                std::string("WELLGen: '") + k + "' is no longer accepted. "
                "Use the structured 'matrices' map keyed by T0..T7. "
                "See docs/generators/WELLGen.md.");
    }
}

}  // namespace

std::unique_ptr<Generator> WELLGen::from_params(const Params& params, int L) {
    int w = (int)params.get_int("w", 32);
    int r = (int)params.get_int("r");
    int p = (int)params.get_int("p");
    int m1 = (int)params.get_int("m1");
    int m2 = (int)params.get_int("m2");
    int m3 = (int)params.get_int("m3");

    reject_legacy_triple(params);
    if (!params.has_struct_map("matrices"))
        throw std::runtime_error(
            "WELLGen: 'matrices' is required (a map keyed by T0..T7). "
            "See docs/generators/WELLGen.md.");

    const StructMap& m = params.get_struct_map("matrices");
    std::vector<MatrixEntry> matrices(8);
    for (int j = 0; j < 8; j++) {
        std::string key = "T" + std::to_string(j);
        auto it = m.find(key);
        if (it == m.end())
            throw std::runtime_error(
                "WELLGen: matrices is missing slot '" + key + "'");
        matrices[j] = decode_matrix_entry(key, it->second);
    }
    // Reject extra slot keys outside T0..T7.
    for (const auto& kv : m) {
        const std::string& k = kv.first;
        if (k.size() != 2 || k[0] != 'T' || k[1] < '0' || k[1] > '7')
            throw std::runtime_error(
                "WELLGen: matrices has unexpected slot key '" + k
                + "' (allowed: T0..T7)");
    }
    return std::make_unique<WELLGen>(w, r, p, m1, m2, m3, matrices, L);
}

std::vector<ParamSpec> WELLGen::param_specs() {
    return {
        {"w",         "int",        true,  true,  32, "",        "", false},
        {"r",         "int",        true,  false, 0,  "",        "", false},
        {"p",         "int",        true,  false, 0,  "",        "", false},
        {"m1",        "int",        false, false, 0,  "range",   "1,r-1", false},
        {"m2",        "int",        false, false, 0,  "range",   "1,r-1", false},
        {"m3",        "int",        false, false, 0,  "range",   "1,r-1", false},
        {"matrices",  "struct_map", false, false, 0,  "none",    "", false},
    };
}

// ── Cost-bounded matrices sampler ──────────────────────────────────────

namespace {

// Per-Mi arg samplers. Each fills the StructEntry with the named args
// expected by `decode_matrix_entry` for that Mi.
//
// Ranges (paper Table I, with degenerate args excluded so search-time
// candidates aren't trivially equivalent to a cheaper Mi):
//   M2/M3/M5: t ∈ [-(w-1), -1] ∪ [1, w-1]   (excludes t=0 collapse)
//   M4:       a ∈ [1, 2^w - 1]              (excludes a=0 half-shift)
//   M5:       b ∈ [1, 2^w - 1]
//   M6:       q,t,s ∈ [0, w-1]; a ∈ [1, 2^w - 1]
int sample_signed_shift(int w, std::mt19937_64& rng) {
    // Pick magnitude in [1, w-1], then sign.
    std::uniform_int_distribution<int> mag(1, w - 1);
    std::uniform_int_distribution<int> sign(0, 1);
    int m = mag(rng);
    return sign(rng) ? m : -m;
}

uint32_t sample_nonzero_u32(std::mt19937_64& rng) {
    std::uniform_int_distribution<uint32_t> d(1, 0xFFFFFFFFu);
    return d(rng);
}

void fill_args_for_Mi(int Mi, int w, std::mt19937_64& rng, StructEntry& e) {
    switch (Mi) {
        case 0:
        case 1:
            break;
        case 2:
        case 3:
            e["t"] = static_cast<int64_t>(sample_signed_shift(w, rng));
            break;
        case 4:
            e["a"] = static_cast<uint64_t>(sample_nonzero_u32(rng));
            break;
        case 5:
            e["t"] = static_cast<int64_t>(sample_signed_shift(w, rng));
            e["b"] = static_cast<uint64_t>(sample_nonzero_u32(rng));
            break;
        case 6: {
            std::uniform_int_distribution<int> qd(0, w - 1);
            e["q"] = static_cast<int64_t>(qd(rng));
            e["t"] = static_cast<int64_t>(qd(rng));
            e["s"] = static_cast<int64_t>(qd(rng));
            e["a"] = static_cast<uint64_t>(sample_nonzero_u32(rng));
            break;
        }
        default: break;
    }
}

}  // namespace

StructMap WELLGen::random_matrices(int w, int max_cost,
                                    std::mt19937_64& rng) {
    if (max_cost <= 0)
        throw std::invalid_argument(
            "WELLGen::random_matrices: max_cost must be > 0 (got "
            + std::to_string(max_cost) + ")");
    if (w != 32)
        throw std::invalid_argument(
            "WELLGen::random_matrices: only w=32 is supported today (got "
            + std::to_string(w) + ")");

    constexpr int kSlots = 8;
    constexpr int kMaxRejectionAttempts = 64;

    std::array<int, kSlots> Mis{};
    bool accepted = false;

    // Phase 1: rejection sampling. 8 i.i.d. uniform Mi draws from {0..6}.
    std::uniform_int_distribution<int> Mi_dist(0, 6);
    for (int attempt = 0; attempt < kMaxRejectionAttempts; ++attempt) {
        int total = 0;
        for (int j = 0; j < kSlots; ++j) {
            Mis[j] = Mi_dist(rng);
            total += static_cost_for_Mi(Mis[j]);
        }
        if (total <= max_cost) { accepted = true; break; }
    }

    // Phase 2: greedy-budgeted fallback. Shuffle slot order, then for
    // each slot pick Mi uniformly from those whose cost ≤ remaining
    // budget. Always succeeds because cost(M0) = 0 ≤ any non-negative
    // remaining budget.
    if (!accepted) {
        std::array<int, kSlots> order{0, 1, 2, 3, 4, 5, 6, 7};
        std::shuffle(order.begin(), order.end(), rng);
        int remaining = max_cost;
        for (int slot : order) {
            std::vector<int> allowed;
            allowed.reserve(7);
            for (int Mi = 0; Mi < 7; ++Mi) {
                if (static_cost_for_Mi(Mi) <= remaining)
                    allowed.push_back(Mi);
            }
            std::uniform_int_distribution<size_t> pick(0, allowed.size() - 1);
            int Mi = allowed[pick(rng)];
            Mis[slot] = Mi;
            remaining -= static_cost_for_Mi(Mi);
        }
    }

    StructMap out;
    for (int j = 0; j < kSlots; ++j) {
        StructEntry e;
        e["M"] = static_cast<int64_t>(Mis[j]);
        fill_args_for_Mi(Mis[j], w, rng, e);
        out["T" + std::to_string(j)] = std::move(e);
    }
    return out;
}

}  // namespace regpoly::core
