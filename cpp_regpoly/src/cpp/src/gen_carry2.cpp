#include "gen_carry2.h"
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <iomanip>

Carry2Gen::Carry2Gen(int w, int r, int p, int m1, int m2, int m3,
                     const std::vector<MatrixEntry>& matrices, int L)
    : Generateur(w * r - p, L),
      w_(w), r_(r), p_(p), m1_(m1), m2_(m2), m3_(m3),
      matrices_(matrices), i_(0), state_bits_(w * r), maskp_(0), umaskp_(0)
{
    state_ = BitVect(state_bits_);
}

std::string Carry2Gen::name() const { return "Carry Generator"; }

// ── Display matching POL output ─────────────────────────────────────────

int Carry2Gen::type_cost(int type) {
    static const int costs[] = {3, 1, 5, 2, 4, 8, 7, 0};
    return (type >= 0 && type < 8) ? costs[type] : 0;
}

std::string Carry2Gen::type_display(int type, const int* pi, const uint64_t* pu) {
    std::ostringstream oss;
    switch (type) {
        case 0: oss << "T0(" << pi[0] << ")"; break;
        case 1: oss << "Identity"; break;
        case 2: oss << "T2(" << std::hex << std::setfill('0') << std::setw(8)
                    << (unsigned)(pu[0] & M32) << std::dec << ")"; break;
        case 3: oss << "T3(" << pi[0] << ")"; break;
        case 4: oss << "T4(" << pi[0] << "," << std::hex << std::setfill('0')
                    << std::setw(8) << (unsigned)(pu[0] & M32) << std::dec << ")"; break;
        case 5: oss << "T5(" << pi[0] << ","
                    << std::hex << std::setfill('0')
                    << std::setw(8) << (unsigned)(pu[0] & M32) << ","
                    << std::setw(8) << (unsigned)(pu[1] & M32) << ","
                    << std::setw(8) << (unsigned)(pu[2] & M32) << std::dec << ")"; break;
        case 6: oss << "T6(" << pi[0] << "," << pi[1] << "," << pi[2] << ")"; break;
        case 7: oss << "ZERO"; break;
        default: oss << "Unknown"; break;
    }
    return oss.str();
}

std::string Carry2Gen::display_str() const {
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
    for (int j = 0; j < 8 && j < (int)matrices_.size(); j++) {
        oss << "\nA_" << j << " = "
            << type_display(matrices_[j].type,
                            matrices_[j].paramsint,
                            matrices_[j].paramsulong);
        cost += type_cost(matrices_[j].type);
    }
    oss << "\nCost = " << cost;
    char buf[32];
    snprintf(buf, sizeof(buf), "%3d %3d %3d", m1_, m2_, m3_);
    oss << "\n" << buf;
    for (int j = 0; j < 8 && j < (int)matrices_.size(); j++) {
        snprintf(buf, sizeof(buf), " %3d", matrices_[j].type);
        oss << buf;
        for (int i = 0; i < 3; i++) {
            snprintf(buf, sizeof(buf), " %3d", matrices_[j].paramsint[i]);
            oss << buf;
        }
        for (int i = 0; i < 3; i++) {
            snprintf(buf, sizeof(buf), " %08x",
                     (unsigned)(matrices_[j].paramsulong[i] & M32));
            oss << buf;
        }
    }
    return oss.str();
}

// ── Core operations (64-bit types, 32-bit semantics via & M32) ──────────

void Carry2Gen::init(const BitVect& init_bv) {
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
uint64_t Carry2Gen::ShiftR(uint64_t v, int s) {
    v &= M32;
    if (s > 0)
        return v >> s;
    else
        return (v << (-s)) & M32;
}

// apply_matrix: all arithmetic in 32-bit semantics.
uint64_t Carry2Gen::apply_matrix(int type, uint64_t v, const int* pi, const uint64_t* pu) {
    v &= M32;
    switch (type) {
        case 0:
            return (v ^ ShiftR(v, pi[0])) & M32;
        case 1:
            return v;
        case 2:
            return (v & 1) ? (((v >> 1) ^ pu[0]) & M32) : (v >> 1);
        case 3:
            return ShiftR(v, pi[0]);
        case 4:
            return (v ^ (ShiftR(v, pi[0]) & pu[0])) & M32;
        case 5: {
            uint64_t cond = v & pu[2];
            uint64_t rot = (((v << pi[0]) | (v >> (32 - pi[0]))) & pu[1]) & M32;
            return (rot ^ (cond ? pu[0] : 0ULL)) & M32;
        }
        case 6:
            return (v ^ ShiftR(v, pi[0]) ^ ShiftR(v, pi[1]) ^ ShiftR(v, pi[2])) & M32;
        case 7:
            return 0;
        default:
            return v;
    }
}

uint64_t Carry2Gen::TMAT(int j, uint64_t val) const {
    return apply_matrix(matrices_[j].type, val,
                        matrices_[j].paramsint, matrices_[j].paramsulong);
}

uint64_t Carry2Gen::V(int idx) const {
    return (uint64_t)state_.get_word(idx, 32);
}

void Carry2Gen::SetV(int idx, uint64_t val) {
    state_.set_word(idx, 32, val & M32);
}

void Carry2Gen::next() {
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

BitVect Carry2Gen::get_output() const {
    BitVect out(L_);
    int nw = (L_ + 31) / 32;
    for (int j = 0; j < nw; j++) {
        uint64_t word = V((i_ + j) % r_);
        out.set_word(j, 32, word);
    }
    return out;
}

void Carry2Gen::get_transition_state(uint64_t* out_words, int out_nwords) const {
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

std::unique_ptr<Generateur> Carry2Gen::copy() const {
    auto g = std::make_unique<Carry2Gen>(w_, r_, p_, m1_, m2_, m3_, matrices_, L_);
    g->state_ = state_.copy();
    g->i_ = i_;
    g->state_bits_ = state_bits_;
    g->maskp_ = maskp_;
    g->umaskp_ = umaskp_;
    return g;
}
