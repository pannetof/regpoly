#include "gen_matsumoto.h"
#include <cstdio>

static int compute_k(int type, int n) {
    if (type == 1 || type == 2)
        return (n + 1) * 32;
    else  // type == 3
        return (n + 3) * 32;
}

MatsumotoGen::MatsumotoGen(int type, int n, int m,
                           const std::vector<int>& paramsint,
                           const std::vector<uint32_t>& paramsunsigned,
                           int L)
    : Generateur(compute_k(type, n), L),
      type_(type), n_(n), m_(m), p_(paramsint), pu_(paramsunsigned)
{
    state_ = BitVect(k_);
}

std::string MatsumotoGen::name() const { return "Matsumoto Generator"; }

std::string MatsumotoGen::display_str() const {
    char buf[512];
    int offset = snprintf(buf, sizeof(buf), "ran%d -> n = %d, m = %d", type_, n_, m_);
    for (size_t i = 0; i < p_.size(); i++) {
        offset += snprintf(buf + offset, sizeof(buf) - offset, ", p%d = %d", (int)i, p_[i]);
    }
    for (size_t i = 0; i < pu_.size(); i++) {
        offset += snprintf(buf + offset, sizeof(buf) - offset, ", pu%d = %08x", (int)i, pu_[i]);
    }
    return std::string(buf);
}

void MatsumotoGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

uint32_t MatsumotoGen::SHIFT(uint32_t v, int i) {
    if (i > 0)
        return v >> i;
    else
        return v << (-i);
}

uint32_t MatsumotoGen::V(int idx) const {
    return (uint32_t)state_.get_word(idx, 32);
}

void MatsumotoGen::SetV(int idx, uint32_t val) {
    state_.set_word(idx, 32, val & 0xFFFFFFFFu);
}

void MatsumotoGen::next() {
    if (type_ == 1) {
        uint32_t z0 = V(n_) ^ SHIFT(V(0), p_[0]) ^ SHIFT(V(0), p_[1]);
        uint32_t z1 = z0 ^ V(m_);
        uint32_t z2 = z1 ^ SHIFT(z1, p_[2]);
        uint32_t z3 = V(0) ^ z2 ^ SHIFT(z2, p_[3]);
        state_.lshift(32);
        SetV(n_ - 1, z3);
        SetV(n_, z2);

    } else if (type_ == 2) {
        uint32_t z0 = V(m_);
        uint32_t z1;
        if (V(n_) & 1) {
            z1 = (V(n_) >> 1) ^ SHIFT(z0, p_[0]) ^ SHIFT(z0, p_[1]) ^ pu_[0];
        } else {
            z1 = (V(n_) >> 1) ^ SHIFT(z0, p_[0]) ^ SHIFT(z0, p_[1]);
        }
        uint32_t z2 = V(0);
        uint32_t z3 = z2 ^ SHIFT(z2, p_[2]);
        uint32_t z4 = z3 ^ SHIFT(z3, p_[3]);
        state_.lshift(32);
        SetV(n_ - 1, z1);
        SetV(n_, z4);

    } else if (type_ == 3) {
        uint32_t v0 = V(0);
        uint32_t z0 = V(n_) ^ v0 ^ SHIFT(V(m_), p_[0]);
        uint32_t z1 = V(n_ + 1) ^ z0 ^ SHIFT(z0, p_[1]);
        uint32_t z2 = V(n_ + 2) ^ z1 ^ SHIFT(z1, p_[2]);
        state_.lshift(32);
        SetV(n_ - 1, v0 ^ z0 ^ z1 ^ z2);
        SetV(n_, z0);
        SetV(n_ + 1, z1);
        SetV(n_ + 2, z2);
    }
}

std::unique_ptr<Generateur> MatsumotoGen::copy() const {
    auto g = std::make_unique<MatsumotoGen>(type_, n_, m_, p_, pu_, L_);
    g->state_ = state_.copy();
    return g;
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generateur> MatsumotoGen::from_params(const Params& params, int L) {
    int type = (int)params.get_int("type");
    int n = (int)params.get_int("n");
    int m = (int)params.get_int("m");
    auto paramsint = params.get_int_vec("paramsint");
    auto paramsunsigned_u64 = params.get_uint_vec("paramsunsigned");
    std::vector<uint32_t> paramsunsigned;
    for (auto v : paramsunsigned_u64)
        paramsunsigned.push_back((uint32_t)v);
    return std::make_unique<MatsumotoGen>(type, n, m, paramsint, paramsunsigned, L);
}

std::vector<ParamSpec> MatsumotoGen::param_specs() {
    return {
        {"type",           "int",      true,  false, 0, "",     "", false},
        {"n",              "int",      true,  false, 0, "",     "", false},
        {"m",              "int",      false, false, 0, "range","1,n-1", false},
        {"paramsint",      "int_vec",  false, false, 0, "none", "", false},
        {"paramsunsigned", "uint_vec", false, false, 0, "none", "", false},
    };
}
