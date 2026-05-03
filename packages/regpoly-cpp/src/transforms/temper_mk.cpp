#include "temper_mk.h"
#include <sstream>
#include <iomanip>

TemperMKTrans::TemperMKTrans(int w, int type, int eta, int mu,
                             int u, int l, uint64_t b, uint64_t c)
    : type_(type), eta_(eta), mu_(mu), u_(u), l_(l), b_(b), c_(c)
{
    w_ = w;
}

std::string TemperMKTrans::name() const {
    if (type_ == 2)
        return "Matsumoto-Kurita Tempering(II)";
    return "Matsumoto-Kurita(I) Tempering";
}

std::string TemperMKTrans::display_str() const {
    int hex_digits = (w_ + 3) / 4;
    std::ostringstream oss;
    if (type_ == 1) {
        oss << "Matsumoto-Kurita Tempering (I)("
            << eta_ << "," << mu_ << ")(w=" << w_ << ")  b = "
            << std::hex << std::setfill('0') << std::setw(hex_digits) << b_
            << "  c = "
            << std::setfill('0') << std::setw(hex_digits) << c_;
    } else {
        oss << "Matsumoto-Nishimura Tempering (II)("
            << u_ << "," << eta_ << "," << mu_ << "," << l_ << ")(w=" << w_ << ")  b = "
            << std::hex << std::setfill('0') << std::setw(hex_digits) << b_
            << "  c = "
            << std::setfill('0') << std::setw(hex_digits) << c_;
    }
    return oss.str();
}

void TemperMKTrans::apply(BitVect& state) const {
    int n = state.nbits();

    BitVect b_n(n);
    BitVect c_n(n);
    for (int i = 0; i < w_; i++) {
        if ((b_ >> (w_ - 1 - i)) & 1ULL)
            b_n.set_bit(i, 1);
    }
    for (int i = 0; i < w_; i++) {
        if ((c_ >> (w_ - 1 - i)) & 1ULL)
            c_n.set_bit(i, 1);
    }

    BitVect state_w = state.copy();
    state_w.and_mask(w_);

    BitVect rest = state.copy();
    rest.and_invmask(w_);

    if (type_ == 2) {
        BitVect shifted = state_w.copy();
        shifted.rshift(u_);
        state_w.xor_with(shifted);
        state_w.and_mask(w_);
    }

    {
        BitVect shifted = state_w.copy();
        shifted.lshift(eta_);
        for (int i = 0; i < shifted.nwords(); i++)
            shifted.data()[i] &= b_n.data()[i];
        state_w.xor_with(shifted);
    }

    {
        BitVect shifted = state_w.copy();
        shifted.lshift(mu_);
        for (int i = 0; i < shifted.nwords(); i++)
            shifted.data()[i] &= c_n.data()[i];
        state_w.xor_with(shifted);
    }

    if (type_ == 2) {
        BitVect shifted = state_w.copy();
        shifted.rshift(l_);
        state_w.xor_with(shifted);
        state_w.and_mask(w_);
    }

    rest.xor_with(state_w);
    state = rest;
}

std::unique_ptr<Transformation> TemperMKTrans::copy() const {
    return std::make_unique<TemperMKTrans>(w_, type_, eta_, mu_, u_, l_, b_, c_);
}

void TemperMKTrans::update(const Params& params) {
    w_ = (int)params.get_int("w", w_);
    b_ = (uint64_t)params.get_int("b", (int64_t)b_);
    c_ = (uint64_t)params.get_int("c", (int64_t)c_);
}

// ── Factory method ─────────────────────────────────────────────────────

std::unique_ptr<Transformation> TemperMKTrans::from_params(
    const std::string& type_name, const Params& params)
{
    int w = (int)params.get_int("w");
    int mk_type = (type_name == "tempMK2") ? 2 : 1;
    int eta = (int)params.get_int("eta");
    int mu = (int)params.get_int("mu");
    int u = (int)params.get_int("u", 0);
    int l = (int)params.get_int("l", 0);
    uint64_t b = (uint64_t)params.get_int("b", 0);
    uint64_t c = (uint64_t)params.get_int("c", 0);
    return std::make_unique<TemperMKTrans>(w, mk_type, eta, mu, u, l, b, c);
}

std::vector<ParamSpec> TemperMKTrans::param_specs() {
    return {
        {"w",   "int", true,  false, 0, "",        "", false},
        {"eta", "int", true,  false, 0, "",        "", false},
        {"mu",  "int", true,  false, 0, "",        "", false},
        {"u",   "int", true,  true,  0, "",        "", false},
        {"l",   "int", true,  true,  0, "",        "", false},
        {"b",   "int", false, false, 0, "bitmask", "w", true},
        {"c",   "int", false, false, 0, "bitmask", "w", true},
    };
}
