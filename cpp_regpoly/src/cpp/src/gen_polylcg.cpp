#include "gen_polylcg.h"
#include <sstream>
#include <iomanip>
#include <vector>

PolyLCG::PolyLCG(int k, const BitVect& poly, int L)
    : Generateur(k, L), poly_(poly.copy()) {}

std::string PolyLCG::name() const { return "Polynomial LCG"; }

std::string PolyLCG::display_str() const {
    // Collect exponents where poly_ has bit set
    // exponents = [k - i - 1 for i in range(k) if poly.get_bit(i) == 1]
    std::vector<int> exponents;
    for (int i = 0; i < k_; i++) {
        if (poly_.get_bit(i) == 1)
            exponents.push_back(k_ - i - 1);
    }

    // line1: " {k} " + " ".join(str(e) for e in exponents) + " "
    std::ostringstream line1;
    line1 << " " << k_ << " ";
    for (size_t i = 0; i < exponents.size(); i++) {
        if (i > 0) line1 << " ";
        line1 << exponents[i];
    }
    line1 << " ";

    // line2: "hexadecimal notation:\n {poly_val:0{hex_digits}x} "
    // Compute poly_val: the integer value of poly_ (nbits = k_)
    // poly_ is a BitVect of k_ bits. We need to reconstruct the integer.
    int hex_digits = (k_ + 3) / 4;

    // Build the integer value from bits
    // bit 0 is MSB, bit k-1 is LSB
    // val = sum of (get_bit(i) << (k - 1 - i)) for i in 0..k-1
    // But that's the same as reading the top k bits as an integer
    // We can use get_word if k <= 64, otherwise we need to be careful
    std::ostringstream hex_ss;
    if (k_ <= 64) {
        uint64_t poly_val = poly_.get_word(0, k_);
        hex_ss << std::hex << std::setfill('0') << std::setw(hex_digits) << poly_val;
    } else {
        // For k > 64, build hex string from groups of 4 bits
        for (int d = 0; d < hex_digits; d++) {
            int nibble = 0;
            for (int b = 0; b < 4; b++) {
                int bit_idx = k_ - hex_digits * 4 + d * 4 + b;
                if (bit_idx >= 0 && bit_idx < k_) {
                    nibble |= (poly_.get_bit(bit_idx) << (3 - b));
                }
            }
            hex_ss << std::hex << nibble;
        }
    }

    std::ostringstream line2;
    line2 << "hexadecimal notation:\n " << hex_ss.str() << " ";

    return line1.str() + "\n" + line2.str();
}

void PolyLCG::init(const BitVect& init_bv) {
    state_ = init_bv.copy();
}

void PolyLCG::next() {
    int xor_flag = state_.get_bit(0);
    state_.lshift(1);
    if (xor_flag)
        state_.xor_with(poly_);
    state_.and_mask(k_);
}

std::unique_ptr<Generateur> PolyLCG::copy() const {
    auto p = std::make_unique<PolyLCG>(k_, poly_, L_);
    p->state_ = state_.copy();
    return p;
}
