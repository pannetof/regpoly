// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#include "digital_net.h"
#include "bitvect.h"

#include <sstream>
#include <stdexcept>

namespace regpoly::core {

namespace {

// F_2 matrix-vector product: result[r] = parity of (row r) AND vec.
// Rows are stored as `BitVect`s of width `m`; `vec` is also width `m`.
// Output is a width-`m` `BitVect`.
BitVect matvec_f2(const std::vector<BitVect>& rows, const BitVect& vec) {
    const int m = static_cast<int>(rows.size());
    BitVect out(m);
    for (int r = 0; r < m; ++r) {
        // XOR rows[r] AND vec, then popcount % 2.
        const int nwords = rows[r].nwords();
        uint64_t acc = 0;
        const uint64_t* a = rows[r].data();
        const uint64_t* b = vec.data();
        for (int w = 0; w < nwords; ++w) {
            acc ^= a[w] & b[w];
        }
        // Parity of acc.
        acc ^= acc >> 32;
        acc ^= acc >> 16;
        acc ^= acc >> 8;
        acc ^= acc >> 4;
        acc ^= acc >> 2;
        acc ^= acc >> 1;
        out.set_bit(r, static_cast<int>(acc & 1ULL));
    }
    return out;
}

}  // namespace

DigitalNet::DigitalNet(int m, int s_max)
    : Generator(m, m), j_(0), s_max_(s_max) {
    if (m <= 0) {
        throw std::invalid_argument(
            "DigitalNet: m must be positive (got " + std::to_string(m) + ")");
    }
    if (m > 64) {
        // v1 cap; the t-value kernel's branch-and-bound is calibrated for
        // m <= 64. Beyond that, the Niederreiter–Pirsic dual method
        // (currently a stub) is the appropriate kernel.
        throw std::invalid_argument(
            "DigitalNet: m exceeds v1 cap of 64 (got " + std::to_string(m) + ")");
    }
    if (s_max < m) {
        std::ostringstream oss;
        oss << "DigitalNet: s_max (" << s_max << ") must be >= m (" << m << "). "
            << "The equidistribution kernel passes indice_max = kg = m to "
               "GaussMatrix::prepare, which calls next() m-1 times; a "
               "DigitalNet's next() throws past s_max-1.";
        throw std::invalid_argument(oss.str());
    }
}

void DigitalNet::init(const BitVect& init_bv) {
    if (init_bv.nbits() != k_) {
        throw std::invalid_argument(
            "DigitalNet::init: seed width " + std::to_string(init_bv.nbits())
            + " does not match m = " + std::to_string(k_));
    }
    // Zero seed is legal for a digital net; index 0 is the first point.
    state_ = init_bv.copy();
    j_ = 1;
}

void DigitalNet::next() {
    if (j_ >= s_max_) {
        std::ostringstream oss;
        oss << "DigitalNet::next: coordinate " << j_ << " is at the s_max ("
            << s_max_ << ") bound; cannot advance further. "
            << "Caller likely passed indice_max > s_max to GaussMatrix::prepare.";
        throw DimensionExceededError(oss.str());
    }
    ++j_;
}

BitVect DigitalNet::get_output() const {
    if (j_ < 1 || j_ > s_max_) {
        std::ostringstream oss;
        oss << "DigitalNet::get_output: coordinate j_ = " << j_
            << " is out of range [1, " << s_max_ << "]. "
            << "Call init(seed) before consuming output.";
        throw DimensionExceededError(oss.str());
    }
    const auto& rows = generating_matrix(j_);
    return matvec_f2(rows, state_);
}

std::optional<std::string>
DigitalNet::compute_default_test_method(const std::string& test_type) const {
    if (test_type == "equidistribution") return std::string("matricial");
    if (test_type == "tvalue")           return std::string("schmid");
    return std::nullopt;
}

}  // namespace regpoly::core
