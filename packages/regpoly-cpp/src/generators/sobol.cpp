// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#include "sobol.h"
#include "bitvect.h"

#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace regpoly::core {

namespace {

// Embedded direction-number table for coordinates j = 2..8 (the first 7
// non-identity coordinates). Values from Joe & Kuo (2003) Table 1;
// reproduced verbatim across QMCPy, Numerical Recipes §7.7,
// MathWorks' sobolset, and Wikipedia's "Sobol" article — these seven
// initial rows are the most-cross-checked values in the Sobol literature.
//
// Layout: kEmbeddedTable[k] is the data for coordinate j = k + 2.
// "a" is the binary encoding of the interior coefficients of the
// primitive polynomial (excluding the leading and trailing 1s); the
// Antonov–Saleev recurrence interprets bit (deg - 1 - i) of `a` as
// the i-th interior coefficient a_i.
const std::vector<SobolDirNumbers> kEmbeddedTable = {
    /* j=2  */ { /*deg=*/1, /*a=*/0, /*m_init=*/{1}                },
    /* j=3  */ { /*deg=*/2, /*a=*/1, /*m_init=*/{1, 3}             },
    /* j=4  */ { /*deg=*/3, /*a=*/1, /*m_init=*/{1, 3, 1}          },
    /* j=5  */ { /*deg=*/3, /*a=*/2, /*m_init=*/{1, 1, 1}          },
    /* j=6  */ { /*deg=*/4, /*a=*/1, /*m_init=*/{1, 1, 3, 3}       },
    /* j=7  */ { /*deg=*/4, /*a=*/4, /*m_init=*/{1, 3, 5, 13}      },
    /* j=8  */ { /*deg=*/5, /*a=*/2, /*m_init=*/{1, 1, 5, 5, 17}   },
};

// Build the m×m identity as a row-major vector<BitVect>: row r has bit r set.
std::vector<BitVect> build_identity_rows(int m) {
    std::vector<BitVect> rows;
    rows.reserve(m);
    for (int r = 0; r < m; ++r) {
        BitVect row(m);
        row.set_bit(r, 1);
        rows.push_back(std::move(row));
    }
    return rows;
}

// Convert m direction-number integers v[1..m] (where v[k] occupies
// the k high-order bits of an m-bit field) into a row-major
// generating matrix: column k-1 of C_j is the m-bit MSB-first
// encoding of v[k].
//
// Result is m row-`BitVect`s of width m. Row r, column k-1 = bit r of v[k].
std::vector<BitVect> matrix_from_v(const std::vector<uint64_t>& v, int m) {
    std::vector<BitVect> rows(m, BitVect(m));
    for (int k = 1; k <= m; ++k) {
        uint64_t vk = v[k];
        // bit r of vk (counted from MSB of the m-bit field) = (vk >> (m-1-r)) & 1.
        for (int r = 0; r < m; ++r) {
            int bit = static_cast<int>((vk >> (m - 1 - r)) & 1ULL);
            if (bit) rows[r].set_bit(k - 1, 1);
        }
    }
    return rows;
}

// Compute direction numbers m_{j,1..M} for coordinate j (j ≥ 2) using
// the Antonov–Saleev recurrence on Joe-Kuo's encoding of the
// primitive polynomial.
//
// Returns m_full[1..M] (index 0 unused) with m_full[k] = m_{j,k}.
// m_full[k] is an integer in [0, 2^k); for k > s, computed from the
// recurrence in §sobol.h.
std::vector<uint64_t> compute_m_sequence(const SobolDirNumbers& row, int M) {
    const int s = row.deg;
    std::vector<uint64_t> m_full(M + 1, 0);
    for (int k = 1; k <= s && k <= M; ++k) {
        m_full[k] = static_cast<uint64_t>(row.m_init[k - 1]);
    }
    for (int k = s + 1; k <= M; ++k) {
        // m_{j,k} = 2·a_1·m_{j,k-1} XOR ... XOR 2^{s-1}·a_{s-1}·m_{j,k-s+1}
        //          XOR 2^s·m_{j,k-s} XOR m_{j,k-s}
        uint64_t acc = 0;
        for (int i = 1; i <= s - 1; ++i) {
            // bit (s-1-i) of `a` indicates whether a_i is set.
            uint64_t a_i = (row.a >> (s - 1 - i)) & 1ULL;
            acc ^= a_i * (m_full[k - i] << i);
        }
        // Last two terms (whose "coefficient" is always 1).
        acc ^= m_full[k - s] << s;
        acc ^= m_full[k - s];
        m_full[k] = acc;
    }
    return m_full;
}

// For coordinate j ≥ 2, build the m×m generating matrix as row-`BitVect`s.
std::vector<BitVect> build_C_j(const SobolDirNumbers& row, int m) {
    auto m_full = compute_m_sequence(row, m);
    std::vector<uint64_t> v(m + 1, 0);
    // v[k] = m_{j,k} << (m - k)  — left-align the k-bit direction
    // number into an m-bit MSB field.
    for (int k = 1; k <= m; ++k) {
        v[k] = m_full[k] << (m - k);
    }
    return matrix_from_v(v, m);
}

}  // namespace

int SobolNet::embedded_table_size() {
    return static_cast<int>(kEmbeddedTable.size());
}

SobolNet::SobolNet(int m, int s_max, const std::vector<SobolDirNumbers>& table)
    : DigitalNet(m, s_max) {
    if (static_cast<int>(table.size()) + 1 < s_max) {
        std::ostringstream oss;
        oss << "SobolNet: table covers coordinates 2.." << (table.size() + 1)
            << " but s_max = " << s_max
            << " requires data through coordinate " << s_max << ".";
        throw std::invalid_argument(oss.str());
    }
    build_matrices(m, s_max, table);
}

SobolNet::SobolNet(int m, int s_max)
    : SobolNet(m, s_max, kEmbeddedTable) {}

void SobolNet::build_matrices(int m, int s_max,
                              const std::vector<SobolDirNumbers>& table) {
    auto mats = std::make_shared<std::vector<std::vector<BitVect>>>();
    mats->reserve(s_max);
    // Coordinate j=1: identity.
    mats->push_back(build_identity_rows(m));
    // Coordinates j=2..s_max: Antonov–Saleev recurrence.
    for (int j = 2; j <= s_max; ++j) {
        const SobolDirNumbers& row = table[j - 2];
        if (row.deg <= 0 || static_cast<int>(row.m_init.size()) < row.deg) {
            std::ostringstream oss;
            oss << "SobolNet: malformed direction-number row for j=" << j
                << " (deg=" << row.deg << ", m_init.size()=" << row.m_init.size()
                << ").";
            throw std::invalid_argument(oss.str());
        }
        mats->push_back(build_C_j(row, m));
    }
    matrices_ = std::shared_ptr<const std::vector<std::vector<BitVect>>>(
        std::move(mats));
}

const std::vector<BitVect>& SobolNet::generating_matrix(int j) const {
    return (*matrices_)[j - 1];
}

std::string SobolNet::display_str() const {
    std::ostringstream oss;
    oss << "SobolNet(m=" << m() << ", s_max=" << s_max() << ")";
    return oss.str();
}

std::unique_ptr<Generator> SobolNet::copy() const {
    return std::unique_ptr<Generator>(new SobolNet(*this));
}

}  // namespace regpoly::core
