// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.
//
// Phase 1 smoke test for the DigitalNet abstract class.
//
// IdentityNet (C_j = I for all j) lets us drive every DigitalNet
// contract — init/next/get_output, the s_max>=m invariant, the
// bounded next() throw — without dragging in a concrete
// implementation. It also exercises the architectural claim that
// the unmodified equidistribution kernel runs on a DigitalNet.

#include <gtest/gtest.h>

#include "digital_net.h"
#include "bitvect.h"
#include "gauss.h"
#include "transformation.h"
#include "equidistribution_runner.h"

#include <climits>
#include <memory>
#include <vector>

using regpoly::core::BitVect;
using regpoly::core::DigitalNet;
using regpoly::core::GaussMatrix;
using regpoly::core::Generator;
using regpoly::core::run_matricial_equidistribution;

namespace {

// IdentityNet — every generating matrix is the m×m identity.
// Outputs at coordinate j (1..s_max) are always equal to the index.
class IdentityNet : public DigitalNet {
public:
    IdentityNet(int m, int s_max)
        : DigitalNet(m, s_max) {
        // Build a single shared "identity" matrix of m rows × m cols.
        // Row r has bit r set (the r-th standard basis row).
        std::vector<BitVect> rows;
        rows.reserve(m);
        for (int r = 0; r < m; ++r) {
            BitVect row(m);
            row.set_bit(r, 1);
            rows.push_back(std::move(row));
        }
        identity_rows_ = std::make_shared<const std::vector<BitVect>>(std::move(rows));
    }

    std::string name() const override { return "IdentityNet"; }
    std::string display_str() const override {
        return "IdentityNet(m=" + std::to_string(m()) +
               ", s_max=" + std::to_string(s_max()) + ")";
    }
    std::unique_ptr<Generator> copy() const override {
        // Default copy-ctor clones state_, j_, s_max_, and the
        // shared_ptr (cheap reference bump on identity_rows_).
        return std::unique_ptr<Generator>(new IdentityNet(*this));
    }

protected:
    const std::vector<BitVect>& generating_matrix(int /*j*/) const override {
        return *identity_rows_;
    }

private:
    std::shared_ptr<const std::vector<BitVect>> identity_rows_;
};

}  // namespace

// ── Construction invariants ────────────────────────────────────────────

TEST(DigitalNet, ConstructorRejectsTooSmallSMax) {
    EXPECT_THROW(IdentityNet(/*m=*/4, /*s_max=*/3), std::invalid_argument);
}

TEST(DigitalNet, ConstructorAcceptsSMaxEqualToM) {
    EXPECT_NO_THROW(IdentityNet(/*m=*/4, /*s_max=*/4));
}

TEST(DigitalNet, ConstructorAcceptsSMaxLargerThanM) {
    EXPECT_NO_THROW(IdentityNet(/*m=*/4, /*s_max=*/10));
}

TEST(DigitalNet, ConstructorRejectsNonPositiveM) {
    EXPECT_THROW(IdentityNet(/*m=*/0, /*s_max=*/0), std::invalid_argument);
    EXPECT_THROW(IdentityNet(/*m=*/-1, /*s_max=*/5), std::invalid_argument);
}

TEST(DigitalNet, ConstructorRejectsExcessiveM) {
    EXPECT_THROW(IdentityNet(/*m=*/65, /*s_max=*/65), std::invalid_argument);
}

// ── j-counter contract ─────────────────────────────────────────────────

TEST(DigitalNet, InitSetsJ1AndStateMatchesSeed) {
    IdentityNet net(4, 4);
    BitVect seed(4);
    seed.set_bit(1, 1);  // index = 0100
    net.init(seed);
    EXPECT_EQ(net.current_j(), 1);
    EXPECT_EQ(net.state().get_bit(1), 1);
    EXPECT_EQ(net.state().get_bit(0), 0);
}

TEST(DigitalNet, InitAcceptsZeroSeed) {
    IdentityNet net(4, 4);
    BitVect seed(4);  // all zeros
    EXPECT_NO_THROW(net.init(seed));
    EXPECT_EQ(net.current_j(), 1);
}

TEST(DigitalNet, InitRejectsWidthMismatch) {
    IdentityNet net(4, 4);
    BitVect bad(5);
    EXPECT_THROW(net.init(bad), std::invalid_argument);
}

TEST(DigitalNet, NextIncrementsCoordinate) {
    IdentityNet net(4, 5);
    BitVect seed(4);
    seed.set_bit(0, 1);
    net.init(seed);
    EXPECT_EQ(net.current_j(), 1);
    net.next();
    EXPECT_EQ(net.current_j(), 2);
    net.next();
    EXPECT_EQ(net.current_j(), 3);
}

TEST(DigitalNet, NextThrowsPastSMax) {
    IdentityNet net(/*m=*/4, /*s_max=*/4);
    BitVect seed(4);
    seed.set_bit(0, 1);
    net.init(seed);
    EXPECT_EQ(net.current_j(), 1);
    net.next();          // j=2
    net.next();          // j=3
    net.next();          // j=4
    EXPECT_THROW(net.next(), DigitalNet::DimensionExceededError);
}

TEST(DigitalNet, GetOutputThrowsWithoutInit) {
    IdentityNet net(4, 4);
    EXPECT_THROW(net.get_output(), DigitalNet::DimensionExceededError);
}

// ── Identity-matrix output semantics ──────────────────────────────────

TEST(DigitalNet, IdentityGetOutputReturnsIndex) {
    IdentityNet net(4, 4);
    BitVect seed(4);
    seed.set_bit(0, 1);   // 1000
    seed.set_bit(2, 1);   // 1010
    net.init(seed);
    BitVect out = net.get_output();
    EXPECT_EQ(out.get_bit(0), 1);
    EXPECT_EQ(out.get_bit(1), 0);
    EXPECT_EQ(out.get_bit(2), 1);
    EXPECT_EQ(out.get_bit(3), 0);
    // Same output at every coordinate (identity).
    net.next();
    BitVect out2 = net.get_output();
    EXPECT_EQ(out2.get_bit(0), 1);
    EXPECT_EQ(out2.get_bit(2), 1);
}

// ── copy() produces an independent clone ──────────────────────────────

TEST(DigitalNet, CopyClonesStateAndCounter) {
    IdentityNet net(4, 4);
    BitVect seed(4);
    seed.set_bit(1, 1);
    net.init(seed);
    net.next();
    EXPECT_EQ(net.current_j(), 2);

    auto clone = net.copy();
    auto* dn = dynamic_cast<DigitalNet*>(clone.get());
    ASSERT_NE(dn, nullptr);
    EXPECT_EQ(dn->current_j(), 2);
    EXPECT_EQ(dn->state().get_bit(1), 1);

    // Advancing the clone leaves the original alone.
    dn->next();
    EXPECT_EQ(dn->current_j(), 3);
    EXPECT_EQ(net.current_j(), 2);
}

// ── Default test methods ──────────────────────────────────────────────

TEST(DigitalNet, DefaultMethodForEquidistributionIsMatricial) {
    IdentityNet net(4, 4);
    auto m = net.default_test_method("equidistribution");
    ASSERT_TRUE(m.has_value());
    EXPECT_EQ(*m, "matricial");
}

TEST(DigitalNet, DefaultMethodForTValueIsSchmid) {
    IdentityNet net(4, 4);
    auto m = net.default_test_method("tvalue");
    ASSERT_TRUE(m.has_value());
    EXPECT_EQ(*m, "schmid");
}

TEST(DigitalNet, DefaultMethodForUnknownTestIsNullopt) {
    IdentityNet net(4, 4);
    auto m = net.default_test_method("nonsense_test");
    EXPECT_FALSE(m.has_value());
}

// ── The architectural claim: GaussMatrix::prepare runs unchanged ──────
//
// For an IdentityNet, every C_j = I, so the M matrix built by
// `GaussMatrix::prepare(..., kg, indice_max=kg, L=m)` should have row
// `gl` = the gl-th column of I, repeated `kg` times. That column is
// the basis vector e_gl, which has bit `gl` set.

TEST(DigitalNet, PrepareBuildsRepeatedIdentityRows) {
    const int m = 4;
    IdentityNet net(m, m);
    std::vector<Generator*> gens = {&net};
    std::vector<int> gen_k = {m};
    std::vector<std::vector<regpoly::core::Transformation*>> trans = {{}};

    GaussMatrix mat = GaussMatrix::prepare(gens, gen_k, trans,
                                           /*kg=*/m,
                                           /*indice_max=*/m,
                                           /*L=*/m);
    EXPECT_EQ(mat.nrows(), m);
    EXPECT_EQ(mat.ncols(), m * m);

    // Row gl, columns [j*m, (j+1)*m): only column j*m + gl is set
    // (the gl-th bit of C_j · e_gl = gl-th bit of e_gl).
    for (int gl = 0; gl < m; ++gl) {
        for (int j = 0; j < m; ++j) {
            for (int b = 0; b < m; ++b) {
                bool expected = (b == gl);
                bool actual = mat.bit_test(gl, j * m + b);
                EXPECT_EQ(actual, expected)
                    << "row=" << gl << " block=" << j << " bit=" << b;
            }
        }
    }
}

// ── The architectural claim part 2: equidistribution runs unchanged ───
//
// For IdentityNet (every C_j = I), the equidistribution dimension k(l)
// at resolution l is 1 — only the first row (C_1's first l rows) gives
// independent vectors; subsequent C_j blocks are duplicates. So the
// matricial kernel should report ecart[l] = floor(kg/l) - 1 for every
// resolution l = 1..m.

TEST(DigitalNet, MatricialEquidistributionRunsAndReportsDegenerate) {
    const int m = 4;
    IdentityNet net(m, m);
    std::vector<int> delta(m + 1, INT_MAX);
    auto res = run_matricial_equidistribution(
        net, /*kg=*/m, /*L=*/m, /*Lmax=*/m, delta, /*mse=*/INT_MAX);

    // At resolution l=1, kg/l = 4 blocks; the rank picks up only the
    // first block. ecart[1] = 4 - 1 = 3.
    EXPECT_EQ(res.ecart[1], 3);
    // At resolution l=2, kg/l = 2 blocks; ecart[2] = 2 - 1 = 1.
    EXPECT_EQ(res.ecart[2], 1);
    // At resolution l=4, kg/l = 1 block; ecart[4] = 1 - 1 = 0.
    EXPECT_EQ(res.ecart[4], 0);
    EXPECT_TRUE(res.verified);
}
