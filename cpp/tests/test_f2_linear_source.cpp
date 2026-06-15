// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

// Phase 2 Step 1 — source-side scaffolding. Validates that:
//   1. F2LinearSource is instantiable through a minimal concrete subclass
//      (which means the additive design of Phase 1 + Step 1 is well-formed).
//   2. get_output() returns RAW bits during the Phase 2 migration window
//      (matches today's Recurrence::get_output() semantics).
//   3. tempered_output() returns chain ∘ raw_output().
//   4. clone_source() returns a typed unique_ptr<F2LinearSource> whose
//      tempering chain is a deep copy of the donor's.
//
// At this phase no production class inherits from F2LinearSource. The
// TestSource fixture is local to this test file.

#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>

#include "bitvect.h"
#include "f2_linear_source.h"
#include "tempering_chain.h"
#include "temper_mk.h"

using namespace regpoly::core;


namespace {

// Build the standard MT19937 tempMK2 step (b=0x9D2C5680, c=0xEFC60000).
std::unique_ptr<Transformation> make_mt_temper() {
    return std::make_unique<TemperMKTrans>(
        /*w=*/32, /*type=*/2,
        /*eta=*/7, /*mu=*/15, /*u=*/11, /*l=*/18,
        /*b=*/0x9D2C5680ULL, /*c=*/0xEFC60000ULL);
}

// Minimal concrete F2LinearSource: state is one 32-bit word; raw_output()
// returns it verbatim. next() and init() are unused by the tests below.
class TestSource : public F2LinearSource {
public:
    TestSource(uint32_t initial)
        : F2LinearSource(/*k=*/32, /*L=*/32), state_(32) {
        state_.set_word(0, 32, initial);
    }

    std::string name() const override { return "TestSource"; }
    std::string display_str() const override { return "TestSource(stub)"; }
    void init(const BitVect& seed) override { state_ = seed.copy(); }
    void next() override { /* no-op for the scaffold tests */ }

    BitVect raw_output() const override { return state_.copy(); }

    std::unique_ptr<ITestable> copy() const override {
        return std::make_unique<TestSource>(*this);
    }

    std::optional<std::string>
    default_test_method(const std::string&) const override {
        return std::nullopt;
    }

private:
    BitVect state_;
};

}  // namespace


TEST(F2LinearSource, GetOutputReturnsTempered) {
    // get_output() returns chain ∘ raw_output(). Compared against an
    // independent oracle (apply the same chain step directly).
    TestSource src(0xDEADBEEF);

    // Empty chain: get_output() == raw_output().
    EXPECT_EQ(src.get_output().get_word(0, 32),
              src.raw_output().get_word(0, 32));

    // Non-empty chain: tempered diverges from raw, matches oracle.
    TemperingChain chain;
    chain.add(make_mt_temper());
    src.set_tempering(std::move(chain));

    uint64_t raw      = src.raw_output().get_word(0, 32);
    uint64_t tempered = src.get_output().get_word(0, 32);
    EXPECT_NE(raw, tempered);

    BitVect direct = src.raw_output();
    auto oracle = make_mt_temper();
    oracle->apply(direct);
    EXPECT_EQ(tempered, direct.get_word(0, 32));
}


TEST(F2LinearSource, CloneSourceIsDeep) {
    TestSource src(0xCAFEBABE);
    TemperingChain chain;
    chain.add(make_mt_temper());
    src.set_tempering(std::move(chain));

    std::unique_ptr<F2LinearSource> clone = src.clone_source();

    // Same outputs (raw + tempered) — same state, same chain.
    EXPECT_EQ(clone->raw_output().get_word(0, 32),
              src.raw_output().get_word(0, 32));
    EXPECT_EQ(clone->get_output().get_word(0, 32),
              src.get_output().get_word(0, 32));

    // Independent storage: chain steps are different C++ objects.
    ASSERT_EQ(clone->tempering().size(), 1);
    ASSERT_EQ(src.tempering().size(), 1);
    EXPECT_NE(&clone->tempering().at(0), &src.tempering().at(0));
}


TEST(F2LinearSource, CopyCtorDeepClonesChain) {
    // Direct copy ctor (used by leaf-class `copy()` overrides via
    // `std::make_unique<MTGen>(*this)`).
    TestSource original(0xFEEDFACE);
    TemperingChain chain;
    chain.add(make_mt_temper());
    original.set_tempering(std::move(chain));

    TestSource copied = original;

    EXPECT_EQ(copied.tempering().size(), 1);
    EXPECT_NE(&copied.tempering().at(0), &original.tempering().at(0));

    EXPECT_EQ(copied.get_output().get_word(0, 32),
              original.get_output().get_word(0, 32));
}
