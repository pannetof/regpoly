// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "equidistribution_method.h"
#include "generator.h"

#include <memory>
#include <string>
#include <vector>

// Polymorphic Test for the seek-search loop.
//
// Replaces the SeekTestKind switch that previously lived in
// seek_search.cpp. Each Test:
//   - knows its type_name() ("equidistribution" / "collision_free" / "tuplets"),
//   - runs against a Generator and writes the relevant subset of
//     SeekIterResult fields,
//   - reports whether the iteration passes its short-circuit predicate.
//
// CollisionFree never short-circuits (passed() always returns true).
// Equidistribution short-circuits when the kernel says verified and
// se > mse. Tuplets short-circuits when the kernel reports
// not-ok. The seek_search loop orchestrates this uniformly.
//
// Adding a new test type means writing one Test subclass and
// registering it via TestRegistry::reg("baz", ...) in factory.cpp.
// No existing dispatch site needs editing.

// Forward declaration to avoid circular include.
struct SeekIterResult;

class Test {
public:
    virtual ~Test() = default;

    // Stable identifier used by the YAML parser, the Python wrapper,
    // and the registry. Matches the canonical test_type names
    // ("equidistribution", "collision_free", "tuplets").
    virtual std::string type_name() const = 0;

    // Run the test against `gen` (kg, L are the active combination's
    // numbers). Updates the relevant subset of `iter`. Tests that
    // require state from a previous test (collision_free reads
    // iter.me_test_L) consult `iter` directly.
    virtual void run(SeekIterResult& iter,
                     const Generator& gen,
                     int kg, int L) const = 0;

    // True iff the iteration must continue past this test. Equivalent
    // to "this test passed". The seek loop short-circuits on the first
    // false. Tests that never short-circuit (collision_free) always
    // return true regardless of result.
    virtual bool passed(const SeekIterResult& iter) const = 0;
};


// ── Concrete tests ────────────────────────────────────────────────────

// Equidistribution test. Owns its method (Matricial, Lattice, …) — the
// method abstraction handles per-algorithm differences; this class
// handles the fields-and-short-circuit logic.
class EquidistributionTestRunner : public Test {
public:
    EquidistributionTestRunner(
        std::unique_ptr<EquidistributionMethod> method,
        int eq_L_max_test,
        std::vector<int> delta,
        int mse);

    std::string type_name() const override { return "equidistribution"; }
    void run(SeekIterResult& iter,
             const Generator& gen,
             int kg, int L) const override;
    bool passed(const SeekIterResult& iter) const override;

    int mse() const { return mse_; }
    std::string method_name() const { return method_->name(); }

private:
    std::unique_ptr<EquidistributionMethod> method_;
    int eq_L_max_test_;
    std::vector<int> delta_;
    int mse_;
};


// Collision-free test. Never short-circuits.
class CollisionFreeTestRunner : public Test {
public:
    std::string type_name() const override { return "collision_free"; }
    void run(SeekIterResult& iter,
             const Generator& gen,
             int kg, int L) const override;
    bool passed(const SeekIterResult&) const override { return true; }
};


// Tuplets test. Short-circuits on tup_is_ok == false.
class TupletsTestRunner : public Test {
public:
    TupletsTestRunner(int d,
                      std::vector<int> h,
                      double threshold,
                      int testtype);

    std::string type_name() const override { return "tuplets"; }
    void run(SeekIterResult& iter,
             const Generator& gen,
             int kg, int L) const override;
    bool passed(const SeekIterResult& iter) const override;

private:
    int d_;
    std::vector<int> h_;
    double threshold_;
    int testtype_;
};
