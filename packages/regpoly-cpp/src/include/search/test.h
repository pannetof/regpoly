// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "equidistribution_method.h"
#include "generator.h"

#include <memory>
#include <string>
#include <vector>

/**
 * @file test.h
 * @brief Polymorphic search-test predicates evaluated per combo.
 * @ingroup core
 *
 * Polymorphic `Test` for the seek-search loop. Replaces the
 * `SeekTestKind` switch that previously lived in `seek_search.cpp`.
 * Each `Test`:
 *  - knows its `type_name()` (`"equidistribution"` /
 *    `"collision_free"` / `"tuplets"`),
 *  - runs against a `Generator` and writes the relevant subset of
 *    `SeekIterResult` fields,
 *  - reports whether the iteration passes its short-circuit predicate.
 *
 * `CollisionFree` never short-circuits (`passed()` always returns
 * true). `Equidistribution` short-circuits when the kernel says
 * verified and `se > mse`. `Tuplets` short-circuits when the kernel
 * reports not-ok. The seek-search loop orchestrates this uniformly.
 *
 * Adding a new test type means writing one `Test` subclass and
 * registering it via `TestRegistry::reg("baz", ...)` in `factory.cpp`.
 * No existing dispatch site needs editing.
 */

// Forward declaration to avoid circular include.
namespace regpoly::core {

struct SeekIterResult;

/**
 * @brief Abstract base for a search predicate evaluated per combo.
 *
 * Subclasses implement one of the canonical test types
 * (`equidistribution`, `collision_free`, `tuplets`); they run against
 * a `Generator`, populate the relevant slots of `SeekIterResult`, and
 * report whether the surrounding seek loop should keep evaluating
 * later tests for this combo.
 *
 * @ingroup core
 */
class Test {
public:
    virtual ~Test() = default;

    /**
     * @brief Stable test identifier (`"equidistribution"`, etc.).
     *
     * Matches the canonical test-type names used by the YAML parser,
     * the Python wrapper, and the registry.
     *
     * @return  Test type name.
     */
    virtual std::string type_name() const = 0;

    /**
     * @brief Run the test against `gen` and update `iter` in place.
     *
     * Tests that require state from a previous test
     * (`collision_free` reads `iter.me_test_L`) consult `iter`
     * directly rather than receiving a separate context object.
     *
     * @param iter  Iteration result being assembled (mutated in place).
     * @param gen   Generator to evaluate.
     * @param kg    Combined state size of the active combination.
     * @param L     Output word width of the active combination.
     */
    virtual void run(SeekIterResult& iter,
                     const Generator& gen,
                     int kg, int L) const = 0;

    /**
     * @brief Short-circuit predicate consulted after `run`.
     *
     * The seek loop stops evaluating further tests for the current
     * combo on the first `false`. Tests that never short-circuit
     * (`collision_free`) always return true regardless of result.
     *
     * @param iter  Iteration result populated by `run`.
     * @return      True iff the iteration may continue past this test.
     */
    virtual bool passed(const SeekIterResult& iter) const = 0;
};


// ── Concrete tests ────────────────────────────────────────────────────

/**
 * @brief Equidistribution search predicate.
 *
 * Owns its `EquidistributionMethod` (Matricial, Lattice, Harase,
 * NotPrimitive, SimdNotPrimitive, Nothing). The method abstraction
 * handles per-algorithm differences; this class handles the field
 * population and short-circuit logic (fails when the kernel is
 * verified and `se > mse`).
 *
 * @see :py:class:`regpoly.analyses.equidistribution_test.EquidistributionTest`
 *
 * @ingroup core
 */
class EquidistributionTestRunner : public Test {
public:
    /**
     * @brief Construct an equidistribution test predicate.
     *
     * @param method        Owning equidistribution method (from `MethodRegistry`).
     * @param eq_L_max_test Maximum resolution to test.
     * @param delta         Per-resolution gap budget (size `eq_L_max_test + 1`).
     * @param mse           Upper bound on the cumulative SE.
     */
    EquidistributionTestRunner(
        std::unique_ptr<EquidistributionMethod> method,
        int eq_L_max_test,
        std::vector<int> delta,
        int mse);

    /** @brief Test type name — `"equidistribution"`. */
    std::string type_name() const override { return "equidistribution"; }
    void run(SeekIterResult& iter,
             const Generator& gen,
             int kg, int L) const override;
    bool passed(const SeekIterResult& iter) const override;

    /** @brief Cumulative-SE cap configured at construction. */
    int mse() const { return mse_; }
    /** @brief Method name reported by the underlying `EquidistributionMethod`. */
    std::string method_name() const { return method_->name(); }

private:
    std::unique_ptr<EquidistributionMethod> method_;
    int eq_L_max_test_;
    std::vector<int> delta_;
    int mse_;
};


/**
 * @brief Collision-free search predicate (never short-circuits).
 *
 * Runs the collision-free rank deficit walk on the active generator
 * and writes the result into `SeekIterResult::cf_*`. `passed()`
 * always returns true — collision-free is informational and never
 * stops the test stack.
 *
 * @see :py:class:`regpoly.analyses.collision_free_test.CollisionFreeTest`
 *
 * @ingroup core
 */
class CollisionFreeTestRunner : public Test {
public:
    /** @brief Test type name — `"collision_free"`. */
    std::string type_name() const override { return "collision_free"; }
    void run(SeekIterResult& iter,
             const Generator& gen,
             int kg, int L) const override;
    /** @brief Never short-circuits; always returns true. */
    bool passed(const SeekIterResult&) const override { return true; }
};


/**
 * @brief Tuplets search predicate.
 *
 * Runs the tuplets uniformity test on the active generator and
 * short-circuits when the kernel reports `tup_is_ok == false`.
 *
 * @see :py:class:`regpoly.analyses.tuplets_test.TupletsTest`
 *
 * @ingroup core
 */
class TupletsTestRunner : public Test {
public:
    /**
     * @brief Construct a tuplets test predicate.
     *
     * @param d          Tuplets dimension.
     * @param h          Tuplets shape parameters (1-indexed).
     * @param threshold  Rejection threshold.
     * @param testtype   Test type (0 = SUM, 1 = MAX).
     */
    TupletsTestRunner(int d,
                      std::vector<int> h,
                      double threshold,
                      int testtype);

    /** @brief Test type name — `"tuplets"`. */
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

}  // namespace regpoly::core
