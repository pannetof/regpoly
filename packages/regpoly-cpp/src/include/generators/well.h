// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <functional>

/**
 * @file well.h
 * @brief Well Equidistributed Long-period Linear generator (Panneton, L'Ecuyer & Matsumoto 2006).
 *
 * `WELLGen` implements the WELL family of F_2-linear word-recurrence
 * generators. The recurrence applies 8 algorithm slots `T0..T7`, each
 * picking from the paper's Table I list of M-class transformations
 * `M0..M6`. The catalog file
 * `docs/library/panneton-lecuyer-matsumoto-2006.yaml` carries the
 * published configurations (WELL512a, WELL19937a, WELL44497a/b, …).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief WELL generator (Panneton, L'Ecuyer & Matsumoto 2006).
 *
 * Structural parameters: `w` (word width), `r` (state size in
 * words), `p` (output-mask bit count), and three offsets `m1, m2,
 * m3` used in the recurrence. The 8 algorithm slots `T0..T7` are
 * described by a structured `matrices: {T0..T7}` map whose entries
 * select an M-class plus its arguments. Registered as `"WELLGen"`
 * (alias `"WELLRNG"`).
 *
 * The flagship WELL19937a uses `w = 32, r = 624, p = 31, m1 = 70,
 * m2 = 179, m3 = 449`, with the matrices map listed in the catalog
 * file.
 *
 * @code{.cpp}
 *   // Construct via the catalog or the YAML config loader — see
 *   // regpoly::library::Catalog::generator() or
 *   // yaml_config::load_seek_config().  Direct factory construction
 *   // requires populating the structured `matrices: {T0..T7}` map:
 *   //   auto gen = create_generator("WELLGen", params, L);
 *   // See docs/library/panneton-lecuyer-matsumoto-2006.yaml for
 *   // ready-to-use parameter sets.
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class WELLGen : public Generator {
public:
    /**
     * @brief One algorithm slot's M-class transformation (paper Table I).
     *
     * `Mi` selects which M-class; the active fields depend on `Mi`:
     *
     *   - M0 = 0       y = 0                            (no args)
     *   - M1 = 1       y = x                            (no args)
     *   - M2(t) = 2    y = x >> t  /  x << -t           (t signed)
     *   - M3(t) = 3    y = x ^ shift(x, t)               (t signed)
     *   - M4(a) = 4    y = (x >> 1) ^ a if LSB(x)        (a 32-bit)
     *   - M5(t, b) = 5 y = x ^ (shift(x, t) & b)         (t signed, b 32-bit)
     *   - M6(q, t, s, a) = 6  rotate-mask-conditional XOR
     *     (q, t, s in [0, w-1], a 32-bit)
     *
     * M6's `d_s` mask and `(1 << t)` test bit are synthesised at
     * runtime from the small-integer fields (matches paper Table I
     * exactly).
     */
    struct MatrixEntry {
        int Mi = 0;        ///< M-class selector (0..6).
        int q  = 0;        ///< M6: rotate amount.
        int t  = 0;        ///< M2/M3/M5: signed shift; M6: bit position 0..w-1.
        int s  = 0;        ///< M6: mask selector 0..w-1.
        uint32_t a = 0;    ///< M4, M6: XOR constant.
        uint32_t b = 0;    ///< M5: mask.
    };

    /**
     * @brief Construct a WELLGen with explicit matrices map.
     *
     * Most callers should go through `create_generator("WELLGen", ...)`
     * rather than this constructor directly.
     *
     * @param w         Word width in bits (only `w = 32` is supported today).
     * @param r         State size in words.
     * @param p         Output mask bit count (lower `p` bits of `state[i+1]`).
     * @param m1        First offset used in the recurrence.
     * @param m2        Second offset used in the recurrence.
     * @param m3        Third offset used in the recurrence.
     * @param matrices  The 8 algorithm slots `T0..T7` in paper order.
     * @param L         Output resolution in bits.
     */
    WELLGen(int w, int r, int p, int m1, int m2, int m3,
              const std::vector<MatrixEntry>& matrices, int L);

    /**
     * @brief Build a WELLGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys w, r, p, m1, m2, m3, matrices.
     * @param L       Output resolution in bits.
     * @return        A constructed WELLGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /**
     * @brief Total cost of this generator's matrices (sum over T0..T7).
     *
     * Costs come from POL / Panneton et al. and are non-monotonic:
     * M0=0, M1=1, M2=2, M3=3, M4=5, M5=4, M6=8.
     */
    int total_cost() const;

    /**
     * @brief Cost of one M-class.
     *
     * Public so the cost-bounded sampler and external tooling (web
     * app, search loop) can compute costs without instantiating a
     * generator.
     *
     * @param Mi  M-class selector in [0, 6].
     * @return    The integer cost.
     */
    static int static_cost_for_Mi(int Mi);

    /**
     * @brief Sample a fresh `matrices` map whose total cost <= `max_cost`.
     *
     * Args per Mi are drawn uniformly from per-class ranges (paper
     * Table I, degenerate values excluded). Algorithm: rejection
     * sampling first (8 i.i.d. uniform Mi draws, accept if sum <= cap,
     * up to 64 attempts), then a greedy-budgeted fallback (shuffle
     * slot order, pick Mi from `{Mi : cost <= remaining_budget}` per
     * slot). The greedy distribution is biased toward expensive Mi
     * early; acceptable because callers run many iterations.
     *
     * @param w         Word width (only `w = 32` is supported today).
     * @param max_cost  Maximum total cost across the 8 slots.
     * @param rng       Pseudo-random generator used for sampling.
     * @return          The sampled `matrices: {T0..T7}` map.
     * @throws std::invalid_argument  If `max_cost <= 0` or `w != 32`.
     */
    static StructMap random_matrices(
        int w, int max_cost, std::mt19937_64& rng);

    /** @brief Family display name — returns the canonical "WELLGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (sizes plus T0..T7 list). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `w * r` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Advance the WELL recurrence by one step (applies all 8 slots). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Current output word truncated to `L()` bits. */
    BitVect get_output() const override;
    /** @brief Pack the canonical state into `out_words` (raw uint64_t form). */
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

private:
    int w_;
    int r_;
    int p_;         // output mask bits
    int m1_, m2_, m3_;
    // 8 algorithm slots T0..T7 in the WELL recurrence; each slot's
    // class is one of M0..M6 from Table I of Panneton et al. (2006).
    std::vector<MatrixEntry> matrices_;
    int i_;         // circular buffer pointer
    int state_bits_; // w * r (full state size for circular buffer)
    uint64_t maskp_;   // UPPER mask: p least significant bits
    uint64_t umaskp_;  // LOWER mask: ~maskp
    static constexpr uint64_t M32 = 0xFFFFFFFFULL;

    static uint64_t ShiftR(uint64_t v, int s);
    static uint64_t apply_matrix(const MatrixEntry& m, uint64_t v);
    uint64_t TMAT(int j, uint64_t val) const;
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);

    static int type_cost(int Mi);
    static std::string type_display(const MatrixEntry& m);
};

}  // namespace regpoly::core
