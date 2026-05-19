// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "transformation.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>

/**
 * @file temper_mk.h
 * @brief Matsumoto–Kurita output tempering (`tempMK*` family).
 *
 * Linear output transformation introduced by Matsumoto & Kurita for
 * the original TGFSR and reused by MT19937 / WELL / MELG. Four
 * `type` variants are dispatched through a single class — the
 * canonical names registered with the transformation factory are
 * `tempMK`, `tempMK2`, `tempMKopt`, and `tempMK2opt`. The `*opt`
 * variants share parameter shape but signal to the search loop that
 * `(b, c)` are bitmask parameters the tempering optimiser may mutate
 * between trials.
 *
 * Bitmask parameters `b` and `c` are flagged `optimizable = true` in
 * `param_specs()`; the structural shifts `eta`, `mu`, `u`, `l` are
 * marked `structural = true` and must be supplied by the caller.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Matsumoto–Kurita-style output tempering (one transformation,
 *        four `type` variants).
 *
 * Implements the per-step linear map applied to the F₂-linear output
 * word: a sequence of left/right shifts of widths `eta` / `mu` /
 * `u` / `l` interleaved with AND-mask folds against the two
 * `w`-bit bitmasks `b` and `c`. The exact composition is selected
 * by the integer `type` field at construction:
 *
 *   - `type = 1` (`tempMK`)    — single-mask MT19937 / WELL form.
 *   - `type = 2` (`tempMK2`)   — double-mask MELG / TinyMT form.
 *   - `*opt` variants          — same shape, marked as targets of
 *                                the bitmask-search loop.
 *
 * Construction goes through the factory (`create_transformation`) and
 * the registry resolves the `type_name` string to a concrete `type`
 * integer. Direct construction sets `type` numerically.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w",   32);
 *   p.set_int("eta", 7);
 *   p.set_int("mu",  15);
 *   p.set_int("u",   11);
 *   p.set_int("l",   18);
 *   p.set_int("b",   0x9D2C5680);
 *   p.set_int("c",   0xEFC60000);
 *   auto t = create_transformation("tempMK2", p);  // MT19937 tempering
 * @endcode
 *
 * @see :py:class:`regpoly.core.transformation.Transformation`  — the Python wrapper.
 * @ingroup core
 */
class TemperMKTrans : public Transformation {
public:
    /**
     * @brief Construct a TemperMK transformation with explicit parameters.
     *
     * Most callers should use `create_transformation("tempMK[2][opt]", ...)`
     * via the factory, which resolves the canonical type name to the
     * integer `type` and validates the parameter shape against
     * `param_specs()`.
     *
     * @param w      Word width in bits.
     * @param type   Variant selector (1 = single-mask, 2 = double-mask).
     * @param eta    First shift amount.
     * @param mu     Second shift amount.
     * @param u      Additional left-shift width.
     * @param l      Additional right-shift width.
     * @param b      First bitmask (`w`-bit, optimizable).
     * @param c      Second bitmask (`w`-bit, optimizable; ignored when type=1).
     */
    TemperMKTrans(int w, int type, int eta, int mu, int u, int l,
                  uint64_t b, uint64_t c);

    /**
     * @brief Factory hook — build from a Params dict + canonical type name.
     * @param type_name  One of `"tempMK"`, `"tempMK2"`, `"tempMKopt"`,
     *                   `"tempMK2opt"`.
     * @param params     Parameter dict with keys `w, eta, mu, u, l, b, c`.
     * @return           Owning pointer to the new TemperMKTrans.
     * @throws std::runtime_error  If a required parameter is missing.
     */
    static std::unique_ptr<Transformation> from_params(
        const std::string& type_name, const Params& params);

    /**
     * @brief Parameter specs declared by this transformation.
     *
     * Returns the seven-entry vector with `w / eta / mu / u / l`
     * structural and `b / c` flagged `optimizable = true` (so the
     * tempering optimiser knows to mutate them).
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Canonical type name (`tempMK` / `tempMK2` / `tempMK*opt`). */
    std::string name() const override;

    /** @brief Human-readable parameter dump for the regpoly-cli `show` mode. */
    std::string display_str() const override;

    /** @brief Apply the in-place linear map to `state` (the output word). */
    void apply(BitVect& state) const override;

    /** @brief Deep clone — independent of `*this`. */
    std::unique_ptr<Transformation> copy() const override;

    /**
     * @brief Mutate bitmask parameters from `params` (optimiser hook).
     *
     * The tempering optimiser calls this between trials to write new
     * `b` / `c` values without rebuilding the transformation. Structural
     * fields (`eta`, `mu`, `u`, `l`) are not touched.
     *
     * @param params  Subset of param dict containing the keys to update.
     */
    void update(const Params& params) override;

    /** @brief Variant selector (1 = single-mask, 2 = double-mask). */
    int type() const { return type_; }
    /** @brief First shift amount. */
    int eta() const { return eta_; }
    /** @brief Second shift amount. */
    int mu() const { return mu_; }
    /** @brief Additional left-shift width. */
    int u() const { return u_; }
    /** @brief Additional right-shift width. */
    int l() const { return l_; }
    /** @brief First bitmask (optimizable). */
    uint64_t b() const { return b_; }
    /** @brief Second bitmask (optimizable). */
    uint64_t c() const { return c_; }

private:
    int type_;
    int eta_, mu_;
    int u_, l_;
    uint64_t b_, c_;
};

}  // namespace regpoly::core
