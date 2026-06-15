// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "f2_linear_source.h"
#include "i_testable.h"
#include "params.h"
#include <cassert>
#include <memory>
#include <optional>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

/**
 * @file generator.h
 * @brief Abstract base class for every F_2-linear PRNG family in REGPOLY.
 *
 * `regpoly::core::Recurrence` (legacy alias: `Recurrence`) is the single
 * interface every concrete recurrence-driven PRNG family (MT, WELL,
 * LFSR, SFMT, dSFMT, MarsaXorshift, Tausworthe, Cellular Automata, …)
 * inherits from. Subclasses implement the pure-virtual contract
 * (`name`, `display_str`, `init`, `next`, `copy`) and may override the
 * SIMD-aware optional methods documented below if they batch their
 * internal recurrence (SFMT, dSFMT, MTGP).
 *
 * `Recurrence` is a sibling of `DigitalNet` under `F2LinearSource` —
 * the rename from `Generator` to `Recurrence` (2026-05) encodes the
 * architectural distinction between "recurrence-driven PRNG" and
 * "digital net" in the type name.
 *
 * Higher-level kernels — equidistribution analyses, tuplets runner,
 * combined-generator iteration — never branch on concrete subclass
 * type; everything goes through this interface plus the
 * `components()` virtual.
 *
 * @see :py:class:`regpoly.core.generator.Recurrence`
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Abstract base class for every recurrence-driven F_2-linear PRNG in REGPOLY.
 *
 * State model: every recurrence carries a `BitVect` of width `k_` that
 * represents the current F_2-linear state. `next()` advances the state
 * one fundamental step (applies the F_2-linear recurrence). The output
 * word is read via `get_output()` and has width `L_` bits.
 * `copy()` returns an independent deep clone — the search loop relies
 * on this for every candidate generator it considers.
 *
 * **Optional SIMD overrides.** Recurrences that batch their internal
 * step (SFMT runs `generate_all()` once per N 128-bit words, dSFMT
 * similarly) override the `simd_*` family of virtuals so the
 * SIMD-aware notprimitive equidistribution kernel can run a PIS
 * reduction without knowing the concrete subclass. The defaults
 * collapse to single-lane behaviour.
 *
 * **Default test method.** `default_test_method()` returns the
 * recommended equidistribution method for a given test type. Its
 * computation can be expensive (Berlekamp–Massey + factorisation of
 * `2^k - 1`) and is therefore memoised per-instance; the cache is
 * **not** thread-safe — callers sharing a single `Recurrence` across
 * threads must serialise. Current callers build fresh per-request
 * combinations, so no instance is shared.
 *
 * @code{.cpp}
 *   // A Recurrence is normally received from the factory or catalog,
 *   // not constructed directly. The exact factory call depends on
 *   // the caller's setup; e.g.
 *   //   auto g = make_via_factory("MTGen", params, L);
 *   using namespace regpoly::core;
 *   std::unique_ptr<Recurrence> g = some_factory_call();
 *   BitVect seed(g->k());
 *   seed.set_bit(0, 1);
 *   g->init(seed);
 *   for (int i = 0; i < 4; ++i) {
 *       g->next();
 *       // Consume g->get_output().
 *   }
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Recurrence`
 *
 * @ingroup core
 */
class Recurrence : public F2LinearSource {
public:
    /**
     * @brief Construct a recurrence with state width `k` and output width `L`.
     *
     * Allocates the internal `state_` `BitVect` of width `k` (all
     * zeros). Subclasses are responsible for completing initialisation
     * in their own constructors and / or via `init()`.
     *
     * `k` and `L` are forwarded to the `F2LinearSource` base, which
     * stores them as `k_` / `L_` and exposes `k()` / `L()` accessors.
     *
     * @param k  State width in bits.
     * @param L  Output word width in bits.
     */
    Recurrence(int k, int L) : F2LinearSource(k, L), state_(k) {}

    // ── Pure virtual interface (makes Recurrence abstract) ────────────

    /** @brief Canonical family name (e.g. `"MTGen"`, `"WELLGen"`). */
    std::string name() const override = 0;

    /**
     * @brief Human-readable parametrised string for diagnostics.
     *
     * Typically includes the family name and the structural parameters
     * (e.g. `"MTGen(k=19937, L=32)"`).
     *
     * @return  Display string.
     */
    std::string display_str() const override = 0;

    /**
     * @brief Initialise the F_2-linear state from `init_bv`.
     *
     * Implementations must accept any non-zero seed of width `k()`
     * and leave the generator ready to emit its first output via
     * the usual `next()` → `get_output()` cycle.
     *
     * @param init_bv  Seed of width `k()` (must be non-zero).
     */
    void init(const BitVect& init_bv) override = 0;

    /**
     * @brief Advance the F_2-linear state by one fundamental step.
     *
     * Concrete subclasses apply their family's recurrence here.
     * For SIMD families that batch internal updates, "one fundamental
     * step" matches what `output_phases()` describes.
     */
    void next() override = 0;

    /**
     * @brief Return an independent deep clone of this recurrence.
     *
     * The clone must produce identical outputs given identical seeds.
     * Used pervasively by the search loop, the combined-generator
     * builder, and the equidistribution kernels. Returns
     * `unique_ptr<ITestable>` to satisfy the polymorphic `ITestable`
     * interface; callers that need typed `Recurrence*` access should
     * use `clone_recurrence()` instead.
     *
     * @return  A `unique_ptr<ITestable>` to a freshly heap-allocated clone.
     */
    std::unique_ptr<ITestable> copy() const override = 0;

    /**
     * @brief Deep clone, typed as `unique_ptr<Recurrence>`.
     *
     * Non-virtual forwarder that wraps the polymorphic `copy()` and
     * recovers the `Recurrence` static type via `static_cast`. Every
     * `Recurrence` subclass's `copy()` returns a `Recurrence`-derived
     * pointer by construction; debug builds verify via `dynamic_cast`.
     *
     * Used by kernels (`me_notprimitive`, `me_harase`, `temper_optimizer`,
     * `combo_enumerator`) that need to clone a recurrence while
     * preserving the `Recurrence` API surface (`state()`, `char_poly()`,
     * `simd_*`).
     */
    std::unique_ptr<Recurrence> clone_recurrence() const {
        auto p = copy();
        assert(dynamic_cast<Recurrence*>(p.get()) != nullptr
               && "Recurrence subclass copy() must return a "
                  "Recurrence-derived pointer");
        return std::unique_ptr<Recurrence>(
            static_cast<Recurrence*>(p.release()));
    }

    // ── Public accessors ────────────────────────────────────────────

    // k() and L() are inherited from F2LinearSource.

    /** @brief Read-only view of the current F_2-linear state. */
    const BitVect& state() const { return state_; }

    // ── Overridable behaviour ───────────────────────────────────────

    /**
     * @brief Return the un-tempered output word (overrides
     *        `F2LinearSource::raw_output`).
     *
     * Default implementation returns the `L` most significant bits of
     * `state_`. Concrete PRNG families that emit something other than
     * the top `L` bits override this. `F2LinearSource::get_output()`
     * wraps this (returns raw during the Phase 2 migration window;
     * will flip to tempered at the end of Phase 2).
     *
     * @return  Un-tempered output as a `BitVect` of width `L()`.
     */
    BitVect raw_output() const override;

    /**
     * @brief Fill `out_words` with the transition-state words used by analyses.
     *
     * The exact serialisation depends on the family; concrete
     * subclasses document their packing. Used by analyses that need
     * the full F_2-linear state laid out as 64-bit words.
     *
     * @param out_words   Output buffer.
     * @param out_nwords  Buffer capacity in 64-bit words.
     */
    virtual void get_transition_state(uint64_t* out_words, int out_nwords) const;

    /**
     * @brief Number of consecutive L-bit output words emitted per state-advance.
     *
     * Number of consecutive `L`-bit output words that come out of the
     * SAME F_2-linear state before the next state-advance (the next
     * `next()` call that actually applies the recurrence). Default = 1:
     * every `next()` advances the state.
     *
     * Generators that batch their internal recurrence (e.g. `SFMTGen`
     * runs `generate_all()` once per N 128-bit words = 4·N L=32-bit
     * words) override this so equidistribution analyses can correctly
     * partition output rows by phase within a state-cycle. See §6.4 of
     * `docs/C.md`. This is consulted only via the abstract interface;
     * `me_notprimitive.cpp` does not branch on generator type.
     *
     * @return  Number of output words per fundamental f-application.
     */
    virtual int output_phases() const { return 1; }

    /**
     * @brief Number of L-bit lanes packed into one SIMD super-word.
     *
     * Default = 1, meaning "no SIMD packing" — the generator emits one
     * independent `L`-bit value per `next()` call. SIMD generators
     * (`SFMTGen`, `dSFMT`, `MTGPGen`) override this: e.g. `SFMTGen-32`
     * returns 4 (a 128-bit internal word split into four 32-bit lanes).
     *
     * Consumed by the `simd_notprimitive` equidistribution method
     * (`test_me_notprimitive_simd` in `me_notprimitive_simd.h`), which
     * assembles `simd_lane_count()` consecutive `next()` outputs into
     * a 128-bit super-word and runs the SIMD-aware PIS reduction from
     * Saito–Matsumoto 2008. No code path in `cpp_regpoly` branches on
     * concrete generator type — this is the only hint the SIMD module
     * needs.
     *
     * @return  Lane count per 128-bit super-word.
     */
    virtual int simd_lane_count() const { return 1; }

    /**
     * @brief Period exponent of the generator.
     *
     * `log2(period + 1)` for full-period generators. Default = `k()`
     * (full-period assumption). Generators whose true F_2-linear state
     * is larger than the period exponent (SFMT, dSFMT, MTGP —
     * non-primitive characteristic polynomial) override this to return
     * the period exponent (e.g. `SFMTGen-19937` returns 19937 even
     * though `k() = 19968`).
     *
     * Consumed by SIMD / notprimitive methods to compute the
     * `d(v) = floor(period_exponent/v) - k(v)` gap, matching the
     * convention used in published equidistribution tables.
     *
     * @return  Period exponent.
     */
    virtual int period_exponent() const { return k_; }

    /**
     * @brief Set the F_2-linear state directly with no f-application side effect.
     *
     * No `generate_all`, no `buf_idx` reset. Default: copy `s` into
     * `state_`. SIMD generators with internal scratch buffers
     * (`SFMTGen`'s `work_`) override this to keep all internal state
     * mirrors in sync. Required by the SIMD-PIS reduction, which XORs
     * basis-vector states between PIS iterations without wanting to
     * advance the recurrence.
     *
     * @param s  New state (width must equal `k()`).
     */
    virtual void set_raw_state(const BitVect& s) {
        state_ = s.copy();
    }

    // ── SIMD-aware primitives for METHOD_SIMD_NOTPRIMITIVE ───────────
    //
    // These primitives let `me_notprimitive_simd.cpp` run a SIMD-aware
    // PIS reduction without knowing the concrete generator type.
    // Defaults collapse to existing single-lane behaviour; SIMD
    // generators (`SFMTGen`) override to do per-word f-applications
    // and lane-rotated super-word reads matching MTToolBox's
    // `AlgorithmSIMDEquidistribution` semantics.

    /**
     * @brief Apply ONE fundamental f-application to the state.
     *
     * For non-SIMD generators this is just `next()`. For SFMT-style
     * generators that batch internal updates, override to update
     * exactly ONE 128-bit word per call (one inner-loop iteration of
     * `generate_all`).
     */
    virtual void simd_advance_one_word() { next(); }

    /**
     * @brief Read the current 128-bit "super-word" for SIMD-PIS.
     *
     * `start_mode` is in `{0, ..., lane_count - 1}` and rotates which
     * lanes come from the just-updated word vs the previous word.
     * Default (non-SIMD) ignores `start_mode` and zero-pads
     * `get_output()` up to 128 bits.
     *
     * @param start_mode  Lane-rotation index (ignored by the default).
     * @return            128-bit super-word as a `BitVect`.
     */
    virtual BitVect simd_read_super_word(int /*start_mode*/) const {
        // Q2 (Phase 2): SIMD-PIS reads raw (un-tempered) bits — the
        // Saito-Matsumoto 2008 reduction operates on the state-space
        // recurrence, which is inherently pre-tempering. Routing through
        // `raw_output()` keeps the contract correct after `get_output()`
        // flips to return tempered bits at the end of Phase 2.
        BitVect out(128);
        BitVect cur = raw_output();
        int n = std::min(L_, cur.nbits());
        for (int b = 0; b < n; b++)
            if (cur.get_bit(b)) out.set_bit(b, 1);
        return out;
    }

    /**
     * @brief XOR another recurrence's state into this state in place.
     *
     * No f-application side effect. Required by SIMD-PIS basis
     * reductions. Default: `state_ ^= other.state_`. SFMT override
     * also XORs the internal `work_` scratch buffer.
     *
     * @param other  Recurrence whose state is XORed into `*this`.
     */
    virtual void simd_add_state(const Recurrence& other) {
        BitVect s = state_.copy();
        s.xor_with(other.state_);
        set_raw_state(s);
    }

    /**
     * @brief Reset the per-word advance index.
     *
     * For SFMT's `pword_idx_`, so a sequence of
     * `simd_advance_one_word()` calls starts from word 0. Default:
     * no-op (non-SIMD generators have no such index).
     */
    virtual void simd_reset_word_index() {}

    /**
     * @brief Number of `simd_advance_one_word()` calls per full f-application.
     *
     * For SFMT this equals `k / 128 = N`, but for dSFMT
     * `k = 128 * (N + 1)` (the lung is part of the F_2-linear state)
     * while one full batch is still N `do_recursion` calls — dSFMT
     * overrides to return `N`. Default = `k / 128` (assumes every
     * 128-bit slot of state corresponds to one per-word advance,
     * which matches SFMT and any future SIMD generator without a
     * separate carry register).
     *
     * @return  Number of per-word advances composing one full step.
     */
    virtual int simd_full_step_words() const { return k_ / 128; }

    // ── Algorithms (use the virtual interface above) ────────────────

    /**
     * @brief Compute the characteristic polynomial of the state-transition.
     *
     * Default: Berlekamp–Massey driven through `next()` and
     * `get_output()`. Overridden by `TGFSRGen`, the MT family,
     * `PolyLCGGen`, and `MELGGen` with closed-form expressions.
     *
     * @return  Characteristic polynomial as a `BitVect`.
     */
    virtual BitVect char_poly() const;

    /**
     * @brief True iff the characteristic polynomial is primitive over GF(2).
     *
     * Equivalent to saying the generator has full period `2^k - 1`.
     * Internally calls `char_poly()` and factors `2^k - 1`, so the
     * call is expensive for large `k`.
     *
     * @return  True iff full-period.
     */
    bool is_full_period() const;

    /**
     * @brief Recommended equidistribution method for `test_type`.
     *
     * Memoised — the underlying computation can call
     * `is_full_period()` (Berlekamp–Massey + factor `2^k - 1`), which
     * is expensive for large `k`. The cache is per-instance and
     * **not** thread-safe; callers that share a single `Recurrence`
     * across threads must serialise access externally. Today's callers
     * build fresh per-request combinations, so no instance is shared.
     *
     * @param test_type  Test family name (e.g. `"equidistribution"`).
     * @return           Method name, or `std::nullopt` if no default
     *                   is defined for that `test_type`.
     */
    std::optional<std::string> default_test_method(const std::string& test_type) const override;

protected:
    /**
     * @brief Compute the default test method (override point for subclasses).
     *
     * Subclasses with a known answer (SIMD families, `CombinedF2LinearSource`)
     * override this; the base implementation handles
     * `"equidistribution"` via the `is_full_period()` / `k()` runtime
     * rule. Callers must go through `default_test_method()` to
     * benefit from memoisation.
     *
     * @param test_type  Test family name.
     * @return           Method name, or `std::nullopt` if undefined.
     */
    virtual std::optional<std::string> compute_default_test_method(const std::string& test_type) const;

public:

    /**
     * @brief Compute the K×K transition matrix.
     *
     * @return  K row `BitVect`s, where row `i` has bit `j` set iff
     *          `A[i][j] = 1`.
     */
    std::vector<BitVect> transition_matrix() const;

    /**
     * @brief Component recurrences (decomposition for combined-aware kernels).
     *
     * The default treats `*this` as a single primitive — exactly the
     * right answer for every concrete subclass EXCEPT
     * `CombinedF2LinearSource`, which overrides to expose its components.
     *
     * Eliminates the `dynamic_cast<CombinedF2LinearSource*>` formerly
     * scattered across `single_gen_adapters.cpp`,
     * `equidistribution_runner.cpp`, and `tuplets_runner.cpp`.
     *
     * Returned pointers alias the recurrence's internal components
     * and are valid for the lifetime of `*this`.
     *
     * @return  Vector of component pointers (`{this}` by default).
     */
    virtual std::vector<Recurrence*> components() const;

protected:
    // k_ and L_ are inherited from F2LinearSource (protected).
    BitVect state_;     ///< Current F_2-linear state (width `k_`).

private:
    mutable std::unordered_map<std::string, std::optional<std::string>>
        default_test_method_cache_;
};

}  // namespace regpoly::core
