// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "params.h"
#include <memory>
#include <optional>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
namespace regpoly::core {

class Generator {
public:
    Generator(int k, int L) : state_(k), k_(k), L_(L) {}
    virtual ~Generator() = default;

    // Pure virtual interface — makes Generator abstract
    virtual std::string name() const = 0;
    virtual std::string display_str() const = 0;
    virtual void init(const BitVect& init_bv) = 0;
    virtual void next() = 0;
    virtual std::unique_ptr<Generator> copy() const = 0;

    // Public accessors
    int k() const { return k_; }
    int L() const { return L_; }
    const BitVect& state() const { return state_; }

    // Overridable behaviour
    virtual BitVect get_output() const;
    virtual void get_transition_state(uint64_t* out_words, int out_nwords) const;

    // Number of consecutive L-bit output words emitted from the SAME
    // F2-linear state before the next state-advance (next() call that
    // actually applies the recurrence).  Default = 1: every next()
    // advances the state.
    //
    // Generators that batch their internal recurrence (e.g. SFMTGen runs
    // generate_all() once per N 128-bit words = 4N L=32-bit words)
    // override this so equidistribution analyses can correctly partition
    // output rows by phase within a state-cycle.  See §6.4 of
    // docs/C.md.  This is consulted only via the abstract interface;
    // me_notprimitive.cpp does not branch on generator type.
    virtual int output_phases() const { return 1; }

    // Number of L-bit lanes packed into one SIMD super-word.  Default 1
    // means "no SIMD packing" — the generator emits one independent
    // L-bit value per next() call.  SIMD generators (SFMTGen, dSFMT, MTGPGen)
    // override this: e.g. SFMTGen-32 returns 4 (a 128-bit internal word
    // split into four 32-bit lanes).  Consumed by the simd_notprimitive
    // equidistribution method (test_me_notprimitive_simd in me_notprimitive_simd.h),
    // which assembles `simd_lane_count()` consecutive next() outputs
    // into a 128-bit super-word and runs the SIMD-aware PIS reduction
    // from Saito–Matsumoto 2008.  No code path in cpp_regpoly branches
    // on concrete generator type — this is the only hint the SIMD
    // module needs.
    virtual int simd_lane_count() const { return 1; }

    // Period exponent: log2(period + 1) for full-period generators.
    // Default = k() (full-period assumption).  Generators whose true
    // F2-linear state is larger than the period exponent (SFMTGen, dSFMT,
    // MTGPGen — non-primitive characteristic polynomial) override this to
    // return the period exponent (e.g. SFMTGen-19937 returns 19937 even
    // though k() = 19968).  Consumed by SIMD/notprimitive methods to
    // compute the d(v) = floor(period_exponent/v) - k(v) gap, matching
    // the convention used in published equidistribution tables.
    virtual int period_exponent() const { return k_; }

    // Set the F2-linear state directly without triggering any f-application
    // side effect (no generate_all, no buf_idx reset).  Default: copy s into
    // state_.  SIMD generators with internal scratch buffers (SFMTGen's work_)
    // override this to keep all internal state mirrors in sync.  Required
    // by the SIMD-PIS reduction, which XORs basis-vector states between
    // PIS iterations without wanting to advance the recurrence.
    virtual void set_raw_state(const BitVect& s) {
        state_ = s.copy();
    }

    // ── SIMD-aware primitives for METHOD_SIMD_NOTPRIMITIVE ─────────────
    //
    // These primitives let me_notprimitive_simd.cpp run a SIMD-aware PIS reduction
    // without knowing the concrete generator type.  Defaults collapse
    // to existing single-lane behaviour; SIMD generators (SFMTGen) override
    // to do per-word f-applications and lane-rotated super-word reads
    // matching MTToolBox's AlgorithmSIMDEquidistribution semantics.

    // Apply ONE fundamental f-application to the state.  For non-SIMD
    // generators this is just next().  For SFMTGen-style generators that
    // batch internal updates, override to update exactly ONE 128-bit
    // word per call (one inner-loop iteration of generate_all).
    virtual void simd_advance_one_word() { next(); }

    // Read the current "super-word" (128-bit BitVect) for SIMD-PIS.
    // start_mode ∈ {0..lane_count-1} rotates which lanes come from the
    // just-updated word vs the previous word.  Default (non-SIMD)
    // ignores start_mode and zero-pads get_output() to 128 bits.
    virtual BitVect simd_read_super_word(int /*start_mode*/) const {
        BitVect out(128);
        BitVect cur = get_output();
        int n = std::min(L_, cur.nbits());
        for (int b = 0; b < n; b++)
            if (cur.get_bit(b)) out.set_bit(b, 1);
        return out;
    }

    // XOR another generator's state into this state, in-place, with NO
    // f-application side effect.  Required by SIMD-PIS basis reductions.
    // Default: state_ ^= other.state_.  SFMTGen override also XORs work_.
    virtual void simd_add_state(const Generator& other) {
        BitVect s = state_.copy();
        s.xor_with(other.state_);
        set_raw_state(s);
    }

    // Reset the per-word advance index (for SFMTGen's pword_idx_) so a
    // sequence of simd_advance_one_word() calls starts from word 0.
    // Default: no-op (non-SIMD generators have no such index).
    virtual void simd_reset_word_index() {}

    // Number of simd_advance_one_word() calls that compose one full
    // f-application (= one generate_all batch).  For SFMTGen this equals
    // k/128 = N, but for dSFMT k = 128*(N+1) (lung is part of the
    // F2-linear state) while one full batch is still N do_recursion
    // calls — dSFMT overrides to return N.  Default = k/128 (assumes
    // every 128-bit slot of state corresponds to one per-word advance,
    // which matches SFMTGen and any future SIMD generator without a
    // separate carry register).
    virtual int simd_full_step_words() const { return k_ / 128; }

    // Algorithms (use the virtual interface above)
    // Default: Berlekamp-Massey. Overridden by TGFSRGen, MT, PolyLCGGen, MELGGen.
    virtual BitVect char_poly() const;

    // Returns true if the characteristic polynomial is primitive over GF(2),
    // meaning the generator has full period 2^k - 1.
    bool is_full_period() const;

    // Returns the recommended method for the given test_type, or
    // std::nullopt if no default is defined for that test_type.
    // Memoised — the underlying computation can call is_full_period()
    // (Berlekamp-Massey + factor 2^k - 1), which is expensive for
    // large k. The cache is per-instance and NOT thread-safe; callers
    // that share a single Generator across threads must serialise
    // access externally. (Today's callers build fresh per-request
    // Combinations, so no instance is shared.)
    std::optional<std::string> default_test_method(const std::string& test_type) const;

protected:
    // Subclasses with a known answer (SIMD families, CombinedGenerator)
    // override this; the base implementation handles "equidistribution"
    // via the is_full_period() / k() runtime rule. Callers must go
    // through default_test_method() to benefit from memoisation.
    virtual std::optional<std::string> compute_default_test_method(const std::string& test_type) const;

public:

    // Compute the K×K transition matrix.
    // Returns K row BitVects, where row[i] has bit j set if A[i][j] = 1.
    std::vector<BitVect> transition_matrix() const;

    // Component decomposition for kernels that operate on a list of
    // primitive generators plus a per-component tempering chain.
    // The default treats `*this` as a single primitive with no
    // tempering — exactly the right answer for every concrete subclass
    // EXCEPT CombinedGenerator, which overrides to expose its components.
    //
    // Eliminates the dynamic_cast<CombinedGenerator*> formerly scattered
    // across single_gen_adapters.cpp, equidistribution_runner.cpp, and
    // tuplets_runner.cpp.
    //
    // Returned pointers alias the Generator's internal components and
    // are valid for the lifetime of `*this`. Forward declaration so
    // generator.h does not need to include transformation.h.
    virtual std::vector<Generator*> components() const;
    virtual std::vector<std::vector<class Transformation*>> tempering_chains() const;

protected:
    BitVect state_;
    int k_;
    int L_;

private:
    mutable std::unordered_map<std::string, std::optional<std::string>>
        default_test_method_cache_;
};

}  // namespace regpoly::core
