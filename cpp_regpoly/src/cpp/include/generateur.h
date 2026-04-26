#pragma once
#include "bitvect.h"
#include "params.h"
#include <memory>
#include <string>
#include <sstream>
#include <vector>
class Generateur {
public:
    Generateur(int k, int L) : state_(k), k_(k), L_(L) {}
    virtual ~Generateur() = default;

    // Pure virtual interface — makes Generateur abstract
    virtual std::string name() const = 0;
    virtual std::string display_str() const = 0;
    virtual void init(const BitVect& init_bv) = 0;
    virtual void next() = 0;
    virtual std::unique_ptr<Generateur> copy() const = 0;

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
    // Generators that batch their internal recurrence (e.g. SFMT runs
    // generate_all() once per N 128-bit words = 4N L=32-bit words)
    // override this so equidistribution analyses can correctly partition
    // output rows by phase within a state-cycle.  See §6.4 of
    // docs/C.md.  This is consulted only via the abstract interface;
    // notprimitive_de.cpp does not branch on generator type.
    virtual int output_phases() const { return 1; }

    // Number of L-bit lanes packed into one SIMD super-word.  Default 1
    // means "no SIMD packing" — the generator emits one independent
    // L-bit value per next() call.  SIMD generators (SFMT, dSFMT, MTGP)
    // override this: e.g. SFMT-32 returns 4 (a 128-bit internal word
    // split into four 32-bit lanes).  Consumed by the simd_notprimitive
    // equidistribution method (test_me_simd_notprimitive in simd_de.h),
    // which assembles `simd_lane_count()` consecutive next() outputs
    // into a 128-bit super-word and runs the SIMD-aware PIS reduction
    // from Saito–Matsumoto 2008.  No code path in cpp_regpoly branches
    // on concrete generator type — this is the only hint the SIMD
    // module needs.
    virtual int simd_lane_count() const { return 1; }

    // Period exponent: log2(period + 1) for full-period generators.
    // Default = k() (full-period assumption).  Generators whose true
    // F2-linear state is larger than the period exponent (SFMT, dSFMT,
    // MTGP — non-primitive characteristic polynomial) override this to
    // return the period exponent (e.g. SFMT-19937 returns 19937 even
    // though k() = 19968).  Consumed by SIMD/notprimitive methods to
    // compute the d(v) = floor(period_exponent/v) - k(v) gap, matching
    // the convention used in published equidistribution tables.
    virtual int period_exponent() const { return k_; }

    // Set the F2-linear state directly without triggering any f-application
    // side effect (no generate_all, no buf_idx reset).  Default: copy s into
    // state_.  SIMD generators with internal scratch buffers (SFMT's work_)
    // override this to keep all internal state mirrors in sync.  Required
    // by the SIMD-PIS reduction, which XORs basis-vector states between
    // PIS iterations without wanting to advance the recurrence.
    virtual void set_raw_state(const BitVect& s) {
        state_ = s.copy();
    }

    // ── SIMD-aware primitives for METHOD_SIMD_NOTPRIMITIVE ─────────────
    //
    // These primitives let simd_de.cpp run a SIMD-aware PIS reduction
    // without knowing the concrete generator type.  Defaults collapse
    // to existing single-lane behaviour; SIMD generators (SFMT) override
    // to do per-word f-applications and lane-rotated super-word reads
    // matching MTToolBox's AlgorithmSIMDEquidistribution semantics.

    // Apply ONE fundamental f-application to the state.  For non-SIMD
    // generators this is just next().  For SFMT-style generators that
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
    // Default: state_ ^= other.state_.  SFMT override also XORs work_.
    virtual void simd_add_state(const Generateur& other) {
        BitVect s = state_.copy();
        s.xor_with(other.state_);
        set_raw_state(s);
    }

    // Reset the per-word advance index (for SFMT's pword_idx_) so a
    // sequence of simd_advance_one_word() calls starts from word 0.
    // Default: no-op (non-SIMD generators have no such index).
    virtual void simd_reset_word_index() {}

    // Algorithms (use the virtual interface above)
    // Default: Berlekamp-Massey. Overridden by TGFSR, MT, PolyLCG, MELG.
    virtual BitVect char_poly() const;

    // Returns true if the characteristic polynomial is primitive over GF(2),
    // meaning the generator has full period 2^k - 1.
    bool is_full_period() const;

    // Compute the K×K transition matrix.
    // Returns K row BitVects, where row[i] has bit j set if A[i][j] = 1.
    std::vector<BitVect> transition_matrix() const;

protected:
    BitVect state_;
    int k_;
    int L_;
};

