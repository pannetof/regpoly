#include "me_harase.h"
#include "me_helpers.h"  // polychar_comb, MeLatResult
#include <climits>
#include <algorithm>
#include <cstring>

// ══════════════════════════════════════════════════════════════════════════
// Harase PIS method — state-space equidistribution computation.
//
// Implements the algorithm from AlgorithmEquidistribution in MTToolBox
// (Saito/Matsumoto/Harase 2011).  Uses generator copies for all state
// operations — Generator is NOT modified.
//
// The basis has v+1 vectors:
//   basis[0..v-1] = "standard basis" (output one bit at position i, then 0)
//   basis[v]      = the actual generator
//
// Each vector is a set of J generator copies (one per component).
// Shift by z = one call to next() on each component.
// XOR = extract transition states, XOR word arrays, init back.
// ══════════════════════════════════════════════════════════════════════════

namespace {

// ── Lattice vector: J generator copies + output word + degree ───────────

struct LatticeVec {
    // J component generators + transformation chains (not owned)
    std::vector<std::unique_ptr<Generator>> gens;
    std::vector<std::vector<Transformation*>> trans;  // raw ptrs
    int J;

    uint64_t next;  // top v bits of combined (tempered) output
    int count;      // number of advances = polynomial degree
    bool zero;      // vector has been absorbed

    LatticeVec() : J(0), next(0), count(0), zero(false) {}

    // Deep copy of this vector (copies all generator states)
    LatticeVec deep_copy() const {
        LatticeVec v;
        v.J = J;
        v.trans = trans;
        for (auto& g : gens) v.gens.push_back(g->copy());
        v.next = next;
        v.count = count;
        v.zero = zero;
        return v;
    }

    // Clone from template generators (deep copy of all components)
    static LatticeVec from_generators(
        const std::vector<Generator*>& tpl,
        const std::vector<std::vector<Transformation*>>& tr)
    {
        LatticeVec v;
        v.J = (int)tpl.size();
        v.trans = tr;
        for (auto* g : tpl) v.gens.push_back(g->copy());
        v.next = 0;
        v.count = 0;
        v.zero = false;
        return v;
    }

    // Set all component states to zero
    void set_zero() {
        for (auto& g : gens) {
            BitVect z(g->k());
            g->init(z);
        }
    }

    // Advance all components by one step, read top v bits of combined output
    void next_state(int v) {
        if (zero) return;
        int zero_count = 0;
        int max_zero = 0;
        for (auto& g : gens) max_zero += g->k();
        max_zero *= 2;

        do {
            // Advance all components
            for (auto& g : gens) g->next();
            count++;

            // Get combined (tempered) output
            next = get_combined_output(v);

            if (next != 0) return;

            zero_count++;
            if (zero_count > max_zero) {
                zero = true;
                return;
            }
        } while (true);
    }

    // XOR this vector's state with src's state
    void add(const LatticeVec& src) {
        // XOR each component's transition state
        for (int j = 0; j < J; j++) {
            int ki = gens[j]->k();
            int nw = (ki + 63) / 64;
            std::vector<uint64_t> state_a(nw, 0), state_b(nw, 0);

            gens[j]->get_transition_state(state_a.data(), nw);
            src.gens[j]->get_transition_state(state_b.data(), nw);

            // XOR
            for (int w = 0; w < nw; w++)
                state_a[w] ^= state_b[w];

            // Init back from XOR'd state
            BitVect bv(ki);
            for (int w = 0; w < std::min(nw, bv.nwords()); w++)
                bv.data()[w] = state_a[w];
            gens[j]->init(bv);
        }
        next ^= src.next;
    }

private:
    // Get top v bits of combined (tempered) output as a uint64_t
    uint64_t get_combined_output(int v) const {
        uint64_t result = 0;
        for (int j = 0; j < J; j++) {
            BitVect out = gens[j]->get_output();
            // Apply transformations
            if (!trans[j].empty()) {
                for (auto* t : trans[j]) t->apply(out);
            }
            // Extract top v bits (bit 0..v-1 in MSB-first = top v bits of first word)
            uint64_t word = out.top_word();
            result ^= word;
        }
        // Mask to top v bits
        if (v < 64)
            result &= ~0ULL << (64 - v);
        return result;
    }
};

// ── Find rightmost set bit position (0 = MSB) ──────────────────────────
// Returns -1 if x == 0.

static inline int calc_1pos(uint64_t x) {
    if (x == 0) return -1;
    return 63 - __builtin_ctzll(x);  // position from MSB
}

}  // namespace

// ══════════════════════════════════════════════════════════════════════════
// test_me_harase — PIS method for equidistribution
// ══════════════════════════════════════════════════════════════════════════

MeLatResult test_me_harase(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    int RES = std::min(maxL, kg);
    int w = L;  // output word width (typically 64)

    MeLatResult result;
    result.ecart.resize(maxL + 1, -1);
    result.se = 0;

    // Allocate basis: v+1 vectors, where v starts at RES
    int v = RES;
    int size = v + 1;
    std::vector<LatticeVec> basis(size);

    // basis[0..v-1]: standard basis vectors.
    // Each outputs a single 1-bit at position i, then zeros.
    // Implemented as: zero-state generator with next = (1 << (63-i)).
    for (int i = 0; i < v; i++) {
        basis[i] = LatticeVec::from_generators(gens, trans);
        basis[i].set_zero();
        basis[i].count = 0;
        basis[i].zero = false;
        basis[i].next = 1ULL << (63 - i);
    }

    // basis[v]: the actual generator, initialized with canonical state e_0
    basis[v] = LatticeVec::from_generators(gens, trans);
    {
        // Init each component with canonical e_0
        int pos = 0;
        for (int j = 0; j < basis[v].J; j++) {
            int ki = basis[v].gens[j]->k();
            BitVect e0(ki);
            e0.set_bit(0, 1);
            basis[v].gens[j]->init(e0);
            pos += ki;
        }
    }
    basis[v].next_state(v);

    // ── Main loop: compute k(v) for v = RES down to 1 ──────────────────

    int maxl = maxL;
    for (int bit_len = RES; bit_len >= 1; bit_len--) {
        // Get equidistribution dimension k(bit_len)
        int pivot_index = calc_1pos(basis[bit_len].next);

        while (!basis[bit_len].zero) {
            if (pivot_index == -1) break;
            if (pivot_index >= bit_len) break;

            // Swap to keep counts balanced (reduce the one with larger count)
            if (basis[bit_len].count > basis[pivot_index].count) {
                std::swap(basis[bit_len], basis[pivot_index]);
            }

            // Gauss elimination: XOR the pivot into the working vector
            basis[bit_len].add(basis[pivot_index]);

            if (basis[bit_len].next == 0) {
                // All top bits cleared — advance the generator
                basis[bit_len].next_state(bit_len);
                pivot_index = calc_1pos(basis[bit_len].next);
            } else {
                pivot_index = calc_1pos(basis[bit_len].next);
            }
        }

        // k(bit_len) = minimum count among the surviving basis vectors
        int min_count = basis[0].count;
        for (int i = 1; i < bit_len; i++) {
            if (min_count > basis[i].count)
                min_count = basis[i].count;
        }

        int t_v = std::min(min_count, kg / bit_len);
        result.ecart[bit_len] = kg / bit_len - t_v;
        result.se += result.ecart[bit_len];

        if (result.ecart[bit_len] > delta[bit_len] || result.se > mse) {
            maxl = bit_len;
            break;
        }

        // ── Adjust for next resolution (bit_len - 1) ───────────────────
        if (bit_len > 1) {
            int new_len = bit_len - 1;
            uint64_t mask = (new_len < 64) ? (~0ULL << (64 - new_len)) : ~0ULL;
            for (int i = 0; i <= bit_len; i++) {
                basis[i].next &= mask;
                if (basis[i].next == 0 && !basis[i].zero) {
                    basis[i].next_state(new_len);
                }
            }
        }
    }

    // Finalize
    result.se = 0;
    for (int l = 1; l <= maxl; l++) {
        if (result.ecart[l] == -1) result.ecart[l] = 0;
        result.se += result.ecart[l];
    }
    for (int l = maxl + 1; l <= maxL; l++) {
        if (result.ecart[l] == -1) result.ecart[l] = INT_MAX;
    }
    return result;
}

// ══════════════════════════════════════════════════════════════════════════
// compute_kv — PIS for a single resolution v.
// Same algorithm as test_me_harase but only for one resolution.
// Cost: O(k/v) generator steps.
// ══════════════════════════════════════════════════════════════════════════

int compute_kv(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int v)
{
    // Build v+1 basis vectors: v standard + 1 generator
    int size = v + 1;
    std::vector<LatticeVec> basis(size);

    // Standard basis: output bit i = 1 once, then 0
    for (int i = 0; i < v; i++) {
        basis[i] = LatticeVec::from_generators(gens, trans);
        basis[i].set_zero();
        basis[i].count = 0;
        basis[i].zero = false;
        basis[i].next = 1ULL << (63 - i);
    }

    // Generator vector: canonical init, advance one step
    basis[v] = LatticeVec::from_generators(gens, trans);
    {
        int pos = 0;
        for (int j = 0; j < basis[v].J; j++) {
            int ki = basis[v].gens[j]->k();
            BitVect e0(ki);
            e0.set_bit(0, 1);
            basis[v].gens[j]->init(e0);
            pos += ki;
        }
    }
    basis[v].next_state(v);

    // PIS reduction at resolution v
    int pivot_index = calc_1pos(basis[v].next);
    while (!basis[v].zero) {
        if (pivot_index == -1 || pivot_index >= v) break;

        if (basis[v].count > basis[pivot_index].count)
            std::swap(basis[v], basis[pivot_index]);

        basis[v].add(basis[pivot_index]);

        if (basis[v].next == 0) {
            basis[v].next_state(v);
            pivot_index = calc_1pos(basis[v].next);
        } else {
            pivot_index = calc_1pos(basis[v].next);
        }
    }

    // k(v) = min count among surviving basis vectors
    int min_count = basis[0].count;
    for (int i = 1; i < v; i++)
        if (min_count > basis[i].count)
            min_count = basis[i].count;

    return std::min(min_count, kg / v);
}

// ══════════════════════════════════════════════════════════════════════════
// PISCache — StackBase strategy for tempering optimization
//
// compute_all() runs the full PIS (v = L..1), saving a snapshot of the
// basis BEFORE each reduction.  restore_and_reduce(v) restores the
// snapshot at v, rebuilds the generator vector (which depends on the
// bitmask), and re-reduces.  Cost: O(k/v) per restore_and_reduce call.
// ══════════════════════════════════════════════════════════════════════════

struct PISCache::Impl {
    std::vector<Generator*> gens;
    std::vector<std::vector<Transformation*>> trans;
    int kg;
    int L;

    // Current basis (L+1 vectors)
    std::vector<LatticeVec> basis;

    // Cached snapshots: stack[v] = basis state BEFORE reducing at resolution v.
    // stack[v] has v+1 vectors (indices 0..v).
    std::vector<std::vector<LatticeVec>> stack;

    void init_basis() {
        int v = L;
        basis.resize(v + 1);

        // Standard basis vectors
        for (int i = 0; i < v; i++) {
            basis[i] = LatticeVec::from_generators(gens, trans);
            basis[i].set_zero();
            basis[i].count = 0;
            basis[i].zero = false;
            basis[i].next = 1ULL << (63 - i);
        }

        // Generator vector
        basis[v] = LatticeVec::from_generators(gens, trans);
        for (int j = 0; j < basis[v].J; j++) {
            int ki = basis[v].gens[j]->k();
            BitVect e0(ki);
            e0.set_bit(0, 1);
            basis[v].gens[j]->init(e0);
        }
        basis[v].next_state(v);
    }

    // Deep-copy the basis vectors 0..bit_len into a snapshot
    std::vector<LatticeVec> snapshot(int bit_len) {
        std::vector<LatticeVec> snap;
        snap.reserve(bit_len + 1);
        for (int i = 0; i <= bit_len; i++)
            snap.push_back(basis[i].deep_copy());
        return snap;
    }

    // Restore from a snapshot
    void restore(int bit_len, const std::vector<LatticeVec>& snap) {
        for (int i = 0; i <= bit_len && i < (int)snap.size(); i++)
            basis[i] = snap[i].deep_copy();
    }

    // PIS reduction at a single resolution
    int reduce_at(int bit_len) {
        int pivot_index = calc_1pos(basis[bit_len].next);
        while (!basis[bit_len].zero) {
            if (pivot_index == -1 || pivot_index >= bit_len) break;
            if (basis[bit_len].count > basis[pivot_index].count)
                std::swap(basis[bit_len], basis[pivot_index]);
            basis[bit_len].add(basis[pivot_index]);
            if (basis[bit_len].next == 0) {
                basis[bit_len].next_state(bit_len);
                pivot_index = calc_1pos(basis[bit_len].next);
            } else {
                pivot_index = calc_1pos(basis[bit_len].next);
            }
        }
        int min_count = basis[0].count;
        for (int i = 1; i < bit_len; i++)
            if (min_count > basis[i].count) min_count = basis[i].count;
        return std::min(min_count, kg / bit_len);
    }

    // Adjust basis from resolution bit_len to bit_len-1
    void adjust(int bit_len) {
        int new_len = bit_len - 1;
        uint64_t mask = (new_len < 64) ? (~0ULL << (64 - new_len)) : ~0ULL;
        for (int i = 0; i <= bit_len; i++) {
            basis[i].next &= mask;
            if (basis[i].next == 0 && !basis[i].zero)
                basis[i].next_state(new_len);
        }
    }

    // Rebuild the generator vector (basis[v]) from scratch for resolution v.
    // Called after a bitmask perturbation to reflect the new tempering.
    void rebuild_generator_vec(int v) {
        basis[v] = LatticeVec::from_generators(gens, trans);
        for (int j = 0; j < basis[v].J; j++) {
            int ki = basis[v].gens[j]->k();
            BitVect e0(ki);
            e0.set_bit(0, 1);
            basis[v].gens[j]->init(e0);
        }
        basis[v].next_state(v);
    }
};

PISCache::PISCache(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L)
    : gens_(gens), trans_(trans), kg_(kg), L_(L),
      impl_(std::make_shared<Impl>())
{
    impl_->gens = gens;
    impl_->trans = trans;
    impl_->kg = kg;
    impl_->L = L;
    impl_->stack.resize(L + 1);
}

std::vector<int> PISCache::compute_all() {
    auto& I = *impl_;
    I.init_basis();

    std::vector<int> ecart(L_ + 1, 0);

    for (int v = L_; v >= 1; v--) {
        // Save snapshot BEFORE reducing at resolution v
        I.stack[v] = I.snapshot(v);

        // Reduce
        int kv = I.reduce_at(v);
        ecart[v] = kg_ / v - kv;

        // Adjust for next resolution
        if (v > 1) I.adjust(v);
    }
    return ecart;
}

int PISCache::restore_and_reduce(int v) {
    auto& I = *impl_;

    // Restore the cached basis at resolution v (before reduction)
    if (v <= L_ && !I.stack[v].empty()) {
        I.restore(v, I.stack[v]);
    } else {
        // No cache — rebuild from scratch
        I.init_basis();
        for (int bl = L_; bl > v; bl--)
            I.adjust(bl);
    }

    // Rebuild the generator vector (depends on current bitmask)
    I.rebuild_generator_vec(v);

    // Reduce at resolution v
    return I.reduce_at(v);
}
