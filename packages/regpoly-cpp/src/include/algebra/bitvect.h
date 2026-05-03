#pragma once
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>

class BitVect {
    std::vector<uint64_t> words_;
    int nbits_;
    static constexpr int WL = 64;

public:
    BitVect() : nbits_(0) {}

    explicit BitVect(int nbits)
        : words_((nbits + WL - 1) / WL, 0ULL), nbits_(nbits) {}

    int nbits() const { return nbits_; }
    int nwords() const { return (int)words_.size(); }
    const uint64_t* data() const { return words_.data(); }
    uint64_t* data() { return words_.data(); }

    uint64_t top_word() const {
        if (words_.empty()) return 0;
        return words_[0];
    }

    // ── Bit access (bit 0 = MSB) ────────────────────────────────────────
    int get_bit(int i) const;
    void set_bit(int i, int v);

    // ── Word access: idx-th w-bit block from the left ────────────────────
    uint64_t get_word(int idx, int w) const;
    void set_word(int idx, int w, uint64_t val);

    // ── Bitwise in-place ─────────────────────────────────────────────────
    void xor_with(const BitVect& other);
    void and_mask(int t);
    void and_invmask(int t);

    // ── Shifts in-place ──────────────────────────────────────────────────
    void lshift(int k);
    void rshift(int k);

    // ── Copy ─────────────────────────────────────────────────────────────
    BitVect copy() const;
    void copy_part_from(const BitVect& other, int l);

    // ── Zero out ─────────────────────────────────────────────────────────
    void zero();

    // ── Set from flat array ──────────────────────────────────────────────
    void set_from_words(const uint64_t* data, int nwords);

private:
    void clear_tail();
};
