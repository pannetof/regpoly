// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>

/**
 * @file bitvect.h
 * @brief F_2 bit-vector primitive — the bottom of the algebra stack.
 * @ingroup core
 *
 * Defines `BitVect`, the F_2 vector primitive every higher-level
 * REGPOLY algebra type (`GaussMatrix`, `PolVect`, `BitMatrix` in the
 * Python wrapper, …) is built on top of. Bit ordering is MSB-first:
 * bit index `0` is the most significant bit of word 0.
 */

namespace regpoly::core {

/**
 * @brief F_2 bit vector with shift / XOR / slice / word-level access.
 *
 * Width-`nbits` vector stored as a `std::vector<uint64_t>` of
 * `ceil(nbits / 64)` words. Bit `0` is the MSB of word `0`. The
 * tail bits beyond `nbits - 1` in the last word are kept zeroed by
 * the mutating operations.
 *
 * @code{.cpp}
 *   using regpoly::core::BitVect;
 *   BitVect a(64);
 *   a.set_bit(0, 1);             // MSB
 *   a.set_bit(63, 1);            // LSB
 *   uint64_t high = a.get_word(0, 32);   // top 32 bits as a uint64_t
 *   BitVect b = a.copy();
 *   b.xor_with(a);               // b becomes zero
 * @endcode
 *
 * @see :py:class:`regpoly.core.bitvect.BitVect`
 *
 * @ingroup core
 */
class BitVect {
    std::vector<uint64_t> words_;
    int nbits_;
    static constexpr int WL = 64;

public:
    /** @brief Construct an empty (width-zero) vector. */
    BitVect() : nbits_(0) {}

    /**
     * @brief Construct an all-zero vector of width `nbits`.
     * @param nbits  Width in bits.
     */
    explicit BitVect(int nbits)
        : words_((nbits + WL - 1) / WL, 0ULL), nbits_(nbits) {}

    /** @brief Width in bits. */
    int nbits() const { return nbits_; }
    /** @brief Number of 64-bit words backing the vector. */
    int nwords() const { return (int)words_.size(); }
    /** @brief Pointer to the underlying word storage (const). */
    const uint64_t* data() const { return words_.data(); }
    /** @brief Pointer to the underlying word storage (mutable). */
    uint64_t* data() { return words_.data(); }

    /** @brief First (MSB) 64-bit word; `0` if the vector is empty. */
    uint64_t top_word() const {
        if (words_.empty()) return 0;
        return words_[0];
    }

    // ── Bit access (bit 0 = MSB) ────────────────────────────────────────
    /**
     * @brief Read the bit at index `i` (MSB-first).
     * @param i  Bit index in `[0, nbits())`.
     * @return   `0` or `1`.
     */
    int get_bit(int i) const;
    /**
     * @brief Write the bit at index `i` (MSB-first).
     * @param i  Bit index in `[0, nbits())`.
     * @param v  `0` or `1`.
     */
    void set_bit(int i, int v);

    // ── Word access: idx-th w-bit block from the left ────────────────────
    /**
     * @brief Read the `idx`-th `w`-bit block (MSB-first).
     * @param idx  Block index (each block is `w` bits wide).
     * @param w    Block width in bits (`<= 64`).
     * @return     The block packed as a `uint64_t` (LSB-aligned).
     */
    uint64_t get_word(int idx, int w) const;
    /**
     * @brief Write the `idx`-th `w`-bit block (MSB-first).
     * @param idx  Block index.
     * @param w    Block width in bits (`<= 64`).
     * @param val  Block value (LSB-aligned).
     */
    void set_word(int idx, int w, uint64_t val);

    // ── Bitwise in-place ─────────────────────────────────────────────────
    /**
     * @brief XOR `other` into `*this`.
     * @param other  Vector to XOR in.
     */
    void xor_with(const BitVect& other);
    /**
     * @brief AND `other` into `*this`.
     * @param other  Vector to AND in.
     */
    void and_with(const BitVect& other);
    /**
     * @brief AND `*this` with the mask `0..0 1..1` (leftmost `t` ones).
     * @param t  Number of leading 1 bits in the mask.
     */
    void and_mask(int t);
    /**
     * @brief AND `*this` with the inverse mask (leftmost `t` zeros).
     * @param t  Number of leading 0 bits in the mask.
     */
    void and_invmask(int t);

    // ── Shifts in-place ──────────────────────────────────────────────────
    /**
     * @brief Logical left shift in place.
     * @param k  Shift amount in bits.
     */
    void lshift(int k);
    /**
     * @brief Logical right shift in place.
     * @param k  Shift amount in bits.
     */
    void rshift(int k);

    // ── Copy ─────────────────────────────────────────────────────────────
    /**
     * @brief Return an independent deep copy.
     * @return  Freshly-allocated copy of `*this`.
     */
    BitVect copy() const;
    /**
     * @brief Copy the leading `l` bits from `other` into `*this`.
     * @param other  Source vector.
     * @param l      Number of leading bits to copy.
     */
    void copy_part_from(const BitVect& other, int l);

    /**
     * @brief Extract a `n_bits`-wide slice starting at `start_bit` (MSB-first).
     *
     * Used by the MarsaXorshift wide path when `w > 64`.
     *
     * @param start_bit  Start bit (`0` = MSB).
     * @param n_bits     Slice width in bits.
     * @return           A new `BitVect` of width `n_bits`.
     */
    BitVect get_slice(int start_bit, int n_bits) const;
    /**
     * @brief Overwrite `slice.nbits()` bits starting at `start_bit`.
     *
     * @param start_bit  Start bit (`0` = MSB).
     * @param slice      Slice to write.
     */
    void    set_slice(int start_bit, const BitVect& slice);

    // ── Zero out ─────────────────────────────────────────────────────────
    /** @brief Zero every bit (preserves width). */
    void zero();

    // ── Set from flat array ──────────────────────────────────────────────
    /**
     * @brief Bulk-load words from a raw buffer (width stays unchanged).
     * @param data    Source word buffer.
     * @param nwords  Number of words to copy.
     */
    void set_from_words(const uint64_t* data, int nwords);

private:
    void clear_tail();
};

}  // namespace regpoly::core
