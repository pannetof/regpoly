#include "bitvect.h"
#include <algorithm>

int BitVect::get_bit(int i) const {
    int w = i / WL;
    int b = WL - 1 - (i % WL);
    return (words_[w] >> b) & 1;
}

void BitVect::set_bit(int i, int v) {
    int w = i / WL;
    int b = WL - 1 - (i % WL);
    if (v)
        words_[w] |= (1ULL << b);
    else
        words_[w] &= ~(1ULL << b);
}

uint64_t BitVect::get_word(int idx, int w) const {
    int bit_pos = idx * w;
    int start_word = bit_pos / WL;
    int start_bit_in_word = bit_pos % WL;
    int avail = WL - start_bit_in_word;

    if (w <= avail) {
        int shift = avail - w;
        return (words_[start_word] >> shift) & ((w == WL) ? ~0ULL : ((1ULL << w) - 1));
    } else {
        int from_first = avail;
        int from_second = w - from_first;
        uint64_t hi = words_[start_word] & ((1ULL << from_first) - 1);
        uint64_t lo = words_[start_word + 1] >> (WL - from_second);
        return (hi << from_second) | lo;
    }
}

void BitVect::set_word(int idx, int w, uint64_t val) {
    int bit_pos = idx * w;
    int start_word = bit_pos / WL;
    int start_bit_in_word = bit_pos % WL;
    int avail = WL - start_bit_in_word;

    uint64_t mask_w = (w == WL) ? ~0ULL : ((1ULL << w) - 1);
    val &= mask_w;

    if (w <= avail) {
        int shift = avail - w;
        uint64_t clear_mask = mask_w << shift;
        words_[start_word] = (words_[start_word] & ~clear_mask) | (val << shift);
    } else {
        int from_first = avail;
        int from_second = w - from_first;
        uint64_t hi_mask = (1ULL << from_first) - 1;
        uint64_t hi_val = val >> from_second;
        words_[start_word] = (words_[start_word] & ~hi_mask) | hi_val;
        uint64_t lo_val = val & ((1ULL << from_second) - 1);
        int lo_shift = WL - from_second;
        uint64_t lo_mask = ((1ULL << from_second) - 1) << lo_shift;
        words_[start_word + 1] = (words_[start_word + 1] & ~lo_mask) | (lo_val << lo_shift);
    }
}

void BitVect::xor_with(const BitVect& other) {
    int n = std::min((int)words_.size(), (int)other.words_.size());
    for (int i = 0; i < n; i++)
        words_[i] ^= other.words_[i];
}

void BitVect::and_mask(int t) {
    if (t <= 0) {
        zero();
        return;
    }
    if (t >= nbits_) return;
    int full_words = t / WL;
    int rem = t % WL;
    if (rem != 0) {
        words_[full_words] &= (~0ULL) << (WL - rem);
        full_words++;
    }
    for (int i = full_words; i < (int)words_.size(); i++)
        words_[i] = 0;
}

void BitVect::and_invmask(int t) {
    if (t <= 0) return;
    if (t >= nbits_) {
        zero();
        return;
    }
    int full_words = t / WL;
    int rem = t % WL;
    for (int i = 0; i < full_words; i++)
        words_[i] = 0;
    if (rem != 0) {
        words_[full_words] &= (~0ULL) >> rem;
    }
}

void BitVect::lshift(int k) {
    if (k <= 0) return;
    if (k >= nbits_) {
        zero();
        return;
    }
    int nw = (int)words_.size();
    int word_shift = k / WL;
    int bit_shift = k % WL;

    if (bit_shift == 0) {
        for (int i = 0; i < nw - word_shift; i++)
            words_[i] = words_[i + word_shift];
        for (int i = nw - word_shift; i < nw; i++)
            words_[i] = 0;
    } else {
        int comp = WL - bit_shift;
        for (int i = 0; i < nw - word_shift - 1; i++)
            words_[i] = (words_[i + word_shift] << bit_shift)
                       | (words_[i + word_shift + 1] >> comp);
        words_[nw - word_shift - 1] = words_[nw - 1] << bit_shift;
        for (int i = nw - word_shift; i < nw; i++)
            words_[i] = 0;
    }
    clear_tail();
}

void BitVect::rshift(int k) {
    if (k <= 0) return;
    if (k >= nbits_) {
        zero();
        return;
    }
    int nw = (int)words_.size();
    int word_shift = k / WL;
    int bit_shift = k % WL;

    if (bit_shift == 0) {
        for (int i = nw - 1; i >= word_shift; i--)
            words_[i] = words_[i - word_shift];
        for (int i = 0; i < word_shift; i++)
            words_[i] = 0;
    } else {
        int comp = WL - bit_shift;
        for (int i = nw - 1; i > word_shift; i--)
            words_[i] = (words_[i - word_shift] >> bit_shift)
                       | (words_[i - word_shift - 1] << comp);
        words_[word_shift] = words_[0] >> bit_shift;
        for (int i = 0; i < word_shift; i++)
            words_[i] = 0;
    }
    clear_tail();
}

BitVect BitVect::copy() const {
    BitVect bv;
    bv.nbits_ = nbits_;
    bv.words_ = words_;
    return bv;
}

void BitVect::copy_part_from(const BitVect& other, int l) {
    if (l <= 0) {
        zero();
        return;
    }
    zero();
    int copy_bits = std::min(l, other.nbits_);
    int full_words = copy_bits / WL;
    int rem = copy_bits % WL;
    for (int i = 0; i < full_words && i < (int)words_.size() && i < (int)other.words_.size(); i++)
        words_[i] = other.words_[i];
    if (rem != 0 && full_words < (int)words_.size() && full_words < (int)other.words_.size()) {
        uint64_t mask = (~0ULL) << (WL - rem);
        words_[full_words] = other.words_[full_words] & mask;
    }
}

void BitVect::zero() {
    for (auto& w : words_)
        w = 0;
}

void BitVect::set_from_words(const uint64_t* data, int nwords) {
    int n = std::min(nwords, (int)words_.size());
    for (int i = 0; i < n; i++)
        words_[i] = data[i];
    for (int i = n; i < (int)words_.size(); i++)
        words_[i] = 0;
    clear_tail();
}

void BitVect::clear_tail() {
    int rem = nbits_ % WL;
    if (rem != 0 && !words_.empty()) {
        words_.back() &= (~0ULL) << (WL - rem);
    }
}
