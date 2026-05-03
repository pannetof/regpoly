#include "mtgp.h"
#include <algorithm>
#include <cstdio>
#include <stdexcept>

MTGPGen::MTGPGen(int mexp, int pos, int sh1, int sh2, uint32_t mask,
           const std::vector<uint32_t>& tbl,
           const std::vector<uint32_t>& tmp_tbl,
           int L)
    : Generator(mexp, L),
      mexp_(mexp), N_(mexp / 32 + 1),
      pos_(pos), sh1_(sh1), sh2_(sh2),
      mask_(mask), tbl_(tbl), tmp_tbl_(tmp_tbl),
      idx_(0), last_output_(0)
{
    if (tbl_.size() != 16)
        throw std::invalid_argument("MTGPGen: tbl must have 16 entries");
    if (tmp_tbl_.size() != 16)
        throw std::invalid_argument("MTGPGen: tmp_tbl must have 16 entries");
    state_ = BitVect(32 * N_);
}

std::string MTGPGen::name() const { return "MTGPGen"; }

std::string MTGPGen::display_str() const {
    char buf[160];
    snprintf(buf, sizeof(buf),
             " MEXP=%d N=%d pos=%d sh1=%d sh2=%d mask=%08x",
             mexp_, N_, pos_, sh1_, sh2_, mask_);
    return std::string(buf);
}

void MTGPGen::init(const BitVect& init_bv) {
    state_ = BitVect(32 * N_);
    state_.copy_part_from(init_bv,
                          std::min(init_bv.nbits(), 32 * N_));
    idx_ = 0;
    last_output_ = 0;
    // Guarantee a non-zero state (minor analogue of the reference
    // implementation's initialisation function).
    bool any = false;
    for (int i = 0; i < N_ && !any; i++)
        if (get_word(i) != 0) any = true;
    if (!any) set_word(0, 1);
}

void MTGPGen::next() {
    // One recursion + tempering step on state[idx_].
    const uint32_t x1 = get_word(idx_);
    const uint32_t x2 = get_word((idx_ + 1) % N_);
    const uint32_t y_in = get_word((idx_ + pos_) % N_);
    const uint32_t t = get_word((idx_ + pos_ - 1 + N_) % N_);

    uint32_t y = (x1 & mask_) ^ x2 ^ (y_in << sh1_);
    const uint32_t mat = tbl_[y & 0x0f];
    const uint32_t new_word = y ^ (y >> sh2_) ^ mat;
    set_word(idx_, new_word);

    // Tempering: output = new_word ^ tmp_tbl[ (t ^ (t>>16)) & 0x0f ]
    const uint32_t tt = t ^ (t >> 16);
    last_output_ = new_word ^ tmp_tbl_[tt & 0x0f];

    idx_ = (idx_ + 1) % N_;
}

std::unique_ptr<Generator> MTGPGen::copy() const {
    auto p = std::make_unique<MTGPGen>(mexp_, pos_, sh1_, sh2_,
                                    mask_, tbl_, tmp_tbl_, L_);
    p->state_ = state_.copy();
    p->idx_ = idx_;
    p->last_output_ = last_output_;
    return p;
}

BitVect MTGPGen::get_output() const {
    // Return the 32-bit tempered output of the most recent next().
    BitVect out(32);
    out.set_word(0, 32, (uint64_t)last_output_);
    return out;
}

std::unique_ptr<Generator> MTGPGen::from_params(
    const Params& params, int L)
{
    const int mexp = (int)params.get_int("mexp");
    const int pos  = (int)params.get_int("pos");
    const int sh1  = (int)params.get_int("sh1");
    const int sh2  = (int)params.get_int("sh2");
    const uint32_t mask = (uint32_t)params.get_int("mask");
    auto tbl = params.get_uint_vec("tbl");
    auto tmp_tbl = params.get_uint_vec("tmp_tbl");
    if (tbl.size() != 16)
        throw std::invalid_argument(
            "MTGPGen: 'tbl' must be a 16-element uint_vec");
    if (tmp_tbl.size() != 16)
        throw std::invalid_argument(
            "MTGPGen: 'tmp_tbl' must be a 16-element uint_vec");
    std::vector<uint32_t> t32(16), tt32(16);
    for (int i = 0; i < 16; i++) {
        t32[i]  = (uint32_t)tbl[i];
        tt32[i] = (uint32_t)tmp_tbl[i];
    }
    return std::make_unique<MTGPGen>(mexp, pos, sh1, sh2, mask,
                                  t32, tt32, L);
}

std::vector<ParamSpec> MTGPGen::param_specs() {
    return {
        {"mexp",    "int",      true,  false, 0, "",     "", false},
        {"pos",     "int",      true,  false, 0, "",     "", false},
        {"sh1",     "int",      true,  false, 0, "",     "", false},
        {"sh2",     "int",      true,  false, 0, "",     "", false},
        {"mask",    "int",      true,  false, 0, "",     "", false},
        {"tbl",     "uint_vec", true,  false, 0, "none", "", false},
        {"tmp_tbl", "uint_vec", true,  false, 0, "none", "", false},
    };
}
