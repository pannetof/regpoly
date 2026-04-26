#pragma once
#include "generateur.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>

// Marsaglia 2003 XORSHIFT128.
//
// 128-bit state held as four 32-bit lanes (x, y, z, w) and a 32-bit
// output.  The recurrence is parametrised by a shift triple (a, b, c)
// with 1 ≤ a, b, c ≤ 31, plus a `pattern` selector that picks between
// two recurrence shapes shipped by MTToolBox samples/XORSHIFT.
//
// pattern = 0 (default — matches xorshift-2.cpp, "canonical" form):
//     t = x ^ (x << a)
//     x = y; y = z; z = w
//     w = (w ^ (w >> c)) ^ (t ^ (t >> b))
//
// pattern = 1 (matches xorshift-5.cpp, search-result form):
//     t = (x ^ (x << b)) ^ ((x ^ (x << b)) >> c)
//     x = y; y = z; z = w
//     w = (w ^ (w << a)) ^ t
//
// State layout in `state_` (BitVect, MSB-first; w first so default
// get_output() returns the L-bit prefix of `w`):
//     bits   0.. 31 = w
//     bits  32.. 63 = z
//     bits  64.. 95 = y
//     bits  96..127 = x
class XorShift128 : public Generateur {
public:
    XorShift128(int a, int b, int c, int L, int pattern = 0);

    static std::unique_ptr<Generateur> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;

    int a() const { return a_; }
    int b() const { return b_; }
    int c() const { return c_; }
    int pattern() const { return pattern_; }

private:
    int a_, b_, c_;
    int pattern_;

    void load_state(uint32_t& x, uint32_t& y, uint32_t& z, uint32_t& w) const;
    void store_state(uint32_t x, uint32_t y, uint32_t z, uint32_t w);
};
