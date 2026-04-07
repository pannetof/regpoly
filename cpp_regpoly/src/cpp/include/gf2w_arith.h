#pragma once
#include <cstdint>
#include <vector>

uint64_t gf2w_multiply_z(uint64_t a, int k, uint64_t modM);
uint64_t gf2w_multiply_poly(uint64_t a, uint64_t b, int w, uint64_t modM);
uint64_t gf2w_multiply_normal(uint64_t a, uint64_t b, int w,
                               const uint64_t* table);
std::vector<uint64_t> gf2w_make_table(int w, uint64_t modM);
