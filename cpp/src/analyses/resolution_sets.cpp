// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "resolution_sets.h"

#include <algorithm>
#include <cmath>

using namespace regpoly::core;


namespace regpoly::core {

namespace {

int isqrt(int n) {
    if (n < 0) return 0;
    int r = static_cast<int>(std::sqrt(static_cast<double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

}  // namespace

std::vector<bool> compute_psi12(int kg, int L) {
    std::vector<bool> psi(L + 1, false);
    int r = isqrt(kg);
    int upper = std::min(r, L);
    for (int l = 1; l <= upper; ++l) psi[l] = true;
    int r2 = isqrt(kg - 1);
    int m = (L > 0) ? kg / L : 2;
    if (m < 2) m = 2;
    for (int t = m; t <= r2; ++t) {
        int idx = std::min(kg / t, L);
        psi[idx] = true;
    }
    return psi;
}

std::vector<bool> compute_phi4(int kg, int L) {
    std::vector<bool> phi(kg + 1, false);
    for (int t = 2; t <= kg; ++t) {
        if (kg % t) {
            int lt = kg / t;
            if (lt < L && (kg % (lt + 1)) && (kg / (t - 1)) > lt)
                phi[t] = true;
        }
    }
    return phi;
}

}  // namespace regpoly::core
