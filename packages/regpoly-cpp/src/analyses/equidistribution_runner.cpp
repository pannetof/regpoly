#include "equidistribution_runner.h"

#include "combined.h"
#include "gauss.h"
#include "resolution_sets.h"
#include "transformation.h"

#include <climits>
#include <utility>
#include <vector>

namespace {

struct UnpackedGenerator {
    std::vector<Generator*> gens;
    std::vector<std::vector<Transformation*>> trans;
    std::vector<int> gen_k;
};

UnpackedGenerator unpack_for_kernel(const Generator& gen) {
    if (auto* cg = dynamic_cast<const CombinedGenerator*>(&gen)) {
        UnpackedGenerator u;
        u.gens = cg->raw_component_pointers();
        u.trans = cg->raw_tempering_pointers();
        u.gen_k.reserve(u.gens.size());
        for (auto* c : u.gens) u.gen_k.push_back(c->k());
        return u;
    }
    UnpackedGenerator u;
    u.gens = {const_cast<Generator*>(&gen)};
    u.trans = {{}};
    u.gen_k = {gen.k()};
    return u;
}

}  // namespace

MatricialEquidResult run_matricial_equidistribution(
    const Generator& gen,
    int kg, int L, int Lmax,
    const std::vector<int>& delta, int mse)
{
    auto u = unpack_for_kernel(gen);

    // 1. Build the GaussMatrix once. dimension_equid is destructive so
    //    each iteration works on a copy.
    GaussMatrix mat_full = GaussMatrix::prepare(
        u.gens, u.gen_k, u.trans, kg, /*indice_max=*/kg, L);

    std::vector<int> ecart(Lmax + 1, -1);
    int se = 0;

    std::vector<bool> psi12 = compute_psi12(kg, Lmax);

    bool verif = false;
    int maxl = Lmax;
    int l = 1;
    while (l <= Lmax) {
        if (ecart[l] == -1 && (psi12[l] || verif)) {
            int t = kg / l;
            GaussMatrix mat = mat_full.copy();
            int t_l = mat.dimension_equid(kg, l, L);
            ecart[l] = t - t_l;
            se += ecart[l];

            if (ecart[l] > delta[l] || se > mse) {
                maxl = l;
                break;
            }

            if (ecart[l] != 0) {
                verif = true;
                if (l != 1)
                    l -= 2;
            } else {
                verif = false;
            }
        }
        ++l;
    }

    se = 0;
    for (int i = 1; i <= maxl; ++i) {
        if (ecart[i] == -1) ecart[i] = 0;
        se += ecart[i];
    }
    for (int i = maxl + 1; i <= Lmax; ++i) {
        if (ecart[i] == -1) ecart[i] = INT_MAX;
    }

    return {std::move(ecart), se, true};
}

CollisionFreeResult run_collision_free(
    const Generator& gen,
    int kg, int L, int L_for_phi4)
{
    auto u = unpack_for_kernel(gen);
    std::vector<bool> phi4 = compute_phi4(kg, L_for_phi4);

    std::vector<int> ecart_cf(kg + 1, 0);
    int secf = 0;

    GaussMatrix mat_full = GaussMatrix::prepare(
        u.gens, u.gen_k, u.trans, kg, /*indice_max=*/kg, L);

    for (int t = kg; t >= 2; --t) {
        if (!phi4[t]) continue;
        int l = kg / t;
        GaussMatrix mat = mat_full.copy();
        int rank = mat.rang_cf(kg, t, l + 1, L);
        int gap = kg - rank;
        ecart_cf[l] = gap;
        secf += gap;
    }

    return {std::move(ecart_cf), secf, true};
}
