#include "factory.h"
#include "gen_polylcg.h"
#include "gen_tausworthe.h"
#include "gen_tgfsr.h"
#include "gen_mt.h"
#include "gen_f2w_base.h"
#include "gen_f2w_polylcg.h"
#include "gen_f2w_lfsr.h"
#include "gen_matsumoto.h"
#include "gen_marsaxorshift.h"
#include "gen_ac1d.h"
#include "gen_carry2.h"
#include "trans_permutation.h"
#include "trans_temper_mk.h"
#include <algorithm>
#include <stdexcept>

#include <pybind11/pybind11.h>
namespace py = pybind11;

void register_generator_types(py::module_& m) {
    py::class_<PolyLCG, Generateur, std::unique_ptr<PolyLCG>>(m, "PolyLCG");
    py::class_<Tausworthe, Generateur, std::unique_ptr<Tausworthe>>(m, "Tausworthe");
    py::class_<TGFSRGen, Generateur, std::unique_ptr<TGFSRGen>>(m, "TGFSRGen");
    py::class_<MersenneTwister, Generateur, std::unique_ptr<MersenneTwister>>(m, "MersenneTwister");
    py::class_<GenF2wBase, Generateur, std::unique_ptr<GenF2wBase>>(m, "GenF2wBase");
    py::class_<GenF2wPolyLCG, GenF2wBase, std::unique_ptr<GenF2wPolyLCG>>(m, "GenF2wPolyLCG");
    py::class_<GenF2wLFSR, GenF2wBase, std::unique_ptr<GenF2wLFSR>>(m, "GenF2wLFSR");
    py::class_<MatsumotoGen, Generateur, std::unique_ptr<MatsumotoGen>>(m, "MatsumotoGen");
    py::class_<MarsaXorshiftGen, Generateur, std::unique_ptr<MarsaXorshiftGen>>(m, "MarsaXorshiftGen");
    py::class_<AC1DGen, Generateur, std::unique_ptr<AC1DGen>>(m, "AC1DGen");
    py::class_<Carry2Gen, Generateur, std::unique_ptr<Carry2Gen>>(m, "Carry2Gen");
}

std::unique_ptr<Generateur> create_generator(
    const std::string& family, const Params& params, int L)
{
    if (family == "polylcg") {
        int k = (int)params.get_int("k");
        auto poly_list = params.get_int_vec("poly");
        BitVect poly_bv(k);
        for (int e : poly_list)
            poly_bv.set_bit(k - e - 1, 1);
        return std::make_unique<PolyLCG>(k, poly_bv, L);

    } else if (family == "taus" || family == "taus2") {
        auto poly_list = params.get_int_vec("poly");
        std::vector<int> Q(poly_list);
        std::sort(Q.begin(), Q.end());
        int k = Q.back();
        int s = (int)params.get_int("s");
        bool quicktaus = params.get_bool("quicktaus", true);
        return std::make_unique<Tausworthe>(k, Q, s, quicktaus, L);

    } else if (family == "tgfsr") {
        int w = (int)params.get_int("w");
        int r = (int)params.get_int("r");
        int m = (int)params.get_int("m");
        uint64_t a_val = (uint64_t)params.get_int("a");
        int k = w * r;
        BitVect a_bv(k);
        if (k > 32) {
            for (int i = 0; i < 32; i++)
                if ((a_val >> (31 - i)) & 1)
                    a_bv.set_bit(i, 1);
        } else {
            for (int i = 0; i < k; i++)
                if ((a_val >> (k - 1 - i)) & 1)
                    a_bv.set_bit(i, 1);
        }
        return std::make_unique<TGFSRGen>(w, r, m, a_bv, std::min(w, L));

    } else if (family == "MT") {
        int w = (int)params.get_int("w");
        int r = (int)params.get_int("r");
        int m = (int)params.get_int("m");
        int p = (int)params.get_int("p", 0);
        uint64_t a = (uint64_t)params.get_int("a");
        return std::make_unique<MersenneTwister>(w, r, m, p, a, L);

    } else if (family == "genf2w") {
        int w = (int)params.get_int("w");
        int r = (int)params.get_int("r");
        uint64_t modM = (uint64_t)params.get_int("modM");
        bool normal_basis = params.get_bool("normal_basis", false);
        int step_count = (int)params.get_int("step", 1);
        auto nocoeff_vals = params.get_int_vec("nocoeff");
        auto coeff_vals = params.get_uint_vec("coeff");
        int nbcoeff = (int)nocoeff_vals.size();
        std::string type_str = params.get_string("type", "polylcg");

        if (type_str == "lfsr") {
            return std::make_unique<GenF2wLFSR>(
                w, r, nbcoeff, nocoeff_vals, coeff_vals,
                modM, normal_basis, step_count, L);
        } else {
            return std::make_unique<GenF2wPolyLCG>(
                w, r, nbcoeff, nocoeff_vals, coeff_vals,
                modM, normal_basis, step_count, L);
        }
    } else if (family == "matsumoto") {
        int type = (int)params.get_int("type");
        int n = (int)params.get_int("n");
        int m = (int)params.get_int("m");
        auto paramsint = params.get_int_vec("paramsint");
        auto paramsunsigned_u64 = params.get_uint_vec("paramsunsigned");
        std::vector<uint32_t> paramsunsigned;
        for (auto v : paramsunsigned_u64)
            paramsunsigned.push_back((uint32_t)v);
        return std::make_unique<MatsumotoGen>(type, n, m, paramsint, paramsunsigned, L);

    } else if (family == "marsaxorshift") {
        int type = (int)params.get_int("type");
        int w = (int)params.get_int("w", 32);
        int r = (int)params.get_int("r", 1);
        int m = (int)params.get_int("m", 0);

        MarsaXorshiftGen::Type1Params t1{};
        MarsaXorshiftGen::Type2xParams t2x;
        std::vector<MarsaXorshiftGen::Tap> taps;
        MarsaXorshiftGen::Type4Params t4;
        std::vector<MarsaXorshiftGen::MiEntry> mi;

        if (type == 1) {
            auto abc = params.get_int_vec("shifts");
            t1.a = abc.size() > 0 ? abc[0] : 0;
            t1.b = abc.size() > 1 ? abc[1] : 0;
            t1.c = abc.size() > 2 ? abc[2] : 0;
        } else if (type >= 21 && type <= 25) {
            t2x.p = params.get_int_vec("p");
            t2x.q = params.get_int_vec("q");
            t2x.p.resize(3, 0);
            t2x.q.resize(3, 0);
        } else if (type == 3) {
            auto tap_pos = params.get_int_vec("tap_positions");
            auto tap_shifts = params.get_int_vec("tap_shifts");
            for (size_t i = 0; i < tap_pos.size(); i++) {
                taps.push_back({tap_pos[i],
                    (i < tap_shifts.size()) ? tap_shifts[i] : 0});
            }
        } else if (type == 4) {
            t4.p = params.get_int_vec("p");
            t4.q = params.get_int_vec("q");
            t4.p.resize(2, 0);
            t4.q.resize(2, 0);
        } else if (type == 100) {
            auto mi_pos = params.get_int_vec("mi_positions");
            auto mi_shifts = params.get_int_vec("mi_shifts");
            auto mi_counts = params.get_int_vec("mi_counts");
            int offset = 0;
            for (size_t i = 0; i < mi_pos.size(); i++) {
                MarsaXorshiftGen::MiEntry entry;
                entry.position = mi_pos[i];
                int count = (i < mi_counts.size()) ? mi_counts[i] : 1;
                for (int j = 0; j < count && offset < (int)mi_shifts.size(); j++) {
                    entry.shifts.push_back(mi_shifts[offset++]);
                }
                mi.push_back(entry);
            }
        }

        return std::make_unique<MarsaXorshiftGen>(
            type, w, r, m, t1, t2x, taps, t4, mi, L);

    } else if (family == "AC1D") {
        int n = (int)params.get_int("n");
        auto flat = params.get_int_vec("matrix");
        std::vector<std::vector<int>> matrix(n, std::vector<int>(n, 0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i * n + j < (int)flat.size())
                    matrix[i][j] = flat[i * n + j];
        return std::make_unique<AC1DGen>(n, matrix, L);

    } else if (family == "carry") {
        int w = (int)params.get_int("w", 32);
        int r = (int)params.get_int("r");
        int p = (int)params.get_int("p");
        int m1 = (int)params.get_int("m1");
        int m2 = (int)params.get_int("m2");
        int m3 = (int)params.get_int("m3");

        auto mat_types = params.get_int_vec("mat_types");
        auto mat_pi = params.get_int_vec("mat_pi");
        auto mat_pu_u64 = params.get_uint_vec("mat_pu");

        std::vector<Carry2Gen::MatrixEntry> matrices(8);
        for (int j = 0; j < 8; j++) {
            matrices[j].type = (j < (int)mat_types.size()) ? mat_types[j] : 1;
            for (int x = 0; x < 3; x++) {
                int idx = j * 3 + x;
                matrices[j].paramsint[x] = (idx < (int)mat_pi.size()) ? mat_pi[idx] : 0;
                matrices[j].paramsulong[x] = (idx < (int)mat_pu_u64.size())
                    ? mat_pu_u64[idx] : 0;
            }
        }
        return std::make_unique<Carry2Gen>(w, r, p, m1, m2, m3, matrices, L);
    }
    throw std::invalid_argument("Unknown generator family: " + family);
}

std::unique_ptr<Transformation> create_transformation(
    const std::string& type, const Params& params)
{
    int w = (int)params.get_int("w");

    if (type == "permut") {
        int p = (int)params.get_int("p");
        int q = (int)params.get_int("q");
        return std::make_unique<PermutationTrans>(w, p, q);

    } else if (type == "tempMK" || type == "tempMK2") {
        int mk_type = (type == "tempMK2") ? 2 : 1;
        int eta = (int)params.get_int("eta");
        int mu = (int)params.get_int("mu");
        int u = (int)params.get_int("u", 0);
        int l = (int)params.get_int("l", 0);
        uint64_t b = (uint64_t)params.get_int("b", 0);
        uint64_t c = (uint64_t)params.get_int("c", 0);
        return std::make_unique<TemperMKTrans>(w, mk_type, eta, mu, u, l, b, c);
    }
    throw std::invalid_argument("Unknown transformation type: " + type);
}
