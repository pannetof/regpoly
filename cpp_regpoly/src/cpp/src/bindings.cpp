#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "bitvect.h"
#include "generateur.h"
#include "transformation.h"
#include "gauss.h"
#include "factory.h"
#include "lattice_polys.h"
#include "harase_lattice.h"
#include "lattice_optimizer.h"

#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>

namespace py = pybind11;

// ═══════════════════════════════════════════════════════════════════════════
// BitVect <-> Python int conversion helpers
// ═══════════════════════════════════════════════════════════════════════════

static py::int_ bitvect_to_pyint(const BitVect& bv) {
    if (bv.nbits() == 0)
        return py::int_(0);

    int nw = bv.nwords();
    int nbits = bv.nbits();

    py::int_ result(0);
    for (int i = 0; i < nw; i++) {
        uint64_t w = bv.data()[i];
        if (i == 0)
            result = py::int_(w);
        else
            result = (result << py::int_(64)) | py::int_(w);
    }

    int tail = nw * 64 - nbits;
    if (tail > 0)
        result = result >> py::int_(tail);

    return result;
}

static BitVect pyint_to_bitvect(int nbits, const py::int_& val) {
    BitVect bv(nbits);
    int nw = bv.nwords();
    int tail = nw * 64 - nbits;

    py::int_ shifted;
        if (tail > 0) shifted = val << py::int_(tail); else shifted = val;

    py::int_ mask64(0xFFFFFFFFFFFFFFFFULL);
    for (int i = nw - 1; i >= 0; i--) {
        bv.data()[i] = (shifted & mask64).cast<uint64_t>();
        shifted = shifted >> py::int_(64);
    }
    return bv;
}

// ═══════════════════════════════════════════════════════════════════════════
// py::dict -> Params conversion
// ═══════════════════════════════════════════════════════════════════════════

static Params dict_to_params(const py::dict& d) {
    Params p;
    for (auto item : d) {
        std::string key = item.first.cast<std::string>();
        py::handle val = item.second;

        // Handle "coeffs" specially: list of {value, position} dicts
        // → flatten to "coeff" (uint_vec) and "nocoeff" (int_vec)
        if (key == "coeffs") {
            auto coeffs = val.cast<py::list>();
            std::vector<uint64_t> coeff_vals;
            std::vector<int> nocoeff_vals;
            for (auto c : coeffs) {
                auto cd = c.cast<py::dict>();
                coeff_vals.push_back(cd["value"].cast<uint64_t>());
                nocoeff_vals.push_back(cd["position"].cast<int>());
            }
            p.set_uint_vec("coeff", coeff_vals);
            p.set_int_vec("nocoeff", nocoeff_vals);
            continue;
        }

        // Try bool first (before int, since Python bool is a subclass of int)
        if (py::isinstance<py::bool_>(val)) {
            p.set_bool(key, val.cast<bool>());
        } else if (py::isinstance<py::int_>(val)) {
            try {
                p.set_int(key, val.cast<int64_t>());
            } catch (...) {
                // Value exceeds int64_t range (e.g. large uint64 bitmask)
                // Store it as int64 with the same bit pattern
                uint64_t uval = val.cast<uint64_t>();
                p.set_int(key, static_cast<int64_t>(uval));
            }
        } else if (py::isinstance<py::str>(val)) {
            p.set_string(key, val.cast<std::string>());
        } else if (py::isinstance<py::list>(val)) {
            // Try as vector<int> first, then vector<uint64_t>
            try {
                p.set_int_vec(key, val.cast<std::vector<int>>());
            } catch (...) {
                p.set_uint_vec(key, val.cast<std::vector<uint64_t>>());
            }
        }
        // Skip other types silently (e.g. "type" key already handled by caller)
    }
    return p;
}

// ═══════════════════════════════════════════════════════════════════════════
// Module definition
// ═══════════════════════════════════════════════════════════════════════════

PYBIND11_MODULE(_regpoly_cpp, m) {
    m.doc() = "C++ acceleration module for regpoly PRNG analysis library";

    // ── BitVect ──────────────────────────────────────────────────────────

    py::class_<BitVect>(m, "BitVect")
        .def(py::init<>())
        .def(py::init<int>(), py::arg("nbits"))
        .def("nbits", &BitVect::nbits)
        .def("nwords", &BitVect::nwords)
        .def("get_bit", &BitVect::get_bit)
        .def("set_bit", &BitVect::set_bit)
        .def("get_word", &BitVect::get_word)
        .def("set_word", &BitVect::set_word)
        .def("xor_with", &BitVect::xor_with)
        .def("and_mask", &BitVect::and_mask)
        .def("and_invmask", &BitVect::and_invmask)
        .def("lshift", &BitVect::lshift)
        .def("rshift", &BitVect::rshift)
        .def("copy", &BitVect::copy)
        .def("copy_part_from", &BitVect::copy_part_from)
        .def("top_word", &BitVect::top_word)
        .def("zero", &BitVect::zero)
        .def("to_int", [](const BitVect& bv) -> py::int_ {
            return bitvect_to_pyint(bv);
        })
        .def_static("from_int", [](int nbits, const py::int_& val) -> BitVect {
            return pyint_to_bitvect(nbits, val);
        }, py::arg("nbits"), py::arg("val"))
        .def("__repr__", [](const BitVect& bv) {
            return "<BitVect nbits=" + std::to_string(bv.nbits()) + ">";
        });

    // ── Generateur ──────────────────────────────────────────────────────

    py::class_<Generateur, std::unique_ptr<Generateur>>(m, "Generateur")
        .def("name", &Generateur::name)
        .def("display_str", &Generateur::display_str)
        .def("k", &Generateur::k)
        .def("L", &Generateur::L)
        .def("init", &Generateur::init)
        .def("next", &Generateur::next)
        .def("char_poly", &Generateur::char_poly)
        .def("is_full_period", &Generateur::is_full_period)
        .def("transition_matrix", &Generateur::transition_matrix)
        .def("get_output", &Generateur::get_output)
        .def("copy", [](const Generateur& g) { return g.copy(); })
        .def("state", [](const Generateur& g) -> BitVect { return g.state().copy(); });

    // ── GaussMatrix ──────────────────────────────────────────────────────

    py::class_<GaussMatrix>(m, "GaussMatrix")
        .def(py::init<int, int>(), py::arg("nrows"), py::arg("ncols"))
        .def("copy", &GaussMatrix::copy)
        .def("nrows", &GaussMatrix::nrows)
        .def("ncols", &GaussMatrix::ncols)
        .def("nwords", &GaussMatrix::nwords)
        .def("bit_test", &GaussMatrix::bit_test)
        .def("swap_rows", &GaussMatrix::swap_rows)
        .def("row_xor", &GaussMatrix::row_xor)
        .def("eliminate_column", &GaussMatrix::eliminate_column)
        .def("eliminate_column_masked", &GaussMatrix::eliminate_column_masked)
        .def("find_pivot", &GaussMatrix::find_pivot)
        .def("dimension_equid", &GaussMatrix::dimension_equid,
             py::arg("kg"), py::arg("l"), py::arg("L"))
        .def("resolution_equid", &GaussMatrix::resolution_equid,
             py::arg("kg"), py::arg("t"), py::arg("L"), py::arg("indices"))
        .def("rang_cf", &GaussMatrix::rang_cf,
             py::arg("kg"), py::arg("t"), py::arg("l"), py::arg("L"))
        .def("set_row_from_words", [](GaussMatrix& mat, int row,
                                       const std::vector<uint64_t>& words) {
            mat.set_row_from_words(row, words.data(), (int)words.size());
        }, py::arg("row"), py::arg("words"))
        .def("set_row_from_int", [](GaussMatrix& mat, int row,
                                     const py::int_& val) {
            int nw = mat.nwords();
            int ncols = mat.ncols();
            int tail = nw * 64 - ncols;
            py::int_ shifted;
        if (tail > 0) shifted = val << py::int_(tail); else shifted = val;
            std::vector<uint64_t> words(nw, 0);
            py::int_ mask64(0xFFFFFFFFFFFFFFFFULL);
            for (int i = nw - 1; i >= 0; i--) {
                words[i] = (shifted & mask64).cast<uint64_t>();
                shifted = shifted >> py::int_(64);
            }
            mat.set_row_from_words(row, words.data(), nw);
        }, py::arg("row"), py::arg("val"));

    // ── Transformation ──────────────────────────────────────────────────

    py::class_<Transformation, std::unique_ptr<Transformation>>(m, "Transformation")
        .def("name", &Transformation::name)
        .def("display_str", &Transformation::display_str)
        .def("apply", &Transformation::apply)
        .def("w", &Transformation::w)
        .def("update", [](Transformation& t, const py::dict& d) {
            t.update(dict_to_params(d));
        })
        .def("copy", [](const Transformation& t) { return t.copy(); });

    // ── prepare_mat ──────────────────────────────────────────────────────

    m.def("prepare_mat",
          [](const py::list& gens_py,
             const std::vector<int>& gen_k,
             const py::list& trans_py,
             int kg, int indice_max, int L) {
        std::vector<Generateur*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generateur*>());
        std::vector<std::vector<Transformation*>> trans;
        for (auto comp_trans : trans_py) {
            std::vector<Transformation*> chain;
            for (auto t : comp_trans)
                chain.push_back(t.cast<Transformation*>());
            trans.push_back(chain);
        }
        return GaussMatrix::prepare(gens, gen_k, trans, kg, indice_max, L);
    }, py::arg("gens"), py::arg("gen_k"), py::arg("trans"),
       py::arg("kg"), py::arg("indice_max"), py::arg("L"));

    // ── test_me_lat (dual lattice method) ──────────────────────────────

    m.def("test_me_lat",
          [](const py::list& gens_py,
             const py::list& trans_py,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        std::vector<Generateur*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generateur*>());
        std::vector<std::vector<Transformation*>> trans;
        for (auto comp_trans : trans_py) {
            std::vector<Transformation*> chain;
            for (auto t : comp_trans)
                chain.push_back(t.cast<Transformation*>());
            trans.push_back(chain);
        }
        auto result = test_me_lat(gens, trans, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gens"), py::arg("trans"),
       py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    // ── test_me_lat_harase (primal lattice, Mulders-Storjohann) ────────

    m.def("test_me_lat_harase",
          [](const py::list& gens_py,
             const py::list& trans_py,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        std::vector<Generateur*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generateur*>());
        std::vector<std::vector<Transformation*>> trans;
        for (auto comp_trans : trans_py) {
            std::vector<Transformation*> chain;
            for (auto t : comp_trans)
                chain.push_back(t.cast<Transformation*>());
            trans.push_back(chain);
        }
        auto result = test_me_lat_harase(gens, trans, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gens"), py::arg("trans"),
       py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    // ── LatticeOptCache (dual lattice StackBase for optimizer) ─────────

    py::class_<LatticeOptCache>(m, "LatticeOptCache")
        .def(py::init([](const py::list& gens_py, const py::list& trans_py,
                         int kg, int L) {
            std::vector<Generateur*> gens;
            for (auto item : gens_py) gens.push_back(item.cast<Generateur*>());
            std::vector<std::vector<Transformation*>> trans;
            for (auto ct : trans_py) {
                std::vector<Transformation*> chain;
                for (auto t : ct) chain.push_back(t.cast<Transformation*>());
                trans.push_back(chain);
            }
            return LatticeOptCache(gens, trans, kg, L);
        }), py::arg("gens"), py::arg("trans"), py::arg("kg"), py::arg("L"))
        .def("compute_all", &LatticeOptCache::compute_all)
        .def("compute_gap", &LatticeOptCache::compute_gap, py::arg("v"))
        .def("refresh_inv_g0", &LatticeOptCache::refresh_inv_g0)
        .def("rebuild", &LatticeOptCache::rebuild)
        .def("reset_step", &LatticeOptCache::reset_step)
        .def("step", &LatticeOptCache::step, py::arg("v"));

    // ── PISCache (StackBase strategy for tempering optimizer) ──────────

    py::class_<PISCache>(m, "PISCache")
        .def(py::init([](const py::list& gens_py, const py::list& trans_py,
                         int kg, int L) {
            std::vector<Generateur*> gens;
            for (auto item : gens_py) gens.push_back(item.cast<Generateur*>());
            std::vector<std::vector<Transformation*>> trans;
            for (auto ct : trans_py) {
                std::vector<Transformation*> chain;
                for (auto t : ct) chain.push_back(t.cast<Transformation*>());
                trans.push_back(chain);
            }
            return PISCache(gens, trans, kg, L);
        }), py::arg("gens"), py::arg("trans"), py::arg("kg"), py::arg("L"))
        .def("compute_all", &PISCache::compute_all)
        .def("restore_and_reduce", &PISCache::restore_and_reduce,
             py::arg("v"))
        .def("kg", &PISCache::kg)
        .def("L", &PISCache::L);

    // ── compute_kv (single-resolution PIS) ────────────────────────────

    m.def("compute_kv",
          [](const py::list& gens_py,
             const py::list& trans_py,
             int kg, int v) -> int {
        std::vector<Generateur*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generateur*>());
        std::vector<std::vector<Transformation*>> trans;
        for (auto comp_trans : trans_py) {
            std::vector<Transformation*> chain;
            for (auto t : comp_trans)
                chain.push_back(t.cast<Transformation*>());
            trans.push_back(chain);
        }
        return compute_kv(gens, trans, kg, v);
    }, py::arg("gens"), py::arg("trans"),
       py::arg("kg"), py::arg("v"));

    // ── Generator subclass registrations ────────────────────────────────
    register_generator_types(m);

    // ── Factory functions ────────────────────────────────────────────────

    m.def("create_generator",
          [](const std::string& family, const py::dict& d, int L) {
        return create_generator(family, dict_to_params(d), L);
    }, py::arg("family"), py::arg("params"), py::arg("L"));

    m.def("create_transformation",
          [](const std::string& type, const py::dict& d) {
        return create_transformation(type, dict_to_params(d));
    }, py::arg("type"), py::arg("params"));

    // ── Primitivity test with precomputed factors ─────────────────────

    m.def("is_primitive_with_factors",
          [](const BitVect& char_poly_bv, int k,
             const std::vector<std::string>& factor_strings) -> bool {
        // Build NTL polynomial
        NTL::GF2X f;
        NTL::SetCoeff(f, k);
        for (int j = 0; j < k; j++)
            if (char_poly_bv.get_bit(j))
                NTL::SetCoeff(f, j);

        // Constant term must be 1
        if (!IsOne(coeff(f, 0)))
            return false;

        // Irreducibility check
        if (NTL::IterIrredTest(f) == 0)
            return false;

        // For each prime factor p of 2^k-1, check x^((2^k-1)/p) != 1
        NTL::GF2XModulus F;
        NTL::build(F, f);
        NTL::ZZ order = NTL::power(NTL::ZZ(2), k) - 1;

        for (const auto& s : factor_strings) {
            NTL::ZZ p = NTL::conv<NTL::ZZ>(s.c_str());
            NTL::ZZ exp = order / p;
            NTL::GF2X r;
            NTL::PowerXMod(r, exp, F);
            if (IsOne(r))
                return false;
        }
        return true;
    }, py::arg("char_poly"), py::arg("k"), py::arg("factors"));

    m.def("is_irreducible",
          [](const BitVect& char_poly_bv, int k) -> bool {
        NTL::GF2X f;
        NTL::SetCoeff(f, k);
        for (int j = 0; j < k; j++)
            if (char_poly_bv.get_bit(j))
                NTL::SetCoeff(f, j);
        if (!IsOne(coeff(f, 0)))
            return false;
        return NTL::IterIrredTest(f) != 0;
    }, py::arg("char_poly"), py::arg("k"));

    m.def("get_param_specs",
          [](const std::string& family) -> py::list {
        auto specs = get_param_specs(family);
        py::list result;
        for (auto& s : specs) {
            py::dict d;
            d["name"]       = s.name;
            d["type"]       = s.type;
            d["structural"] = s.structural;
            d["has_default"] = s.has_default;
            d["default"]    = s.default_val;
            d["rand_type"]  = s.rand_type;
            d["rand_args"]  = s.rand_args;
            result.append(d);
        }
        return result;
    }, py::arg("family"));
}
