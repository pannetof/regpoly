// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "bitvect.h"
#include "combination.h"
#include "combined.h"
#include "equidistribution_method.h"
#include "equidistribution_runner.h"
#include "generator.h"
#include "transformation.h"
#include "gauss.h"
#include "factory.h"
#include "me_helpers.h"
#include "me_harase.h"
#include "me_notprimitive.h"
#include "me_notprimitive_simd.h"
#include "primitive_search.h"
#include "primitivity.h"
#include "resolution_sets.h"
#include "search_types.h"
#include "seek_search.h"
#include "tempering_search.h"
#include "catalog.h"
#include "legacy_reader.h"
#include "temper_optimizer.h"
#include "tempering_optimizer.h"
#include "tausworthe.h"
#include "tuplets_runner.h"

#include <NTL/GF2X.h>
#include <cctype>
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
            // Hex/decimal literals stored as strings in YAML (e.g. MT's
            // 'a' = '0x9908B0DF') must convert to int — generator
            // factories call params.get_int(key). Mirrors the same
            // conversion done by catalog.cpp's String case so the
            // Python and direct-YAML paths agree.
            std::string s = val.cast<std::string>();
            bool parsed = false;
            try {
                if (s.size() > 2 && s[0] == '0'
                    && (s[1] == 'x' || s[1] == 'X')) {
                    uint64_t u = std::stoull(s.substr(2), nullptr, 16);
                    p.set_int(key, static_cast<int64_t>(u));
                    parsed = true;
                } else if (!s.empty()
                           && (std::isdigit(static_cast<unsigned char>(s[0]))
                               || s[0] == '-' || s[0] == '+')) {
                    p.set_int(key, std::stoll(s));
                    parsed = true;
                }
            } catch (...) {
                parsed = false;
            }
            if (!parsed) p.set_string(key, s);
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
// Params -> py::dict conversion (mirror of dict_to_params)
// ═══════════════════════════════════════════════════════════════════════════

static py::dict params_to_dict(const Params& p) {
    py::dict d;
    for (const auto& kv : p.ints())      d[py::str(kv.first)] = py::int_(kv.second);
    for (const auto& kv : p.bools())     d[py::str(kv.first)] = py::bool_(kv.second);
    for (const auto& kv : p.strings())   d[py::str(kv.first)] = py::str(kv.second);
    for (const auto& kv : p.int_vecs())  d[py::str(kv.first)] = py::cast(kv.second);
    for (const auto& kv : p.uint_vecs()) d[py::str(kv.first)] = py::cast(kv.second);
    return d;
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

    // ── Generator ──────────────────────────────────────────────────────

    py::class_<Generator, std::unique_ptr<Generator>>(m, "Generator")
        .def("name", &Generator::name)
        .def("display_str", &Generator::display_str)
        .def("k", &Generator::k)
        .def("L", &Generator::L)
        .def("init", &Generator::init)
        .def("next", &Generator::next)
        .def("char_poly", &Generator::char_poly)
        .def("is_full_period", &Generator::is_full_period)
        .def("transition_matrix", &Generator::transition_matrix)
        .def("get_output", &Generator::get_output)
        .def("copy", [](const Generator& g) { return g.copy(); })
        .def("state", [](const Generator& g) -> BitVect { return g.state().copy(); })
        .def("default_test_method", &Generator::default_test_method, py::arg("test_type"));

    // Backwards-compat: legacy French class name still resolves to the same type.
    m.attr("Generateur") = m.attr("Generator");

    // ── CombinedGenerator ───────────────────────────────────────────────
    //
    // First-class Generator subclass that owns J components and presents
    // itself as a single Generator. The factory used here transfers
    // ownership: the input lists must contain pybind11-owned unique_ptrs
    // (the wrapped C++ objects move into the CombinedGenerator and the
    // Python wrapper around the source object becomes inert).

    py::class_<CombinedGenerator, Generator, std::unique_ptr<CombinedGenerator>>(
        m, "CombinedGenerator")
        .def(py::init([](py::list py_components, int Lmax) {
                 std::vector<std::unique_ptr<Generator>> comps;
                 comps.reserve(py_components.size());
                 for (auto h : py_components) {
                     auto& g = h.cast<Generator&>();
                     comps.push_back(g.copy());
                 }
                 return std::make_unique<CombinedGenerator>(
                     std::move(comps), Lmax);
             }),
             py::arg("components"), py::arg("Lmax"),
             "Construct from a list of Generator objects (each is copied) "
             "and a max output resolution Lmax. The combined output is the "
             "XOR of the components' L-bit outputs.")
        .def(py::init([](py::list py_components,
                          py::list py_tempering,
                          int Lmax) {
                 if (py::len(py_components) != py::len(py_tempering))
                     throw std::invalid_argument(
                         "components and tempering chains must have equal "
                         "length");
                 std::vector<std::unique_ptr<Generator>> comps;
                 std::vector<CombinedGenerator::ComponentTempering> chains;
                 comps.reserve(py_components.size());
                 chains.reserve(py_tempering.size());
                 for (size_t j = 0; j < py_components.size(); ++j) {
                     auto& g = py_components[j].cast<Generator&>();
                     comps.push_back(g.copy());
                     CombinedGenerator::ComponentTempering chain;
                     for (auto th : py_tempering[j].cast<py::list>()) {
                         auto& t = th.cast<Transformation&>();
                         chain.push_back(t.copy());
                     }
                     chains.push_back(std::move(chain));
                 }
                 return std::make_unique<CombinedGenerator>(
                     std::move(comps), std::move(chains), Lmax);
             }),
             py::arg("components"), py::arg("tempering_chains"),
             py::arg("Lmax"),
             "Construct from a list of Generator objects, a parallel list of "
             "per-component tempering chains (each a list of Transformation "
             "objects), and a max output resolution Lmax.")
        .def("J", &CombinedGenerator::J)
        .def("prefix_k", &CombinedGenerator::prefix_k);

    // ── Component + Combination (Phase 2.4b-pre) ───────────────────────
    //
    // Iterator over the cartesian product of J component generator
    // pools, with identity-uniqueness and shared-pool C(n,k) selection.
    // Replaces (in C++) the iteration engine of
    // regpoly.core.{component,combination}.

    py::class_<Component, std::shared_ptr<Component>>(m, "Component")
        .def(py::init<>())
        .def("nb_gen",   &Component::nb_gen)
        .def("nb_trans", &Component::nb_trans)
        .def("current_gen", &Component::current_gen)
        .def("set_current_gen", &Component::set_current_gen)
        .def("add_gen",
             [](Component& c, const Generator& g) { c.add_gen(g); },
             py::arg("gen"),
             "Append a deep copy of `gen` to this component's pool.")
        .def("add_trans",
             [](Component& c, const Transformation& t) { c.add_trans(t); },
             py::arg("trans"),
             "Append a deep copy of `trans` to the tempering chain.")
        .def("share_pool_with", &Component::share_pool_with,
             py::arg("other"),
             "Reference `other`'s pool by shared_ptr (enables C(n,k) "
             "selection across slots).")
        .def("active_gen",
             [](Component& c) -> Generator* { return &c.active_gen(); },
             py::return_value_policy::reference_internal)
        .def("gen_at",
             [](Component& c, int i) -> Generator* { return &c.gen_at(i); },
             py::arg("i"),
             py::return_value_policy::reference_internal)
        .def("trans_at",
             [](Component& c, int i) -> Transformation* {
                 return &c.trans_at(i);
             },
             py::arg("i"),
             py::return_value_policy::reference_internal)
        .def("display", &Component::display);

    py::class_<Combination, std::shared_ptr<Combination>>(m, "Combination")
        .def(py::init<int, int>(), py::arg("J"), py::arg("Lmax"))
        .def_property_readonly("J",    &Combination::J)
        .def_property_readonly("Lmax", &Combination::Lmax)
        .def_property_readonly("k_g",  &Combination::k_g)
        .def_property_readonly("L",    &Combination::L)
        .def("component",
             [](Combination& c, int j) -> Component* { return &c.component(j); },
             py::arg("j"),
             py::return_value_policy::reference_internal)
        .def("__getitem__",
             [](Combination& c, int j) -> Generator* { return &c.at(j); },
             py::arg("j"),
             py::return_value_policy::reference_internal)
        .def("reset",     &Combination::reset)
        .def("next",      &Combination::next)
        .def("exhausted", &Combination::exhausted);

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
        std::vector<Generator*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generator*>());
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
        std::vector<Generator*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generator*>());
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

    // ── test_me_harase (primal lattice, Mulders-Storjohann) ────────

    m.def("test_me_harase",
          [](const py::list& gens_py,
             const py::list& trans_py,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        std::vector<Generator*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generator*>());
        std::vector<std::vector<Transformation*>> trans;
        for (auto comp_trans : trans_py) {
            std::vector<Transformation*> chain;
            for (auto t : comp_trans)
                chain.push_back(t.cast<Transformation*>());
            trans.push_back(chain);
        }
        auto result = test_me_harase(gens, trans, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gens"), py::arg("trans"),
       py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    // ── test_me_notprimitive (matricial DE on the BM-recovered
    //                         invariant subspace; no full-period
    //                         assumption) ───────────────────────────────

    m.def("test_me_notprimitive",
          [](const py::list& gens_py,
             const py::list& trans_py,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        std::vector<Generator*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generator*>());
        std::vector<std::vector<Transformation*>> trans;
        for (auto comp_trans : trans_py) {
            std::vector<Transformation*> chain;
            for (auto t : comp_trans)
                chain.push_back(t.cast<Transformation*>());
            trans.push_back(chain);
        }
        auto result = test_me_notprimitive(gens, trans, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gens"), py::arg("trans"),
       py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    // ── test_me_notprimitive_simd (SIMD-aware variant for SFMTGen etc.) ─

    m.def("test_me_notprimitive_simd",
          [](const py::list& gens_py,
             const py::list& trans_py,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        std::vector<Generator*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generator*>());
        std::vector<std::vector<Transformation*>> trans;
        for (auto comp_trans : trans_py) {
            std::vector<Transformation*> chain;
            for (auto t : comp_trans)
                chain.push_back(t.cast<Transformation*>());
            trans.push_back(chain);
        }
        auto result = test_me_notprimitive_simd(gens, trans, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gens"), py::arg("trans"),
       py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    // ── Phase 1 single-Generator& adapter overloads ────────────────────
    //
    // New public API: every kernel that today takes (gens, trans) lists
    // also accepts a single Generator& (typically a CombinedGenerator).
    // Suffix `_gen` distinguishes the new shape from the legacy
    // list-based bindings; Phase 2 will rewrite the Python wrapper layer
    // to call these and the legacy ones can then retire.

    m.def("test_me_lat_gen",
          [](const Generator& gen,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_lat(gen, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("test_me_harase_gen",
          [](const Generator& gen,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_harase(gen, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("test_me_notprimitive_gen",
          [](const Generator& gen,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_notprimitive(gen, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("test_me_notprimitive_simd_gen",
          [](const Generator& gen,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_notprimitive_simd(gen, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("compute_kv_gen",
          [](const Generator& gen, int kg, int v) -> int {
        return compute_kv(gen, kg, v);
    }, py::arg("gen"), py::arg("kg"), py::arg("v"));

    // ── Resolution-set helpers (Phase 2.1) ─────────────────────────────
    //
    // Returns a list of bools indexed 0..L (psi12) or 0..kg (phi4),
    // where out[r] == True iff resolution r is in the set.

    m.def("compute_psi12", &compute_psi12,
          py::arg("kg"), py::arg("L"),
          "Resolutions l in {1..L} whose dimension-equidistribution gap "
          "must be tested for a combined generator with state size kg.");

    m.def("compute_phi4", &compute_phi4,
          py::arg("kg"), py::arg("L"),
          "Dimensions t in {2..kg} whose collision-free rank must be "
          "checked for a combined generator with state size kg.");

    // ── Equidistribution / collision-free runners (Phase 2.3) ──────────
    //
    // Free functions that own the outer test-orchestration loop in C++.
    // The Python EquidistributionTest / CollisionFreeTest classes are
    // now thin wrappers around these.

    m.def("run_matricial_equidistribution",
          [](const Generator& gen, int kg, int L, int Lmax,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto r = run_matricial_equidistribution(gen, kg, L, Lmax, delta, mse);
        py::dict d;
        d["ecart"] = r.ecart;
        d["se"] = r.se;
        d["verified"] = r.verified;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("Lmax"),
       py::arg("delta"), py::arg("mse"));

    m.def("run_collision_free",
          [](const Generator& gen, int kg, int L, int L_for_phi4) -> py::dict {
        auto r = run_collision_free(gen, kg, L, L_for_phi4);
        py::dict d;
        d["ecart_cf"] = r.ecart_cf;
        d["secf"] = r.secf;
        d["verified"] = r.verified;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("L_for_phi4"));

    m.def("run_tuplets",
          [](const Generator& gen, int kg, int L, int tupd,
             const std::vector<int>& tuph, double threshold,
             int testtype) -> py::dict {
        auto r = run_tuplets(gen, kg, L, tupd, tuph, threshold, testtype);
        py::dict d;
        d["tupd"] = r.tupd;
        d["tuph"] = r.tuph;
        d["gap"] = r.gap;
        d["DELTA"] = r.DELTA;
        d["pourcentage"] = r.pourcentage;
        d["firstpart_max"] = r.firstpart_max;
        d["firstpart_sum"] = r.firstpart_sum;
        d["secondpart_max"] = r.secondpart_max;
        d["secondpart_sum"] = r.secondpart_sum;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("tupd"),
       py::arg("tuph"), py::arg("threshold"), py::arg("testtype"));

    // ── TemperOptCache (dual lattice StackBase for optimizer) ─────────

    py::class_<TemperOptCache>(m, "TemperOptCache")
        .def(py::init([](const py::list& gens_py, const py::list& trans_py,
                         int kg, int L) {
            std::vector<Generator*> gens;
            for (auto item : gens_py) gens.push_back(item.cast<Generator*>());
            std::vector<std::vector<Transformation*>> trans;
            for (auto ct : trans_py) {
                std::vector<Transformation*> chain;
                for (auto t : ct) chain.push_back(t.cast<Transformation*>());
                trans.push_back(chain);
            }
            return TemperOptCache(gens, trans, kg, L);
        }), py::arg("gens"), py::arg("trans"), py::arg("kg"), py::arg("L"))
        // Phase 1 single-Generator& constructor.
        .def(py::init([](const Generator& gen, int kg, int L) {
            return TemperOptCache(gen, kg, L);
        }), py::arg("gen"), py::arg("kg"), py::arg("L"))
        .def("compute_all", &TemperOptCache::compute_all)
        .def("compute_gap", &TemperOptCache::compute_gap, py::arg("v"))
        .def("refresh_inv_g0", &TemperOptCache::refresh_inv_g0)
        .def("rebuild", &TemperOptCache::rebuild)
        .def("reset_step", &TemperOptCache::reset_step)
        .def("step", &TemperOptCache::step, py::arg("v"));

    // ── PISCache (StackBase strategy for tempering optimizer) ──────────

    py::class_<PISCache>(m, "PISCache")
        .def(py::init([](const py::list& gens_py, const py::list& trans_py,
                         int kg, int L) {
            std::vector<Generator*> gens;
            for (auto item : gens_py) gens.push_back(item.cast<Generator*>());
            std::vector<std::vector<Transformation*>> trans;
            for (auto ct : trans_py) {
                std::vector<Transformation*> chain;
                for (auto t : ct) chain.push_back(t.cast<Transformation*>());
                trans.push_back(chain);
            }
            return PISCache(gens, trans, kg, L);
        }), py::arg("gens"), py::arg("trans"), py::arg("kg"), py::arg("L"))
        // Phase 1 single-Generator& constructor.
        .def(py::init([](const Generator& gen, int kg, int L) {
            return PISCache(gen, kg, L);
        }), py::arg("gen"), py::arg("kg"), py::arg("L"))
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
        std::vector<Generator*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<Generator*>());
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

    // ── tausworthe_random_poly: sample an admissible polynomial ──────

    m.def("tausworthe_random_poly",
          [](int k, int nb_terms, bool quicktaus, int L, int s) {
        return TauswortheGen::random_poly(k, nb_terms, quicktaus, L, s);
    }, py::arg("k"), py::arg("nb_terms"), py::arg("quicktaus"),
       py::arg("L"), py::arg("s") = 0,
       "Sample a random admissible TauswortheGen polynomial.  Returns the "
       "sorted exponent list [0, q_1, ..., q_{t-2}, k].  Throws if the "
       "(k, nb_terms, quicktaus, L, s) combination is inadmissible.");

    // ── random_param: dispatch a non-generic rand_type to its family ──
    //
    // parametric.py owns the generic samplers (bitmask, range,
    // poly_exponents, bitmask_vec) and hands anything else off here.
    // Currently only TauswortheGen registers samplers; adding a new
    // family means one more if-branch below plus a static
    // `generate_random` method on that family.
    m.def("random_param",
          [](const std::string& rand_type,
             const std::string& rand_args,
             const py::dict& params_dict,
             int L) -> py::tuple {
        Params p = dict_to_params(params_dict);
        RandomParamResult r;
        if (rand_type == "tausworthe_s"
            || rand_type == "tausworthe_poly") {
            r = TauswortheGen::generate_random(
                rand_type, rand_args, p, L);
        } else {
            throw std::invalid_argument(
                "random_param: unsupported rand_type '"
                + rand_type + "'");
        }
        py::object value = r.is_vec
            ? py::cast(r.vec_val)
            : py::object(py::int_(r.int_val));
        py::dict side;
        for (const auto& kv : r.side_ints)
            side[py::str(kv.first)] = py::int_(kv.second);
        return py::make_tuple(value, side);
    }, py::arg("rand_type"), py::arg("rand_args"),
       py::arg("params"), py::arg("L"),
       "Sample a random value for a family-specific rand_type.  "
       "Returns (value, side_effects) where side_effects is a dict of "
       "extra params to splice into the caller's bag (e.g. the `s` "
       "paired with a freshly-sampled TauswortheGen poly).");

    // ── Exhaustive-search enumerator ──────────────────────────────────

    py::class_<GenEnumerator, std::shared_ptr<GenEnumerator>>(m, "GenEnumerator")
        .def("size", [](const GenEnumerator& e) {
            // Lift decimal-string count to a Python int of arbitrary precision.
            return py::module_::import("builtins").attr("int")(e.size_dec());
        })
        .def("at", [](const GenEnumerator& e, const py::object& idx) {
            return params_to_dict(e.at(py::str(idx).cast<std::string>()));
        }, py::arg("idx"))
        .def("axes", [](const GenEnumerator& e) {
            py::list out;
            auto to_int = py::module_::import("builtins").attr("int");
            for (const auto& a : e.axes()) {
                py::dict d;
                d["name"]     = a.name;
                d["size"]     = to_int(a.size_dec);
                d["describe"] = a.describe;
                out.append(d);
            }
            return out;
        });

    m.def("make_gen_enumerator",
          [](const std::string& family, const py::dict& d, int L) -> py::object {
        auto e = make_gen_enumerator(family, dict_to_params(d), L);
        if (!e) return py::none();
        std::shared_ptr<GenEnumerator> shared(std::move(e));
        return py::cast(shared);
    }, py::arg("family"), py::arg("resolved"), py::arg("L"),
       "Build the exhaustive-search enumerator for a family.  Returns "
       "None when the family has no registered enumerator.  Throws "
       "std::invalid_argument('needs_<reason>') when the resolved "
       "inputs are insufficient.");

    m.def("family_is_enumerable",
          [](const std::string& family) {
        return family_is_enumerable(family);
    }, py::arg("family"),
       "True iff the family has a registered exhaustive enumerator.");

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

    // ── Primitivity (Phase 2.4): full-period testing in C++ ────────────
    //
    // Replaces packages/regpoly/src/regpoly/search/primitivity.py. The
    // Cunningham-style factor table is embedded via
    // src/algebra/primitive_factors_data.cpp.

    m.def("is_mersenne_prime_exponent", &is_mersenne_prime_exponent,
          py::arg("k"),
          "True iff 2^k - 1 is a Mersenne prime (single-factor primitivity).");

    m.def("get_primitive_factors_for_k",
          [](int k) -> py::object {
        auto facs = get_primitive_factors_for_k(k);
        if (!facs) return py::none();
        return py::cast(*facs);
    }, py::arg("k"),
       "Sorted list of prime factors of 2^k - 1 as decimal strings, or "
       "None when the factorisation is unavailable.");

    m.def("is_full_period",
          [](const Generator& gen) -> bool { return is_full_period(gen); },
          py::arg("gen"),
          "True iff the generator's characteristic polynomial is "
          "primitive (period 2^k - 1).");

    // ── Search drivers (Phase 2.4) ─────────────────────────────────────
    //
    // Free functions whose loops own the iteration, randomization,
    // factory call and primitivity check. Python registers callbacks
    // for the per-hit and per-progress events; the rest of the
    // orchestration (YAML, .partial recovery, dedup, file merge) stays
    // in Python for now.

    py::class_<SearchProgress>(m, "SearchProgress")
        .def_readonly("tries", &SearchProgress::tries)
        .def_readonly("elapsed_seconds", &SearchProgress::elapsed_seconds);

    py::class_<TestedGenerator>(m, "TestedGenerator")
        .def_readonly("family", &TestedGenerator::family)
        .def_readonly("k", &TestedGenerator::k)
        .def_readonly("L", &TestedGenerator::L)
        .def_readonly("tries_at_hit", &TestedGenerator::tries_at_hit)
        .def_property_readonly("params", [](const TestedGenerator& tg) {
            return params_to_dict(tg.params);
        });

    py::class_<PrimitiveSearchConfig>(m, "PrimitiveSearchConfig")
        .def(py::init<>())
        .def_readwrite("family", &PrimitiveSearchConfig::family)
        .def_readwrite("L", &PrimitiveSearchConfig::L)
        .def_property("structural_params",
            [](const PrimitiveSearchConfig& c) {
                return params_to_dict(c.structural_params);
            },
            [](PrimitiveSearchConfig& c, const py::dict& d) {
                c.structural_params = dict_to_params(d);
            })
        .def_property("fixed_params",
            [](const PrimitiveSearchConfig& c) {
                return params_to_dict(c.fixed_params);
            },
            [](PrimitiveSearchConfig& c, const py::dict& d) {
                c.fixed_params = dict_to_params(d);
            })
        .def_readwrite("max_tries", &PrimitiveSearchConfig::max_tries)
        .def_readwrite("max_seconds", &PrimitiveSearchConfig::max_seconds)
        .def_readwrite("progress_interval",
                       &PrimitiveSearchConfig::progress_interval)
        .def_readwrite("random_seed", &PrimitiveSearchConfig::random_seed);

    m.def("run_primitive_search",
          [](const PrimitiveSearchConfig& cfg,
             const py::function& on_hit,
             const py::function& on_progress) -> int64_t {
        OnHitFn hit = [&on_hit](const TestedGenerator& tg) {
            on_hit(tg);
        };
        OnProgressFn prog = [&on_progress](const SearchProgress& sp) {
            on_progress(sp);
        };
        return run_primitive_search(cfg, hit, prog);
    }, py::arg("config"), py::arg("on_hit"), py::arg("on_progress"),
       "Run the full-period search loop in C++. Invokes on_hit(tg) "
       "for each hit and on_progress(sp) every progress_interval "
       "tries plus once at completion. Returns the total tries "
       "executed.");

    // ── TemperingOptimizerDriver (Phase 2.4d) ─────────────────────────
    //
    // Replaces the recursive optimize(v) loop in
    // regpoly.search.tempering_optimizer.TemperingOptimizer.run_once.
    // The Python wrapper still computes safe_masks (a small structural
    // computation that depends on `mu`/width); the driver consumes
    // them and owns the hot perturbation + cache.step inner loop.

    py::class_<TemperingOptimizerConfig>(m, "TemperingOptimizerConfig")
        .def(py::init<>())
        .def_readwrite("max_essais",  &TemperingOptimizerConfig::max_essais)
        .def_readwrite("delta",       &TemperingOptimizerConfig::delta)
        .def_readwrite("mse",         &TemperingOptimizerConfig::mse)
        .def_readwrite("n_restarts",  &TemperingOptimizerConfig::n_restarts)
        .def_readwrite("random_seed", &TemperingOptimizerConfig::random_seed);

    py::class_<TemperingOptResult>(m, "TemperingOptResult")
        .def_readonly("se",               &TemperingOptResult::se)
        .def_readonly("gaps",             &TemperingOptResult::gaps)
        .def_readonly("elapsed_seconds",  &TemperingOptResult::elapsed_seconds)
        .def_readonly("essais",           &TemperingOptResult::essais);

    // ParamLocators are passed in as a list of (cpp_trans, name, width,
    // current_value) tuples. The driver mutates current_value in-place
    // and writes the value back via cpp_trans.update({name: value}); on
    // return, the locator list reflects the best-found values, and the
    // Python caller can sync those back into its own _params dicts.
    auto build_locators = [](const py::list& tuples)
        -> std::vector<TemperParamLocator> {
        std::vector<TemperParamLocator> out;
        out.reserve(tuples.size());
        for (auto h : tuples) {
            auto t = h.cast<py::tuple>();
            TemperParamLocator loc;
            loc.trans = t[0].cast<Transformation*>();
            loc.param_name = t[1].cast<std::string>();
            loc.width = t[2].cast<int>();
            loc.current_value = t[3].cast<int64_t>();
            out.push_back(loc);
        }
        return out;
    };

    auto export_locators = [](const std::vector<TemperParamLocator>& locs)
        -> py::list {
        py::list out;
        for (const auto& loc : locs)
            out.append(py::int_(loc.current_value));
        return out;
    };

    m.def("run_tempering_optimizer_once",
          [build_locators, export_locators](
              const TemperingOptimizerConfig& cfg,
              TemperOptCache& cache,
              const py::list& param_tuples,
              const std::vector<std::vector<uint64_t>>& safe_masks) -> py::tuple {
        auto params = build_locators(param_tuples);
        auto result = run_tempering_optimizer_once(
            cfg, cache, params, safe_masks);
        return py::make_tuple(result, export_locators(params));
    }, py::arg("config"), py::arg("cache"), py::arg("params"),
       py::arg("safe_masks"),
       "Single recursive optimization pass. Returns "
       "(TemperingOptResult, [final_value_per_locator]).");

    // ── SeekDriver (Phase 2.4b) ────────────────────────────────────────
    //
    // Drives the equidistribution / collision-free / tuplets search
    // loop in C++. The Python Seek.run() shrinks to: build a C++
    // Combination from the Python Combination, build a list of
    // SeekTestSpec, register on_prep / on_iter / on_progress
    // callbacks for re-randomization / persistence / progress, call
    // run_seek_search.

    py::enum_<SeekTestKind>(m, "SeekTestKind")
        // Canonical (post-R3/R4): test type only — method is in
        // SeekTestSpec.method_name.
        .value("Equidistribution",  SeekTestKind::Equidistribution)
        .value("CollisionFree",     SeekTestKind::CollisionFree)
        .value("Tuplets",           SeekTestKind::Tuplets)
        // Deprecated aliases. Translate internally to
        // (Equidistribution, method_name=...) inside run_seek_search.
        // Kept for backward compatibility with existing Python callers.
        .value("EquidistributionMatricial",
               SeekTestKind::EquidistributionMatricial)
        .value("EquidistributionLattice",
               SeekTestKind::EquidistributionLattice)
        .value("EquidistributionHarase",
               SeekTestKind::EquidistributionHarase)
        .value("EquidistributionNotPrimitive",
               SeekTestKind::EquidistributionNotPrimitive)
        .value("EquidistributionSimdNotPrimitive",
               SeekTestKind::EquidistributionSimdNotPrimitive)
        .value("EquidistributionNothing",
               SeekTestKind::EquidistributionNothing);

    // Exposed so the Python wrapper and tests can share the canonical
    // vocabulary with C++ instead of maintaining a parallel string map.
    m.def("equidistribution_method_names",
          []() { return MethodRegistry::names(); },
          "List the equidistribution method names known to "
          "MethodRegistry. The Python YAML parser consults this list "
          "rather than maintaining its own parallel enum.");
    m.def("has_equidistribution_method",
          [](const std::string& name) { return MethodRegistry::has(name); },
          py::arg("name"));

    py::class_<SeekTestSpec>(m, "SeekTestSpec")
        .def(py::init<>())
        .def_readwrite("kind",            &SeekTestSpec::kind)
        .def_readwrite("method_name",     &SeekTestSpec::method_name)
        .def_readwrite("eq_L_max_test",   &SeekTestSpec::eq_L_max_test)
        .def_readwrite("eq_delta",        &SeekTestSpec::eq_delta)
        .def_readwrite("eq_mse",          &SeekTestSpec::eq_mse)
        .def_readwrite("tup_d",           &SeekTestSpec::tup_d)
        .def_readwrite("tup_h",           &SeekTestSpec::tup_h)
        .def_readwrite("tup_threshold",   &SeekTestSpec::tup_threshold)
        .def_readwrite("tup_testtype",    &SeekTestSpec::tup_testtype);

    py::class_<SeekIterResult>(m, "SeekIterResult")
        .def_readonly("selected",          &SeekIterResult::selected)
        .def_readonly("me_ran",            &SeekIterResult::me_ran)
        .def_readonly("me_verified",       &SeekIterResult::me_verified)
        .def_readonly("me_is_me",          &SeekIterResult::me_is_me)
        .def_readonly("me_se",             &SeekIterResult::me_se)
        .def_readonly("me_test_L",         &SeekIterResult::me_test_L)
        .def_readonly("me_ecart",          &SeekIterResult::me_ecart)
        .def_readonly("tup_ran",           &SeekIterResult::tup_ran)
        .def_readonly("tup_verified",      &SeekIterResult::tup_verified)
        .def_readonly("tup_is_ok",         &SeekIterResult::tup_is_ok)
        .def_readonly("tup_firstpart_max", &SeekIterResult::tup_firstpart_max)
        .def_readonly("tup_firstpart_sum", &SeekIterResult::tup_firstpart_sum)
        .def_readonly("tup_secondpart_max",&SeekIterResult::tup_secondpart_max)
        .def_readonly("tup_secondpart_sum",&SeekIterResult::tup_secondpart_sum)
        .def_readonly("cf_ran",            &SeekIterResult::cf_ran)
        .def_readonly("cf_verified",       &SeekIterResult::cf_verified)
        .def_readonly("cf_secf",           &SeekIterResult::cf_secf)
        .def_readonly("cf_ecart_cf",       &SeekIterResult::cf_ecart_cf);

    py::class_<SeekProgress>(m, "SeekProgress")
        .def_readonly("nbgen",          &SeekProgress::nbgen)
        .def_readonly("nb_select",      &SeekProgress::nb_select)
        .def_readonly("nb_me",          &SeekProgress::nb_me)
        .def_readonly("elapsed_seconds",&SeekProgress::elapsed_seconds);

    py::class_<SeekResult>(m, "SeekResult")
        .def_readonly("nbgen",          &SeekResult::nbgen)
        .def_readonly("nb_select",      &SeekResult::nb_select)
        .def_readonly("nb_me",          &SeekResult::nb_me)
        .def_readonly("elapsed_seconds",&SeekResult::elapsed_seconds);

    m.def("run_seek_search",
          [](Combination& comb,
             const std::vector<SeekTestSpec>& tests,
             int nbtries,
             int progress_interval,
             const py::object& on_prep,
             const py::object& on_iter,
             const py::object& on_progress) -> SeekResult {
        SeekOnPrepFn prep_fn = nullptr;
        if (!on_prep.is_none()) {
            prep_fn = [&on_prep](Combination& c, bool is_retry) {
                on_prep(py::cast(&c, py::return_value_policy::reference),
                        is_retry);
            };
        }
        SeekOnIterFn iter_fn = nullptr;
        if (!on_iter.is_none()) {
            iter_fn = [&on_iter](Combination& c, const SeekIterResult& r) {
                on_iter(py::cast(&c, py::return_value_policy::reference), r);
            };
        }
        SeekOnProgressFn prog_fn = nullptr;
        if (!on_progress.is_none()) {
            prog_fn = [&on_progress](const SeekProgress& p) {
                on_progress(p);
            };
        }
        return run_seek_search(comb, tests, nbtries, progress_interval,
                               prep_fn, iter_fn, prog_fn);
    }, py::arg("combination"), py::arg("tests"), py::arg("nbtries"),
       py::arg("progress_interval"),
       py::arg("on_prep")     = py::none(),
       py::arg("on_iter")     = py::none(),
       py::arg("on_progress") = py::none(),
       "Run the seek search loop in C++. Iterates the combination, "
       "runs the configured tests in order, and emits callbacks for "
       "selections + periodic progress. Returns a SeekResult.");

    // ── TemperingSearchDriver (Phase 2.4c) ─────────────────────────────
    //
    // Drives the per-combo / per-try search loop in C++. Per-try work
    // (re-randomize tempering params, optimize, run test) lives in the
    // Python on_try callback because randomize_params lives on Python
    // Transformations.

    py::class_<TemperingSearchConfig>(m, "TemperingSearchConfig")
        .def(py::init<>())
        .def_readwrite("nb_tries",          &TemperingSearchConfig::nb_tries)
        .def_readwrite("progress_interval", &TemperingSearchConfig::progress_interval);

    py::class_<TemperingTryOutcome>(m, "TemperingTryOutcome")
        .def(py::init<>())
        .def_readwrite("got_result", &TemperingTryOutcome::got_result)
        .def_readwrite("score",      &TemperingTryOutcome::score);

    py::class_<TemperingSearchResult>(m, "TemperingSearchResult")
        .def_readonly("nbgen",           &TemperingSearchResult::nbgen)
        .def_readonly("nb_with_result",  &TemperingSearchResult::nb_with_result)
        .def_readonly("elapsed_seconds", &TemperingSearchResult::elapsed_seconds);

    m.def("run_tempering_search",
          [](Combination& comb,
             const TemperingSearchConfig& cfg,
             const py::object& on_combo_start,
             const py::object& on_try,
             const py::object& on_combo_done,
             const py::object& on_progress) -> TemperingSearchResult {
        if (on_try.is_none()) {
            throw std::invalid_argument(
                "run_tempering_search: on_try callback is required");
        }
        TempSearchOnComboStartFn start_fn = nullptr;
        if (!on_combo_start.is_none()) {
            start_fn = [&on_combo_start](Combination& c, int idx) {
                on_combo_start(
                    py::cast(&c, py::return_value_policy::reference), idx);
            };
        }
        TempSearchOnTryFn try_fn =
            [&on_try](Combination& c, int combo_idx, int try_idx,
                      bool is_first) -> TemperingTryOutcome {
                py::object r = on_try(
                    py::cast(&c, py::return_value_policy::reference),
                    combo_idx, try_idx, is_first);
                return r.cast<TemperingTryOutcome>();
            };
        TempSearchOnComboDoneFn done_fn = nullptr;
        if (!on_combo_done.is_none()) {
            done_fn = [&on_combo_done](Combination& c, int idx,
                                       int best_score, int best_try) {
                on_combo_done(
                    py::cast(&c, py::return_value_policy::reference),
                    idx, best_score, best_try);
            };
        }
        TempSearchOnProgressFn prog_fn = nullptr;
        if (!on_progress.is_none()) {
            prog_fn = [&on_progress](const SearchProgress& p) {
                on_progress(p);
            };
        }
        return run_tempering_search(
            comb, cfg, start_fn, try_fn, done_fn, prog_fn);
    }, py::arg("combination"), py::arg("config"),
       py::arg("on_combo_start") = py::none(),
       py::arg("on_try"),
       py::arg("on_combo_done")  = py::none(),
       py::arg("on_progress")    = py::none(),
       "Run the tempering search loop in C++. Per combo, fires on_try "
       "nb_tries times (Python re-randomizes + optimizes + tests + "
       "returns a TemperingTryOutcome). Tracks the best score per combo, "
       "fires on_combo_done with the best result, advances the "
       "Combination, and emits on_progress every progress_interval tries. "
       "Returns a TemperingSearchResult.");

    m.def("run_tempering_optimizer_minimize",
          [build_locators, export_locators](
              const TemperingOptimizerConfig& cfg,
              TemperOptCache& cache,
              const py::list& param_tuples,
              const std::vector<std::vector<uint64_t>>& safe_masks) -> py::tuple {
        auto params = build_locators(param_tuples);
        auto result = run_tempering_optimizer_minimize(
            cfg, cache, params, safe_masks);
        return py::make_tuple(result, export_locators(params));
    }, py::arg("config"), py::arg("cache"), py::arg("params"),
       py::arg("safe_masks"),
       "Iterative delta-tightening loop (n_restarts > 1). Returns "
       "(TemperingOptResult, [final_value_per_locator]).");

    auto specs_to_list = [](const std::vector<ParamSpec>& specs) -> py::list {
        py::list result;
        for (auto& s : specs) {
            py::dict d;
            d["name"]        = s.name;
            d["type"]        = s.type;
            d["structural"]  = s.structural;
            d["has_default"] = s.has_default;
            d["default"]     = s.default_val;
            d["rand_type"]   = s.rand_type;
            d["rand_args"]   = s.rand_args;
            d["optimizable"] = s.optimizable;
            result.append(d);
        }
        return result;
    };

    m.def("get_gen_param_specs",
          [&specs_to_list](const std::string& family) -> py::list {
        return specs_to_list(get_gen_param_specs(family));
    }, py::arg("family"));

    m.def("get_trans_param_specs",
          [&specs_to_list](const std::string& type) -> py::list {
        return specs_to_list(get_trans_param_specs(type));
    }, py::arg("type"));

    // ── Catalog (Phase 3.2) ────────────────────────────────────────────
    //
    // The C++ Catalog reads docs/library/*.yaml and exposes Paper /
    // CatalogGenerator / Author records. The Python regpoly.library
    // package is a thin shim that re-exports these types and adds
    // Path-typed source_path / dict-typed params.

    using namespace regpoly_catalog;

    // ParamValue → native Python value (int / str / bool / list[int]).
    auto pv_to_py = [](const ParamValue& v) -> py::object {
        switch (v.kind) {
            case ParamKind::Int:    return py::int_(v.int_val);
            case ParamKind::String: return py::str(v.string_val);
            case ParamKind::Bool:   return py::bool_(v.bool_val);
            case ParamKind::IntList: return py::cast(v.int_list_val);
        }
        return py::none();
    };
    // ParamMap → dict.
    auto pmap_to_py = [pv_to_py](const ParamMap& m) -> py::dict {
        py::dict d;
        for (const auto& kv : m) d[py::str(kv.first)] = pv_to_py(kv.second);
        return d;
    };
    // TemperingStep → dict {"type": ..., other params...}.
    auto step_to_py = [pmap_to_py](const TemperingStep& t) -> py::dict {
        py::dict d = pmap_to_py(t.params);
        d["type"] = py::str(t.type);
        return d;
    };
    // Component → dict matching Python's _normalize_components shape.
    auto comp_to_py = [pmap_to_py, step_to_py](
        const regpoly_catalog::Component& c) -> py::dict {
        py::dict d;
        d["family"] = py::str(c.family);
        d["L"] = py::int_(c.L);
        d["params"] = pmap_to_py(c.params);
        py::list temper;
        for (const auto& s : c.tempering) temper.append(step_to_py(s));
        d["tempering"] = temper;
        return d;
    };

    py::class_<Author>(m, "Author")
        .def(py::init<>())
        .def_readwrite("family", &Author::family)
        .def_readwrite("given",  &Author::given)
        .def("display", &Author::display)
        .def("short",   &Author::short_name);

    py::class_<CatalogGenerator>(m, "CatalogGenerator")
        .def_readonly("id",        &CatalogGenerator::id)
        .def_readonly("display",   &CatalogGenerator::display)
        .def_readonly("family",    &CatalogGenerator::family)
        .def_readonly("target",    &CatalogGenerator::target)
        .def_readonly("combined",  &CatalogGenerator::combined)
        .def_readonly("Lmax",      &CatalogGenerator::Lmax)
        .def_readonly("notes_md",  &CatalogGenerator::notes_md)
        .def_readonly("starred",   &CatalogGenerator::starred)
        .def_readonly("errors",    &CatalogGenerator::errors)
        .def_property_readonly("valid", &CatalogGenerator::valid)
        .def_property_readonly("components",
            [comp_to_py](const CatalogGenerator& g) -> py::list {
                py::list out;
                for (const auto& c : g.components) out.append(comp_to_py(c));
                return out;
            });

    py::class_<Paper>(m, "Paper")
        .def_readonly("id",          &Paper::id)
        .def_readonly("authors",     &Paper::authors)
        .def_readonly("year",        &Paper::year)
        .def_readonly("title",       &Paper::title)
        .def_readonly("venue",       &Paper::venue)
        .def_readonly("volume",      &Paper::volume)
        .def_readonly("issue",       &Paper::issue)
        .def_readonly("pages",       &Paper::pages)
        .def_readonly("doi",         &Paper::doi)
        .def_readonly("pdf",         &Paper::pdf)
        .def_readonly("bibkey",      &Paper::bibkey)
        .def_readonly("abstract_md", &Paper::abstract_md)
        .def_readonly("notes_md",    &Paper::notes_md)
        .def_readonly("tags",        &Paper::tags)
        .def_readonly("starred",     &Paper::starred)
        .def_readonly("deferred",    &Paper::deferred)
        .def_readonly("generators",  &Paper::generators)
        .def_readonly("source_path", &Paper::source_path)
        .def_readonly("source_mtime",&Paper::source_mtime)
        .def_readonly("errors",      &Paper::errors)
        .def_property_readonly("valid",            &Paper::valid)
        .def("author_list_short",  &Paper::author_list_short)
        .def("display",            &Paper::display)
        .def("acmtrans_citation",  &Paper::acmtrans_citation);

    py::class_<Catalog>(m, "Catalog")
        .def(py::init<std::string>(), py::arg("library_dir"))
        .def("load", &Catalog::load)
        .def("reload_if_stale", &Catalog::reload_if_stale)
        .def("library_dir", &Catalog::library_dir)
        .def("papers",
             [](const Catalog& c, py::object starred, py::object tag,
                bool include_invalid) {
                 Catalog::PapersFilter f;
                 if (!starred.is_none()) f.starred = starred.cast<bool>();
                 if (!tag.is_none()) f.tag = tag.cast<std::string>();
                 f.include_invalid = include_invalid;
                 return c.papers(f);
             },
             py::arg("starred") = py::none(),
             py::arg("tag") = py::none(),
             py::arg("include_invalid") = false)
        .def("paper",
             [](const Catalog& c, const std::string& id) -> py::object {
                 auto p = c.paper(id);
                 if (!p.has_value()) return py::none();
                 return py::cast(*p);
             }, py::arg("paper_id"))
        .def("generator",
             [](const Catalog& c, const std::string& gid) -> py::object {
                 auto loc = c.generator(gid);
                 if (!loc.has_value()) return py::none();
                 return py::make_tuple(loc->first, loc->second);
             }, py::arg("gen_id"))
        .def("all_generators",
             [](const Catalog& c, py::object family) {
                 std::optional<std::string> f;
                 if (!family.is_none()) f = family.cast<std::string>();
                 auto results = c.all_generators(f);
                 py::list out;
                 for (auto& pg : results) {
                     out.append(py::make_tuple(pg.first, pg.second));
                 }
                 return out;
             }, py::arg("family") = py::none());

    m.def("config_hash",
          [](const std::string& family, const py::dict& params,
             const py::list& tempering) -> std::string {
              auto py_to_pv = [](const py::handle& v) -> ParamValue {
                  if (py::isinstance<py::bool_>(v))
                      return ParamValue::make_bool(v.cast<bool>());
                  if (py::isinstance<py::int_>(v))
                      return ParamValue::make_int(v.cast<int64_t>());
                  if (py::isinstance<py::str>(v))
                      return ParamValue::make_string(v.cast<std::string>());
                  if (py::isinstance<py::list>(v)) {
                      try {
                          return ParamValue::make_int_list(
                              v.cast<std::vector<int64_t>>());
                      } catch (...) {
                          return ParamValue::make_string(py::str(v).cast<std::string>());
                      }
                  }
                  return ParamValue::make_string(py::str(v).cast<std::string>());
              };
              ParamMap pm;
              for (auto item : params) {
                  pm.emplace(item.first.cast<std::string>(),
                             py_to_pv(item.second));
              }
              std::vector<TemperingStep> steps;
              for (const auto& tn : tempering) {
                  py::dict td = tn.cast<py::dict>();
                  TemperingStep st;
                  if (td.contains("type")) {
                      st.type = td["type"].cast<std::string>();
                  }
                  for (auto kv : td) {
                      auto k = kv.first.cast<std::string>();
                      if (k == "type") continue;
                      st.params.emplace(k, py_to_pv(kv.second));
                  }
                  steps.push_back(std::move(st));
              }
              return regpoly_catalog::config_hash(family, pm, steps);
          },
          py::arg("family"), py::arg("params"), py::arg("tempering"),
          "Stable short hash of one component config.");

    // ── Legacy reader (Phase 3.3) ──────────────────────────────────────
    //
    // The C++ regpoly_legacy::read_*_specs returns Params records. We
    // surface them to Python as (family, params_dict) tuples so the
    // Python wrapper can call Generator.create / Transformation.create
    // and preserve its _params dict for downstream randomisation.

    m.def("legacy_read_generator_specs",
          [](const std::string& filename, int L) -> py::list {
              auto specs = regpoly_legacy::read_generator_specs(filename, L);
              py::list out;
              for (const auto& s : specs) {
                  out.append(py::make_tuple(
                      py::str(s.family), params_to_dict(s.params)));
              }
              return out;
          }, py::arg("filename"), py::arg("L"),
          "Parse a legacy .dat generator file and return a list of "
          "(family, params_dict) tuples. Used by the Python "
          "regpoly.io.legacy_reader shim to round-trip via "
          "Generator.create() so the Python wrapper's _params dict "
          "stays populated for randomisation paths.");

    m.def("legacy_read_transformation_specs",
          [](const std::string& filename) -> py::tuple {
              auto r = regpoly_legacy::read_transformation_specs(filename);
              py::list specs;
              for (const auto& s : r.specs) {
                  specs.append(py::make_tuple(
                      py::str(s.trans_type), params_to_dict(s.params)));
              }
              return py::make_tuple(specs, r.mk_opt);
          }, py::arg("filename"),
          "Parse a legacy .dat transformations file. Returns a "
          "(specs_list, mk_opt) tuple; specs_list is a list of "
          "(type, params_dict) tuples.");
}
