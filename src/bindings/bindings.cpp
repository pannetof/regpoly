// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "bitvect.h"
#include "bm.h"
#include "combo_enumerator.h"
#include "combined_f2_linear_source.h"
#include "equidistribution_method.h"
#include "equidistribution_runner.h"
#include "tvalue_runner.h"
#include "digital_net.h"
#include "sobol.h"
#include "niederreiter_f2.h"
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
#include "temper_optimizer.h"
#include "tempering_optimizer.h"
#include "random_samplers.h"
#include "tausworthe.h"
#include "param_spec.h"
#include "tuplets_runner.h"
#include "well.h"

#include <NTL/GF2X.h>
#include <cctype>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>
#include <random>

using namespace regpoly::core;
using namespace regpoly::internal;
using namespace regpoly::random;


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
// py::dict -> ParamBag conversion
// ═══════════════════════════════════════════════════════════════════════════

// Convert a Python dict-of-scalars into a StructEntry (one row of a
// StructMap). Every value must be int / bool / str / hex-string. Used
// by the dict_to_params struct_map branch below.
static StructEntry py_dict_to_struct_entry(const py::handle& val,
                                            const std::string& slot_key) {
    if (!py::isinstance<py::dict>(val))
        throw std::runtime_error(
            "matrices['" + slot_key + "'] must be a dict (got "
            + std::string(py::str(val.get_type().attr("__name__"))) + ")");
    StructEntry e;
    auto d = val.cast<py::dict>();
    for (auto kv : d) {
        std::string k = kv.first.cast<std::string>();
        py::handle v = kv.second;
        // Order matters: bool is a subclass of int in Python.
        if (py::isinstance<py::bool_>(v)) {
            e[k] = v.cast<bool>();
        } else if (py::isinstance<py::int_>(v)) {
            // Use unsigned representation for non-negative values that
            // fit in 32 bits unsigned but not int64_t. Default to int64.
            try {
                e[k] = v.cast<int64_t>();
            } catch (...) {
                e[k] = v.cast<uint64_t>();
            }
        } else if (py::isinstance<py::str>(v)) {
            e[k] = v.cast<std::string>();
        } else {
            throw std::runtime_error(
                "matrices['" + slot_key + "']." + k
                + " must be int / bool / str");
        }
    }
    return e;
}

static ParamBag dict_to_params(const py::dict& d) {
    ParamBag p;
    for (auto item : d) {
        std::string key = item.first.cast<std::string>();
        py::handle val = item.second;

        // Reject the legacy WELL flat triple at the binding layer with
        // a clear migration pointer. Some other generators may legitimately
        // use these names in a non-WELL context (none today), so the check
        // is by literal key only.
        if (key == "mat_types" || key == "mat_pi" || key == "mat_pu") {
            throw std::runtime_error(
                "WELLGen: '" + key + "' is no longer accepted. "
                "Use the structured 'matrices' map keyed by T0..T7. "
                "See docs/generators/WELLGen.md.");
        }

        // Structured dict-of-dicts (e.g. WELL `matrices`). Recognise by
        // shape: a top-level py::dict whose keys are str and whose
        // values are themselves py::dict of scalars.
        if (py::isinstance<py::dict>(val)) {
            StructMap m;
            auto outer = val.cast<py::dict>();
            for (auto sub : outer) {
                std::string slot = sub.first.cast<std::string>();
                m[slot] = py_dict_to_struct_entry(sub.second, slot);
            }
            p.set_struct_map(key, std::move(m));
            continue;
        }

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
// ParamBag -> py::dict conversion (mirror of dict_to_params)
// ═══════════════════════════════════════════════════════════════════════════

static py::object scalar_to_py(const ParamScalar& s) {
    if (auto pi = std::get_if<int64_t>(&s)) return py::int_(*pi);
    if (auto pu = std::get_if<uint64_t>(&s)) return py::int_(*pu);
    if (auto ps = std::get_if<std::string>(&s)) return py::str(*ps);
    if (auto pb = std::get_if<bool>(&s)) return py::bool_(*pb);
    return py::none();
}

static py::dict params_to_dict(const ParamBag& p) {
    py::dict d;
    for (const auto& kv : p.ints())      d[py::str(kv.first)] = py::int_(kv.second);
    for (const auto& kv : p.bools())     d[py::str(kv.first)] = py::bool_(kv.second);
    for (const auto& kv : p.strings())   d[py::str(kv.first)] = py::str(kv.second);
    for (const auto& kv : p.int_vecs())  d[py::str(kv.first)] = py::cast(kv.second);
    for (const auto& kv : p.uint_vecs()) d[py::str(kv.first)] = py::cast(kv.second);
    for (const auto& kv : p.struct_maps()) {
        py::dict outer;
        for (const auto& slot : kv.second) {
            py::dict inner;
            for (const auto& arg : slot.second) {
                inner[py::str(arg.first)] = scalar_to_py(arg.second);
            }
            outer[py::str(slot.first)] = inner;
        }
        d[py::str(kv.first)] = outer;
    }
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

    // ── ITestable / F2LinearSource (Phase 3 polymorphic interface) ─────
    //
    // Pybind11 needs both ancestor classes bound to wire the inheritance
    // chain so `Recurrence` subclasses dispatch correctly when callers
    // pass them where `const F2LinearSource&` is expected.

    py::class_<ITestable>(m, "ITestable");

    // F2LinearSource exposes the shared surface (k/L/name/init/next/
    // get_output/copy/default_test_method) so both `Recurrence`
    // subclasses (PRNGs) and `DigitalNet` subclasses inherit them via
    // Python's MRO. Recurrence-specific methods (char_poly, state,
    // transition_matrix, is_full_period, components) stay on Recurrence.
    py::class_<F2LinearSource, ITestable>(m, "F2LinearSource")
        .def("name", &F2LinearSource::name)
        .def("display_str", &F2LinearSource::display_str)
        .def("k", &F2LinearSource::k)
        .def("L", &F2LinearSource::L)
        .def("init", &F2LinearSource::init)
        .def("next", &F2LinearSource::next)
        .def("get_output", &F2LinearSource::get_output)
        .def("copy", [](const F2LinearSource& s) { return s.copy(); })
        // Releases the GIL: the call can take several seconds for
        // non-Mersenne K (Pollard's rho on the cofactor of 2^K-1).
        // Without the release the whole FastAPI event loop stalls
        // for the duration — the library detail page would render
        // as "nothing at all" because every other request queues
        // behind this one.
        .def("default_test_method", &F2LinearSource::default_test_method,
             py::arg("test_type"),
             py::call_guard<py::gil_scoped_release>());

    // ── Recurrence ──────────────────────────────────────────────────────

    py::class_<Recurrence, F2LinearSource, std::unique_ptr<Recurrence>>(m, "Recurrence")
        .def("char_poly", &Recurrence::char_poly)
        .def("is_full_period", &Recurrence::is_full_period)
        .def("transition_matrix", &Recurrence::transition_matrix)
        .def("state", [](const Recurrence& g) -> BitVect { return g.state().copy(); })
        // Polymorphic unpack — a primitive Recurrence returns `{this}`;
        // CombinedF2LinearSource overrides to expose its J components.
        // Returned pointers alias the Recurrence's internal state and
        // are valid for its lifetime.
        .def("components", &Recurrence::components,
             py::return_value_policy::reference_internal);

    // ── CombinedF2LinearSource ───────────────────────────────────────────────
    //
    // Composition: J Recurrence components XOR'd together. Inherits
    // `ITestable` directly (NOT `Recurrence`) — by design, a combined
    // source is composition, not a single recurrence. Lattice / matricial
    // kernels walk the J components via `components()` rather than
    // pretending the wrapper has a combined state.
    //
    // Construction factories transfer ownership: each py::list of
    // components is consumed by clone_recurrence().

    py::class_<CombinedF2LinearSource, ITestable,
               std::unique_ptr<CombinedF2LinearSource>>(
        m, "CombinedF2LinearSource")
        .def(py::init([](py::list py_components, int Lmax) {
                 std::vector<std::unique_ptr<F2LinearSource>> comps;
                 comps.reserve(py_components.size());
                 for (auto h : py_components) {
                     auto& g = h.cast<F2LinearSource&>();
                     comps.push_back(g.clone_source());
                 }
                 return std::make_unique<CombinedF2LinearSource>(
                     std::move(comps), Lmax);
             }),
             py::arg("components"), py::arg("Lmax"),
             "Construct from a list of F2LinearSource objects (Recurrence "
             "PRNGs or DigitalNet sources; each is deep-copied) and a max "
             "output resolution Lmax. The combined output is the XOR of "
             "the components' tempered L-bit outputs.")
        .def(py::init([](py::list py_components,
                          py::list py_tempering,
                          int Lmax) {
                 if (py::len(py_components) != py::len(py_tempering))
                     throw std::invalid_argument(
                         "components and tempering chains must have equal "
                         "length");
                 std::vector<std::unique_ptr<F2LinearSource>> comps;
                 std::vector<TemperingChain> chains;
                 comps.reserve(py_components.size());
                 chains.reserve(py_tempering.size());
                 for (size_t j = 0; j < py_components.size(); ++j) {
                     auto& g = py_components[j].cast<F2LinearSource&>();
                     comps.push_back(g.clone_source());
                     TemperingChain chain;
                     for (auto th : py_tempering[j].cast<py::list>()) {
                         auto& t = th.cast<Transformation&>();
                         chain.add(t.copy());
                     }
                     chains.push_back(std::move(chain));
                 }
                 return std::make_unique<CombinedF2LinearSource>(
                     std::move(comps), std::move(chains), Lmax);
             }),
             py::arg("components"), py::arg("tempering_chains"),
             py::arg("Lmax"),
             "Construct from a list of F2LinearSource objects (Recurrence "
             "PRNGs or DigitalNet sources), a parallel list of per-component "
             "tempering chains (each a list of Transformation objects), and "
             "a max output resolution Lmax.")
        .def("k", &CombinedF2LinearSource::k)
        .def("L", &CombinedF2LinearSource::L)
        .def("J", &CombinedF2LinearSource::J)
        .def("name", &CombinedF2LinearSource::name)
        .def("display_str", &CombinedF2LinearSource::display_str)
        .def("init", &CombinedF2LinearSource::init)
        .def("next", &CombinedF2LinearSource::next)
        .def("get_output", &CombinedF2LinearSource::get_output)
        .def("copy",
             [](const CombinedF2LinearSource& s) { return s.copy(); })
        .def("default_test_method",
             &CombinedF2LinearSource::default_test_method,
             py::arg("test_type"),
             py::call_guard<py::gil_scoped_release>())
        .def("components", &CombinedF2LinearSource::components,
             py::return_value_policy::reference_internal)
        .def("prefix_k", &CombinedF2LinearSource::prefix_k);

    // ── F2LinearSourcePool + ComboEnumerator (Phase 5.2 rename) ────────────
    //
    // Iterator over the cartesian product of J source pools, with
    // identity-uniqueness and shared-pool C(n,k) selection. Replaces
    // (in C++) the iteration engine of
    // regpoly.core.{component,combination}.

    py::class_<F2LinearSourcePool, std::shared_ptr<F2LinearSourcePool>>(
        m, "F2LinearSourcePool")
        .def(py::init<>())
        .def("nb_sources", &F2LinearSourcePool::nb_sources)
        .def("nb_trans",   &F2LinearSourcePool::nb_trans)
        .def("current_index", &F2LinearSourcePool::current_index)
        .def("set_current_index", &F2LinearSourcePool::set_current_index)
        .def("add_source",
             [](F2LinearSourcePool& p, const F2LinearSource& s) { p.add_source(s); },
             py::arg("src"),
             "Append a deep copy of `src` to this pool. `src` may be "
             "any F2LinearSource subclass (Recurrence PRNG or DigitalNet).")
        .def("add_trans",
             [](F2LinearSourcePool& p, const Transformation& t) { p.add_trans(t); },
             py::arg("trans"),
             "Append a deep copy of `trans` to the tempering chain.")
        .def("copy_pool_from", &F2LinearSourcePool::copy_pool_from,
             py::arg("other"),
             "Reference `other`'s pool by shared_ptr (enables C(n,k) "
             "selection across slots).")
        .def("active_source",
             [](F2LinearSourcePool& p) -> F2LinearSource* { return &p.active_source(); },
             py::return_value_policy::reference_internal)
        .def("source_at",
             [](F2LinearSourcePool& p, int i) -> F2LinearSource* { return &p.source_at(i); },
             py::arg("i"),
             py::return_value_policy::reference_internal)
        .def("trans_at",
             [](F2LinearSourcePool& p, int i) -> Transformation* {
                 return &p.trans_at(i);
             },
             py::arg("i"),
             py::return_value_policy::reference_internal)
        .def("display", &F2LinearSourcePool::display);

    py::class_<ComboEnumerator, std::shared_ptr<ComboEnumerator>>(m, "ComboEnumerator")
        .def(py::init<int, int>(), py::arg("J"), py::arg("Lmax"))
        .def_property_readonly("J",    &ComboEnumerator::J)
        .def_property_readonly("Lmax", &ComboEnumerator::Lmax)
        .def_property_readonly("k_g",  &ComboEnumerator::k_g)
        .def_property_readonly("L",    &ComboEnumerator::L)
        .def("pool",
             [](ComboEnumerator& c, int j) -> F2LinearSourcePool* {
                 return &c.pool(j);
             },
             py::arg("j"),
             py::return_value_policy::reference_internal)
        .def("__getitem__",
             [](ComboEnumerator& c, int j) -> F2LinearSource* { return &c.at(j); },
             py::arg("j"),
             py::return_value_policy::reference_internal)
        .def("reset",     &ComboEnumerator::reset)
        .def("next",      &ComboEnumerator::next)
        .def("exhausted", &ComboEnumerator::exhausted);

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

    // Phase 2 Step 6b: chains live on each component intrinsically; the
    // legacy `trans` parameter is no longer needed. The `trans_py` kwarg
    // is still accepted (and ignored) for backwards compatibility, since
    // chains should already be installed on the passed-in generators via
    // `set_tempering()`.
    m.def("prepare_mat",
          [](const py::list& gens_py,
             const std::vector<int>& gen_k,
             const py::list& /*trans_py — ignored, chains on sources*/,
             int kg, int indice_max, int L) {
        std::vector<const F2LinearSource*> gens;
        for (auto item : gens_py)
            gens.push_back(item.cast<const F2LinearSource*>());
        return GaussMatrix::prepare(gens, gen_k, kg, indice_max, L);
    }, py::arg("gens"), py::arg("gen_k"), py::arg("trans"),
       py::arg("kg"), py::arg("indice_max"), py::arg("L"));

    // Polymorphic prepare_mat: takes a single F2LinearSource& and unpacks
    // its sources() internally. Chains are read intrinsically by the
    // kernel via `tempered_output()`. The Python `AbstractTest` base uses
    // this to feed the matricial kernels without going through the
    // list-shape factory.
    m.def("prepare_mat_from_gen",
          [](const ITestable& gen, int indice_max) {
        auto srcs = gen.sources();
        std::vector<int> gen_k;
        gen_k.reserve(srcs.size());
        for (auto* c : srcs) gen_k.push_back(c->k());
        return GaussMatrix::prepare(srcs, gen_k, gen.k(), indice_max, gen.L());
    }, py::arg("gen"), py::arg("indice_max"));

    // ── test_me_lat (dual lattice method) ──────────────────────────────

    // Lattice-family equidistribution kernels: every binding takes a
    // CombinedF2LinearSource. Wrap a bare Generator (J=1) via
    // `regpoly.make_combined(gen, Lmax=...)` before invoking.

    m.def("test_me_lat",
          [](const CombinedF2LinearSource& cs,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_lat(cs, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("cs"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("test_me_harase",
          [](const CombinedF2LinearSource& cs,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_harase(cs, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("cs"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("test_me_notprimitive",
          [](const CombinedF2LinearSource& cs,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_notprimitive(cs, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("cs"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("test_me_notprimitive_simd",
          [](const CombinedF2LinearSource& cs,
             int kg, int L, int maxL,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto result = test_me_notprimitive_simd(cs, kg, L, maxL, delta, mse);
        py::dict d;
        d["ecart"] = result.ecart;
        d["se"] = result.se;
        return d;
    }, py::arg("cs"), py::arg("kg"), py::arg("L"), py::arg("maxL"),
       py::arg("delta"), py::arg("mse"));

    m.def("compute_kv",
          [](const CombinedF2LinearSource& cs, int kg, int v) -> int {
        return compute_kv(cs, kg, v);
    }, py::arg("cs"), py::arg("kg"), py::arg("v"));

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

    m.def("test_me_matricial",
          [](const ITestable& gen, int kg, int L, int Lmax,
             const std::vector<int>& delta, int mse) -> py::dict {
        auto r = test_me_matricial(gen, kg, L, Lmax, delta, mse);
        py::dict d;
        d["ecart"] = r.ecart;
        d["se"] = r.se;
        d["verified"] = r.verified;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("Lmax"),
       py::arg("delta"), py::arg("mse"));

    // ── t-value profile (Phase 4) ──────────────────────────────────────
    //
    // Schmid-style primal enumeration over compositions; dual method
    // is a registered name that raises (placeholder for the
    // Niederreiter–Pirsic kernel that lands in a follow-up).
    m.def("run_tvalue_profile",
          [](const ITestable& gen, int kg, int m, int s_max,
             const std::vector<int>& delta, int max_t_sum,
             const std::string& method) -> py::dict {
        regpoly::core::TValueResult r;
        if (method == "schmid") {
            r = regpoly::core::run_tvalue_profile_schmid(
                gen, kg, m, s_max, delta, max_t_sum);
        } else if (method == "niederreiter_pirsic") {
            r = regpoly::core::run_tvalue_profile_dual(
                gen, kg, m, s_max, delta, max_t_sum);
        } else {
            throw std::invalid_argument(
                "run_tvalue_profile: unknown method '" + method
                + "'. Expected 'schmid' or 'niederreiter_pirsic'.");
        }
        py::dict d;
        d["tvals"] = r.tvals;
        d["se"] = r.se;
        d["verified"] = r.verified;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("m"), py::arg("s_max"),
       py::arg("delta"), py::arg("max_t_sum"), py::arg("method") = "schmid");

    // ── DigitalNet family (Phase 1) + concrete nets (Phases 2-3) ───────
    //
    // Expose the abstract DigitalNet so Python can dynamic_cast / hold
    // refs polymorphically, plus the two concrete first-pass nets.
    py::class_<regpoly::core::DigitalNet, regpoly::core::F2LinearSource>(m, "DigitalNet")
        .def("s_max", &regpoly::core::DigitalNet::s_max)
        .def("m", &regpoly::core::DigitalNet::m)
        .def("current_j", &regpoly::core::DigitalNet::current_j);

    py::class_<regpoly::core::SobolNet, regpoly::core::DigitalNet>(m, "SobolNet")
        .def(py::init<int, int>(),
             py::arg("m"), py::arg("s_max"),
             "Construct a Sobol net using the embedded Joe-Kuo 2003 "
             "trusted table (covers j=2..8 in v1).")
        .def_static("embedded_table_size",
                    &regpoly::core::SobolNet::embedded_table_size);

    py::class_<regpoly::core::NiederreiterF2Net,
               regpoly::core::DigitalNet>(m, "NiederreiterF2Net")
        .def(py::init<int, int>(),
             py::arg("m"), py::arg("s_max"),
             "Construct a Niederreiter F_2 net with irreducibles "
             "generated on the fly.")
        .def("irreducible",
             &regpoly::core::NiederreiterF2Net::irreducible)
        .def("irreducible_degree",
             &regpoly::core::NiederreiterF2Net::irreducible_degree);

    m.def("run_collision_free",
          [](const ITestable& gen, int kg, int L, int L_for_phi4) -> py::dict {
        auto r = run_collision_free(gen, kg, L, L_for_phi4);
        py::dict d;
        d["ecart_cf"] = r.ecart_cf;
        d["secf"] = r.secf;
        d["verified"] = r.verified;
        return d;
    }, py::arg("gen"), py::arg("kg"), py::arg("L"), py::arg("L_for_phi4"));

    m.def("run_tuplets",
          [](const ITestable& gen, int kg, int L, int tupd,
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

    // ── TemperingOptimizerCache (dual lattice StackBase for optimizer) ─────────

    py::class_<TemperingOptimizerCache>(m, "TemperingOptimizerCache")
        .def(py::init([](const CombinedF2LinearSource& cs, int kg, int L) {
            return TemperingOptimizerCache(cs, kg, L);
        }), py::arg("cs"), py::arg("kg"), py::arg("L"),
        py::keep_alive<1, 2>())  // cs must outlive the cache
        .def("compute_all", &TemperingOptimizerCache::compute_all)
        .def("compute_gap", &TemperingOptimizerCache::compute_gap, py::arg("v"))
        .def("refresh_inv_g0", &TemperingOptimizerCache::refresh_inv_g0)
        .def("rebuild", &TemperingOptimizerCache::rebuild)
        .def("reset_step", &TemperingOptimizerCache::reset_step)
        .def("step", &TemperingOptimizerCache::step, py::arg("v"));

    // ── HaraseRankCache (StackBase strategy for tempering optimizer) ──────────

    py::class_<HaraseRankCache>(m, "HaraseRankCache")
        .def(py::init([](const CombinedF2LinearSource& cs, int kg, int L) {
            return HaraseRankCache(cs, kg, L);
        }), py::arg("cs"), py::arg("kg"), py::arg("L"))
        .def("compute_all", &HaraseRankCache::compute_all)
        .def("restore_and_reduce", &HaraseRankCache::restore_and_reduce,
             py::arg("v"))
        .def("kg", &HaraseRankCache::kg)
        .def("L", &HaraseRankCache::L);

    // ── tausworthe_random_poly: sample an admissible polynomial ──────

    m.def("tausworthe_random_poly",
          [](int k, int nb_terms, bool quicktaus, int L, int s) {
        return TauswortheGen::random_poly(k, nb_terms, quicktaus, L, s);
    }, py::arg("k"), py::arg("nb_terms"), py::arg("quicktaus"),
       py::arg("L"), py::arg("s") = 0,
       "Sample a random admissible TauswortheGen polynomial.  Returns the "
       "sorted exponent list [0, q_1, ..., q_{t-2}, k].  Throws if the "
       "(k, nb_terms, quicktaus, L, s) combination is inadmissible.");

    // ── well_random_matrices / well_total_cost ────────────────────────
    //
    // Free-function bindings around the WELL cost-bounded sampler. Used
    // by the Python-side primitive search worker (web app) and by the
    // `regpoly.well` Python module. The C++ search driver hits the same
    // `WELLGen::random_matrices` static directly without going through
    // these helpers.

    m.def("well_random_matrices",
          [](int w, int max_cost, uint64_t seed) -> py::dict {
        std::mt19937_64 rng{seed};
        StructMap sm = WELLGen::random_matrices(w, max_cost, rng);
        py::dict outer;
        for (const auto& slot : sm) {
            py::dict inner;
            for (const auto& arg : slot.second) {
                inner[py::str(arg.first)] = scalar_to_py(arg.second);
            }
            outer[py::str(slot.first)] = inner;
        }
        return outer;
    }, py::arg("w"), py::arg("max_cost"), py::arg("seed") = 0,
       "Sample a WELL `matrices` map (slots T0..T7) whose total per-Mi "
       "cost is <= max_cost. Uses rejection sampling with a "
       "greedy-budgeted fallback. `seed` is consumed by a per-call "
       "std::mt19937_64 for reproducibility. Throws "
       "std::invalid_argument if max_cost <= 0 or w != 32.");

    m.def("well_total_cost",
          [](const py::dict& matrices) -> int {
        int sum = 0;
        for (auto kv : matrices) {
            std::string slot = kv.first.cast<std::string>();
            if (!py::isinstance<py::dict>(kv.second))
                throw std::runtime_error(
                    "well_total_cost: matrices['" + slot
                    + "'] must be a dict");
            auto inner = kv.second.cast<py::dict>();
            if (!inner.contains("M"))
                throw std::runtime_error(
                    "well_total_cost: matrices['" + slot
                    + "'] is missing required key 'M'");
            int Mi = py::cast<int>(inner["M"]);
            sum += WELLGen::static_cost_for_Mi(Mi);
        }
        return sum;
    }, py::arg("matrices"),
       "Sum of per-Mi costs across the slots in a `matrices` dict. "
       "Each slot value must be a dict carrying an integer 'M' key "
       "selecting the M-class (paper Table I).");

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
        ParamBag p = dict_to_params(params_dict);
        RandomParamResult r;
        if (rand_type == "tausworthe_s"
            || rand_type == "tausworthe_poly") {
            r = TauswortheGen::generate_random(
                rand_type, rand_args, p, L);
        } else {
            // Fall back to the generic sampler (handles bitmask, range,
            // poly_exponents, bitmask_vec, irreducible_gf2 — anything
            // listed in random_samplers.cpp's `sample_generic_into`).
            // The Python `parametric.generate_random` shadow-implements
            // the older subset for speed and reaches this binding for
            // anything it doesn't recognise locally.
            ParamSpec spec;
            spec.name = "__random_param_tmp";
            spec.rand_type = rand_type;
            spec.rand_args = rand_args;
            if (!regpoly::random::sample_generic_into(spec, p))
                throw std::invalid_argument(
                    "random_param: unsupported rand_type '"
                    + rand_type + "'");
            // Pull the sampled value back out of the temporary slot.
            auto iv_it = p.int_vecs().find(spec.name);
            auto uv_it = p.uint_vecs().find(spec.name);
            auto i_it = p.ints().find(spec.name);
            if (iv_it != p.int_vecs().end()) {
                r.is_vec = true;
                r.vec_val.assign(iv_it->second.begin(),
                                 iv_it->second.end());
            } else if (uv_it != p.uint_vecs().end()) {
                r.is_vec = true;
                r.vec_val.reserve(uv_it->second.size());
                for (uint64_t v : uv_it->second)
                    r.vec_val.push_back(static_cast<int64_t>(v));
            } else if (i_it != p.ints().end()) {
                r.is_vec = false;
                r.int_val = i_it->second;
            } else {
                throw std::runtime_error(
                    "random_param: sampler returned no value for '"
                    + rand_type + "'");
            }
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

    py::class_<ParameterSpaceEnumerator, std::shared_ptr<ParameterSpaceEnumerator>>(m, "ParameterSpaceEnumerator")
        .def("size", [](const ParameterSpaceEnumerator& e) {
            // Lift decimal-string count to a Python int of arbitrary precision.
            return py::module_::import("builtins").attr("int")(e.size_dec());
        })
        .def("at", [](const ParameterSpaceEnumerator& e, const py::object& idx) {
            return params_to_dict(e.at(py::str(idx).cast<std::string>()));
        }, py::arg("idx"))
        .def("axes", [](const ParameterSpaceEnumerator& e) {
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

    m.def("make_parameter_space_enumerator",
          [](const std::string& family, const py::dict& d, int L) -> py::object {
        auto e = make_parameter_space_enumerator(family, dict_to_params(d), L);
        if (!e) return py::none();
        std::shared_ptr<ParameterSpaceEnumerator> shared(std::move(e));
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

    // ── Recurrence subclass registrations ────────────────────────────────
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
          [](const Recurrence& gen) -> bool { return is_full_period(gen); },
          py::arg("gen"),
          "True iff the generator's characteristic polynomial is "
          "primitive (period 2^k - 1).");

    m.def("is_full_period",
          [](const BitVect& char_poly, int k) -> bool {
              return is_full_period(char_poly, k);
          },
          py::arg("char_poly"), py::arg("k"),
          "Variant: test primitivity of a candidate char-poly directly.  "
          "Bit j of char_poly is the coefficient of x^j (leading x^k is "
          "implicit).");

    m.def("packed_bm",
          [](const ITestable& gen, const BitVect& init_state, int K,
             int bit_idx) -> std::pair<int, BitVect> {
              BitVect mp;
              int L = packed_bm(gen, init_state, K, &mp, bit_idx);
              return {L, mp};
          },
          py::arg("gen"), py::arg("init_state"), py::arg("K"),
          py::arg("bit_idx") = 0,
          "Run packed Berlekamp-Massey over F_2 on output bit `bit_idx` of "
          "`gen` after init(init_state).  Returns (linear_complexity, "
          "min_poly_BitVect).  min_poly is K bits, MSB-first (same "
          "convention as Recurrence::char_poly()).");

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
        .def_readwrite("max_cost", &PrimitiveSearchConfig::max_cost)
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

    py::class_<TemperingOptimizerResult>(m, "TemperingOptimizerResult")
        .def_readonly("se",               &TemperingOptimizerResult::se)
        .def_readonly("gaps",             &TemperingOptimizerResult::gaps)
        .def_readonly("elapsed_seconds",  &TemperingOptimizerResult::elapsed_seconds)
        .def_readonly("essais",           &TemperingOptimizerResult::essais);

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
              TemperingOptimizerCache& cache,
              const py::list& param_tuples,
              const std::vector<std::vector<uint64_t>>& safe_masks) -> py::tuple {
        auto params = build_locators(param_tuples);
        auto result = run_tempering_optimizer_once(
            cfg, cache, params, safe_masks);
        return py::make_tuple(result, export_locators(params));
    }, py::arg("config"), py::arg("cache"), py::arg("params"),
       py::arg("safe_masks"),
       "Single recursive optimization pass. Returns "
       "(TemperingOptimizerResult, [final_value_per_locator]).");

    // ── SeekDriver (Phase 2.4b) ────────────────────────────────────────
    //
    // Drives the equidistribution / collision-free / tuplets search
    // loop in C++. The Python Seek.run() shrinks to: build a C++
    // ComboEnumerator from the Python ComboEnumerator, build a list of
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
          []() { return EquidistributionMethodRegistry::names(); },
          "List the equidistribution method names known to "
          "EquidistributionMethodRegistry. The Python YAML parser consults this list "
          "rather than maintaining its own parallel enum.");
    m.def("has_equidistribution_method",
          [](const std::string& name) { return EquidistributionMethodRegistry::has(name); },
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
          [](ComboEnumerator& comb,
             const std::vector<SeekTestSpec>& tests,
             int nbtries,
             int progress_interval,
             const py::object& on_prep,
             const py::object& on_iter,
             const py::object& on_progress) -> SeekResult {
        SeekOnPrepFn prep_fn = nullptr;
        if (!on_prep.is_none()) {
            prep_fn = [&on_prep](ComboEnumerator& c, bool is_retry) {
                on_prep(py::cast(&c, py::return_value_policy::reference),
                        is_retry);
            };
        }
        SeekOnIterFn iter_fn = nullptr;
        if (!on_iter.is_none()) {
            iter_fn = [&on_iter](ComboEnumerator& c, const SeekIterResult& r) {
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

    py::class_<TemperingTryResult>(m, "TemperingTryResult")
        .def(py::init<>())
        .def_readwrite("got_result", &TemperingTryResult::got_result)
        .def_readwrite("score",      &TemperingTryResult::score);

    py::class_<TemperingSearchResult>(m, "TemperingSearchResult")
        .def_readonly("nbgen",           &TemperingSearchResult::nbgen)
        .def_readonly("nb_with_result",  &TemperingSearchResult::nb_with_result)
        .def_readonly("elapsed_seconds", &TemperingSearchResult::elapsed_seconds);

    m.def("run_tempering_search",
          [](ComboEnumerator& comb,
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
            start_fn = [&on_combo_start](ComboEnumerator& c, int idx) {
                on_combo_start(
                    py::cast(&c, py::return_value_policy::reference), idx);
            };
        }
        TempSearchOnTryFn try_fn =
            [&on_try](ComboEnumerator& c, int combo_idx, int try_idx,
                      bool is_first) -> TemperingTryResult {
                py::object r = on_try(
                    py::cast(&c, py::return_value_policy::reference),
                    combo_idx, try_idx, is_first);
                return r.cast<TemperingTryResult>();
            };
        TempSearchOnComboDoneFn done_fn = nullptr;
        if (!on_combo_done.is_none()) {
            done_fn = [&on_combo_done](ComboEnumerator& c, int idx,
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
       "returns a TemperingTryResult). Tracks the best score per combo, "
       "fires on_combo_done with the best result, advances the "
       "ComboEnumerator, and emits on_progress every progress_interval tries. "
       "Returns a TemperingSearchResult.");

    m.def("run_tempering_optimizer_minimize",
          [build_locators, export_locators](
              const TemperingOptimizerConfig& cfg,
              TemperingOptimizerCache& cache,
              const py::list& param_tuples,
              const std::vector<std::vector<uint64_t>>& safe_masks) -> py::tuple {
        auto params = build_locators(param_tuples);
        auto result = run_tempering_optimizer_minimize(
            cfg, cache, params, safe_masks);
        return py::make_tuple(result, export_locators(params));
    }, py::arg("config"), py::arg("cache"), py::arg("params"),
       py::arg("safe_masks"),
       "Iterative delta-tightening loop (n_restarts > 1). Returns "
       "(TemperingOptimizerResult, [final_value_per_locator]).");

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

    using namespace regpoly::library;

    // ParamValue → native Python value (int / str / bool / list[int]
    // / nested dict-of-dicts for StructMap-shaped fields like WELL
    // `matrices`). Without the StructMap case, Python sees `matrices:
    // None` and any downstream Recurrence.create call rejects the
    // entry — see docs/generators/WELLGen.md.
    auto pv_to_py = [](const ParamValue& v) -> py::object {
        switch (v.kind) {
            case ParamKind::Int:    return py::int_(v.int_val);
            case ParamKind::String: return py::str(v.string_val);
            case ParamKind::Bool:   return py::bool_(v.bool_val);
            case ParamKind::IntList: return py::cast(v.int_list_val);
            case ParamKind::StructMap: {
                py::dict outer;
                for (const auto& slot : v.struct_map_val) {
                    py::dict inner;
                    for (const auto& arg : slot.second) {
                        inner[py::str(arg.first)] = scalar_to_py(arg.second);
                    }
                    outer[py::str(slot.first)] = inner;
                }
                return std::move(outer);
            }
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
        const regpoly::library::CatalogComponent& c) -> py::dict {
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
              return regpoly::library::config_hash(family, pm, steps);
          },
          py::arg("family"), py::arg("params"), py::arg("tempering"),
          "Stable short hash of one component config.");
}
