// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Pybind11 registration for Generator subclass types.
//
// The `py::class_<Cls, Generator>(m, name)` instantiations are
// templated on the concrete C++ type, so they must be enumerated in
// source — there is no way to call them generically from a runtime
// registry. Adding a new family means appending ONE line to
// `register_generator_types` below.
//
// Aliases for renamed families are pulled from GeneratorRegistry, so
// they live alongside the canonical registration in factory.cpp and
// do NOT need to be duplicated here.

#include "factory.h"
#include "generator_registry.h"

#include "cellular_automata.h"
#include "dsfmt.h"
#include "f2w_base.h"
#include "f2w_lfsr.h"
#include "f2w_polylcg.h"
#include "marsaxorshift.h"
#include "melg.h"
#include "mt.h"
#include "mtgp.h"
#include "polylcg.h"
#include "rmt64.h"
#include "sfmt.h"
#include "tausworthe.h"
#include "tgfsr.h"
#include "tinymt32.h"
#include "well.h"
#include "xoroshiro.h"
#include "xoshiro.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace {

// Forward declaration: defined alongside the canonical generator
// registrations in factory.cpp; calling it here ensures the registry
// is populated before we walk it for aliases.
void ensure_registered() {
    // Trigger the lazy initialiser inside factory.cpp by asking for any
    // canonical entry. lookup() throws if the registry is empty, but
    // create_generator with any known name has the same effect. The
    // simplest portable hook is to ask through the public factory:
    // any get_gen_param_specs call performs register_all_generators().
    (void)get_gen_param_specs("PolyLCGGen");
}

// One-line helper for the common case (subclass of Generator).
template <typename Cls>
void bind_gen(py::module_& m, const char* name) {
    py::class_<Cls, Generator, std::unique_ptr<Cls>>(m, name);
}

}  // namespace

void register_generator_types(py::module_& m) {
    // ── Concrete classes ───────────────────────────────────────────────
    bind_gen<PolyLCGGen>      (m, "PolyLCGGen");
    bind_gen<TauswortheGen>   (m, "TauswortheGen");
    bind_gen<TGFSRGen>        (m, "TGFSRGen");
    bind_gen<MTGen>           (m, "MTGen");
    py::class_<F2wBaseGen, Generator, std::unique_ptr<F2wBaseGen>>(m, "F2wBaseGen");
    py::class_<F2wPolyLCGGen, F2wBaseGen, std::unique_ptr<F2wPolyLCGGen>>(m, "F2wPolyLCGGen");
    py::class_<F2wLFSRGen,    F2wBaseGen, std::unique_ptr<F2wLFSRGen>>   (m, "F2wLFSRGen");
    bind_gen<MarsaXorshiftGen>(m, "MarsaXorshiftGen");
    bind_gen<CellularAutomataGen>(m, "CellularAutomataGen");
    py::class_<WELLGen, Generator, std::unique_ptr<WELLGen>>(m, "WELLGen")
        .def("total_cost", &WELLGen::total_cost,
             "Sum of per-Mi costs across the 8 algorithm slots T0..T7.");
    bind_gen<MELGGen>         (m, "MELGGen");
    bind_gen<SFMTGen>         (m, "SFMTGen");
    bind_gen<DSFMTGen>        (m, "DSFMTGen");
    bind_gen<MTGPGen>         (m, "MTGPGen");
    bind_gen<TinyMT32Gen>     (m, "TinyMT32Gen");
    bind_gen<RMT64Gen>        (m, "RMT64Gen");
    bind_gen<XoroshiroGen>    (m, "XoroshiroGen");
    bind_gen<XoshiroGen>      (m, "XoshiroGen");

    // F2wBase legacy alias is a class, not a registry alias, so it is
    // declared here directly (the registry handles only concrete
    // creatable generators).
    m.attr("GenF2wBase") = m.attr("F2wBaseGen");

    // ── Aliases — pulled from the central registry ─────────────────────
    // factory.cpp's register_all_generators() declares every family
    // rename via GR::reg_alias. We mirror those into the Python module
    // by walking the registry rather than maintaining a duplicate list.
    ensure_registered();
    for (const auto& [alias, canonical] : GeneratorRegistry::aliases()) {
        if (py::hasattr(m, canonical.c_str())) {
            m.attr(alias.c_str()) = m.attr(canonical.c_str());
        }
    }
}
