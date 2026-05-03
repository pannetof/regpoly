// Phase 2.4: pybind11 registrations for the generator family classes.
// Split out of factory.cpp so the factory entry points
// (create_generator, get_gen_param_specs, …) can live in regpoly_core
// and be linked from the SearchDriver family without dragging
// pybind11 in.

#include "factory.h"

#include "ac1d.h"
#include "dsfmt.h"
#include "f2w_base.h"
#include "f2w_lfsr.h"
#include "f2w_polylcg.h"
#include "marsaxorshift.h"
#include "matsumoto.h"
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
#include "xorshift128.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;

// Each family is registered under its new canonical name; old names
// remain accessible via aliases set up below so existing YAML and
// SQLite rows continue to dispatch.
void register_generator_types(py::module_& m) {
    py::class_<PolyLCGGen, Generator, std::unique_ptr<PolyLCGGen>>(m, "PolyLCGGen");
    py::class_<TauswortheGen, Generator, std::unique_ptr<TauswortheGen>>(m, "TauswortheGen");
    py::class_<TGFSRGen, Generator, std::unique_ptr<TGFSRGen>>(m, "TGFSRGen");
    py::class_<MTGen, Generator, std::unique_ptr<MTGen>>(m, "MTGen");
    py::class_<F2wBaseGen, Generator, std::unique_ptr<F2wBaseGen>>(m, "F2wBaseGen");
    py::class_<F2wPolyLCGGen, F2wBaseGen, std::unique_ptr<F2wPolyLCGGen>>(m, "F2wPolyLCGGen");
    py::class_<F2wLFSRGen, F2wBaseGen, std::unique_ptr<F2wLFSRGen>>(m, "F2wLFSRGen");
    py::class_<MatsumotoGen, Generator, std::unique_ptr<MatsumotoGen>>(m, "MatsumotoGen");
    py::class_<MarsaXorshiftGen, Generator, std::unique_ptr<MarsaXorshiftGen>>(m, "MarsaXorshiftGen");
    py::class_<AC1DGen, Generator, std::unique_ptr<AC1DGen>>(m, "AC1DGen");
    py::class_<WELLGen, Generator, std::unique_ptr<WELLGen>>(m, "WELLGen");
    py::class_<MELGGen, Generator, std::unique_ptr<MELGGen>>(m, "MELGGen");
    py::class_<SFMTGen, Generator, std::unique_ptr<SFMTGen>>(m, "SFMTGen");
    py::class_<DSFMTGen, Generator, std::unique_ptr<DSFMTGen>>(m, "DSFMTGen");
    py::class_<MTGPGen, Generator, std::unique_ptr<MTGPGen>>(m, "MTGPGen");
    py::class_<XorShift128Gen, Generator, std::unique_ptr<XorShift128Gen>>(m, "XorShift128Gen");
    py::class_<TinyMT32Gen, Generator, std::unique_ptr<TinyMT32Gen>>(m, "TinyMT32Gen");
    py::class_<RMT64Gen, Generator, std::unique_ptr<RMT64Gen>>(m, "RMT64Gen");

    // Backwards-compat: expose every renamed family under its old Python-visible
    // name as well, pointing at the same class object.
    m.attr("PolyLCG")         = m.attr("PolyLCGGen");
    m.attr("Tausworthe")      = m.attr("TauswortheGen");
    m.attr("TGFSR")           = m.attr("TGFSRGen");
    m.attr("MersenneTwister") = m.attr("MTGen");
    m.attr("GenF2wBase")      = m.attr("F2wBaseGen");
    m.attr("GenF2wLFSR")      = m.attr("F2wLFSRGen");
    m.attr("GenF2wPolyLCG")   = m.attr("F2wPolyLCGGen");
    m.attr("WELLRNG")         = m.attr("WELLGen");
    m.attr("MELG")            = m.attr("MELGGen");
    m.attr("SFMT")            = m.attr("SFMTGen");
    m.attr("dSFMTGen")        = m.attr("DSFMTGen");
    m.attr("MTGP")            = m.attr("MTGPGen");
    m.attr("XorShift128")     = m.attr("XorShift128Gen");
    m.attr("TinyMT32")        = m.attr("TinyMT32Gen");
    m.attr("RMT64")           = m.attr("RMT64Gen");
}
