#include "factory.h"
#include "polylcg.h"
#include "tausworthe.h"
#include "tgfsr.h"
#include "mt.h"
#include "f2w_base.h"
#include "f2w_polylcg.h"
#include "f2w_lfsr.h"
#include "matsumoto.h"
#include "marsaxorshift.h"
#include "ac1d.h"
#include "well.h"
#include "melg.h"
#include "sfmt.h"
#include "dsfmt.h"
#include "mtgp.h"
#include "xorshift128.h"
#include "tinymt32.h"
#include "rmt64.h"
#include "permutation.h"
#include "temper_mk.h"
#include "lag_mask.h"
#include <stdexcept>

// Phase 2.4: factory.cpp is now pure C++ — pybind11 wiring lives in
// src/bindings/factory_bindings.cpp so factory functions can be linked
// into regpoly_core (and used by the SearchDriver family).

std::unique_ptr<Generator> create_generator(
    const std::string& family, const Params& params, int L)
{
    if (family == "PolyLCG"         || family == "PolyLCGGen")        return PolyLCGGen::from_params(params, L);
    if (family == "Tausworthe"      || family == "TauswortheGen")     return TauswortheGen::from_params(params, L);
    if (family == "TGFSR"           || family == "TGFSRGen")          return TGFSRGen::from_params(params, L);
    if (family == "MersenneTwister" || family == "MTGen")             return MTGen::from_params(params, L);
    if (family == "GenF2wPolyLCG"   || family == "F2wPolyLCGGen")     return F2wPolyLCGGen::from_params(params, L);
    if (family == "GenF2wLFSR"      || family == "F2wLFSRGen")        return F2wLFSRGen::from_params(params, L);
    if (family == "MatsumotoGen")                                     return MatsumotoGen::from_params(params, L);
    if (family == "MarsaXorshiftGen")                                 return MarsaXorshiftGen::from_params(params, L);
    if (family == "AC1DGen")                                          return AC1DGen::from_params(params, L);
    if (family == "WELLRNG"         || family == "WELLGen")           return WELLGen::from_params(params, L);
    if (family == "MELG"            || family == "MELGGen")           return MELGGen::from_params(params, L);
    if (family == "SFMT"            || family == "SFMTGen")           return SFMTGen::from_params(params, L);
    if (family == "dSFMTGen"        || family == "DSFMTGen")          return DSFMTGen::from_params(params, L);
    if (family == "MTGP"            || family == "MTGPGen")           return MTGPGen::from_params(params, L);
    if (family == "XorShift128"     || family == "XorShift128Gen")    return XorShift128Gen::from_params(params, L);
    if (family == "TinyMT32"        || family == "TinyMT32Gen")       return TinyMT32Gen::from_params(params, L);
    if (family == "RMT64"           || family == "RMT64Gen")          return RMT64Gen::from_params(params, L);
    throw std::invalid_argument("Unknown generator family: " + family);
}

std::unique_ptr<Transformation> create_transformation(
    const std::string& type, const Params& params)
{
    if (type == "permut")                        return PermutationTrans::from_params(params);
    if (type == "tempMK" || type == "tempMK2")   return TemperMKTrans::from_params(type, params);
    if (type == "laggedTempering")               return LaggedTempering::from_params(params);
    throw std::invalid_argument("Unknown transformation type: " + type);
}

std::vector<ParamSpec> get_gen_param_specs(const std::string& family)
{
    if (family == "PolyLCG"         || family == "PolyLCGGen")        return PolyLCGGen::param_specs();
    if (family == "Tausworthe"      || family == "TauswortheGen")     return TauswortheGen::param_specs();
    if (family == "TGFSR"           || family == "TGFSRGen")          return TGFSRGen::param_specs();
    if (family == "MersenneTwister" || family == "MTGen")             return MTGen::param_specs();
    if (family == "GenF2wPolyLCG"   || family == "F2wPolyLCGGen"
                                    || family == "GenF2wLFSR"
                                    || family == "F2wLFSRGen")        return F2wPolyLCGGen::param_specs();
    if (family == "MatsumotoGen")                                     return MatsumotoGen::param_specs();
    if (family == "MarsaXorshiftGen")                                 return MarsaXorshiftGen::param_specs();
    if (family == "AC1DGen")                                          return AC1DGen::param_specs();
    if (family == "WELLRNG"         || family == "WELLGen")           return WELLGen::param_specs();
    if (family == "MELG"            || family == "MELGGen")           return MELGGen::param_specs();
    if (family == "SFMT"            || family == "SFMTGen")           return SFMTGen::param_specs();
    if (family == "dSFMTGen"        || family == "DSFMTGen")          return DSFMTGen::param_specs();
    if (family == "MTGP"            || family == "MTGPGen")           return MTGPGen::param_specs();
    if (family == "XorShift128"     || family == "XorShift128Gen")    return XorShift128Gen::param_specs();
    if (family == "TinyMT32"        || family == "TinyMT32Gen")       return TinyMT32Gen::param_specs();
    if (family == "RMT64"           || family == "RMT64Gen")          return RMT64Gen::param_specs();
    throw std::invalid_argument("Unknown generator family: " + family);
}

std::vector<ParamSpec> get_trans_param_specs(const std::string& type)
{
    if (type == "permut")                        return PermutationTrans::param_specs();
    if (type == "tempMK" || type == "tempMK2")   return TemperMKTrans::param_specs();
    if (type == "laggedTempering")               return LaggedTempering::param_specs();
    throw std::invalid_argument("Unknown transformation type: " + type);
}

// ── Exhaustive-search enumerator registry ─────────────────────────────────

bool family_is_enumerable(const std::string& family)
{
    if (family == "Tausworthe" || family == "TauswortheGen") return true;
    return false;
}

std::unique_ptr<GenEnumerator> make_gen_enumerator(
    const std::string& family, const Params& resolved, int L)
{
    if (family == "Tausworthe" || family == "TauswortheGen")
        return TauswortheGen::make_enumerator(resolved, L);
    return nullptr;
}
