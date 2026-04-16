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
#include "gen_wellrng.h"
#include "gen_melg.h"
#include "trans_permutation.h"
#include "trans_temper_mk.h"
#include "trans_lag_mask.h"
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
    py::class_<WELLRNG, Generateur, std::unique_ptr<WELLRNG>>(m, "WELLRNG");
    py::class_<MELG, Generateur, std::unique_ptr<MELG>>(m, "MELG");
}

std::unique_ptr<Generateur> create_generator(
    const std::string& family, const Params& params, int L)
{
    if (family == "PolyLCG")          return PolyLCG::from_params(params, L);
    if (family == "Tausworthe")       return Tausworthe::from_params(params, L);
    if (family == "TGFSRGen")         return TGFSRGen::from_params(params, L);
    if (family == "MersenneTwister")  return MersenneTwister::from_params(params, L);
    if (family == "GenF2wPolyLCG")    return GenF2wPolyLCG::from_params(params, L);
    if (family == "GenF2wLFSR")       return GenF2wLFSR::from_params(params, L);
    if (family == "MatsumotoGen")     return MatsumotoGen::from_params(params, L);
    if (family == "MarsaXorshiftGen") return MarsaXorshiftGen::from_params(params, L);
    if (family == "AC1DGen")          return AC1DGen::from_params(params, L);
    if (family == "WELLRNG")          return WELLRNG::from_params(params, L);
    if (family == "MELG")             return MELG::from_params(params, L);
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
    if (family == "PolyLCG")          return PolyLCG::param_specs();
    if (family == "Tausworthe")       return Tausworthe::param_specs();
    if (family == "TGFSRGen")         return TGFSRGen::param_specs();
    if (family == "MersenneTwister")  return MersenneTwister::param_specs();
    if (family == "GenF2wPolyLCG" || family == "GenF2wLFSR")
                                      return GenF2wPolyLCG::param_specs();
    if (family == "MatsumotoGen")     return MatsumotoGen::param_specs();
    if (family == "MarsaXorshiftGen") return MarsaXorshiftGen::param_specs();
    if (family == "AC1DGen")          return AC1DGen::param_specs();
    if (family == "WELLRNG")          return WELLRNG::param_specs();
    if (family == "MELG")             return MELG::param_specs();
    throw std::invalid_argument("Unknown generator family: " + family);
}

std::vector<ParamSpec> get_trans_param_specs(const std::string& type)
{
    if (type == "permut")                        return PermutationTrans::param_specs();
    if (type == "tempMK" || type == "tempMK2")   return TemperMKTrans::param_specs();
    if (type == "laggedTempering")               return LaggedTempering::param_specs();
    throw std::invalid_argument("Unknown transformation type: " + type);
}
