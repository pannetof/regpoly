// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "factory.h"
#include "generator_registry.h"
#include "transformation_registry.h"

#include "cellular_automata.h"
#include "dsfmt.h"
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

#include "lag_mask.h"
#include "permutation.h"
#include "temper_mk.h"

#include <stdexcept>

using namespace regpoly::core;


// ── Central registration of every concrete Generator and Transformation ─
//
// Adding a new generator family means:
//   1. Write foo.h / foo.cpp with the Generator subclass.
//   2. Add ONE line in this file: `GR::reg("FooGen", ...)`.
//   3. Add ONE line in factory_bindings.cpp: a `py::class_<FooGen,
//      Generator>(m, "FooGen")` registration.
// (Plus aliases here if the family was renamed.)
//
// The factory entry points below (create_generator, get_gen_param_specs,
// family_is_enumerable, make_gen_enumerator, create_transformation,
// get_trans_param_specs) are pure registry lookups. There is no
// per-family if-chain to maintain.

namespace regpoly::core {

namespace {

using GR = GeneratorRegistry;

void register_all_generators() {
    static const int once = []{
        GR::reg("PolyLCGGen",
                &PolyLCGGen::from_params, &PolyLCGGen::param_specs);
        GR::reg_alias("PolyLCG", "PolyLCGGen");

        GR::reg("TauswortheGen",
                &TauswortheGen::from_params, &TauswortheGen::param_specs,
                /*bind=*/{},
                &TauswortheGen::make_enumerator);
        GR::reg_alias("Tausworthe", "TauswortheGen");

        GR::reg("TGFSRGen",
                &TGFSRGen::from_params, &TGFSRGen::param_specs);
        GR::reg_alias("TGFSR", "TGFSRGen");

        GR::reg("MTGen",
                &MTGen::from_params, &MTGen::param_specs);
        GR::reg_alias("MersenneTwister", "MTGen");

        GR::reg("F2wPolyLCGGen",
                &F2wPolyLCGGen::from_params, &F2wPolyLCGGen::param_specs);
        GR::reg_alias("GenF2wPolyLCG", "F2wPolyLCGGen");

        GR::reg("F2wLFSRGen",
                &F2wLFSRGen::from_params, &F2wLFSRGen::param_specs);
        GR::reg_alias("GenF2wLFSR", "F2wLFSRGen");

        GR::reg("MarsaXorshiftGen",
                &MarsaXorshiftGen::from_params, &MarsaXorshiftGen::param_specs,
                /*bind=*/{},
                &MarsaXorshiftGen::make_enumerator);

        GR::reg("CellularAutomataGen",
                &CellularAutomataGen::from_params,
                &CellularAutomataGen::param_specs);
        GR::reg_alias("CA", "CellularAutomataGen");

        GR::reg("WELLGen",
                &WELLGen::from_params, &WELLGen::param_specs);
        GR::reg_alias("WELLRNG", "WELLGen");

        GR::reg("MELGGen",
                &MELGGen::from_params, &MELGGen::param_specs);
        GR::reg_alias("MELG", "MELGGen");

        GR::reg("SFMTGen",
                &SFMTGen::from_params, &SFMTGen::param_specs);
        GR::reg_alias("SFMT", "SFMTGen");

        GR::reg("DSFMTGen",
                &DSFMTGen::from_params, &DSFMTGen::param_specs);
        GR::reg_alias("dSFMTGen", "DSFMTGen");

        GR::reg("MTGPGen",
                &MTGPGen::from_params, &MTGPGen::param_specs);
        GR::reg_alias("MTGP", "MTGPGen");

        GR::reg("TinyMT32Gen",
                &TinyMT32Gen::from_params, &TinyMT32Gen::param_specs);
        GR::reg_alias("TinyMT32", "TinyMT32Gen");

        GR::reg("RMT64Gen",
                &RMT64Gen::from_params, &RMT64Gen::param_specs);
        GR::reg_alias("RMT64", "RMT64Gen");

        GR::reg("XoroshiroGen",
                &XoroshiroGen::from_params, &XoroshiroGen::param_specs,
                /*bind=*/{},
                &XoroshiroGen::make_enumerator);
        GR::reg_alias("Xoroshiro", "XoroshiroGen");

        GR::reg("XoshiroGen",
                &XoshiroGen::from_params, &XoshiroGen::param_specs,
                /*bind=*/{},
                &XoshiroGen::make_enumerator);
        GR::reg_alias("Xoshiro", "XoshiroGen");
        return 0;
    }();
    (void)once;
}

using TR = TransformationRegistry;

void register_all_transformations() {
    static const int once = []{
        TR::reg("permut",
                [](const Params& p){ return PermutationTrans::from_params(p); },
                &PermutationTrans::param_specs);
        TR::reg("tempMK",
                [](const Params& p){ return TemperMKTrans::from_params("tempMK", p); },
                &TemperMKTrans::param_specs);
        TR::reg("tempMK2",
                [](const Params& p){ return TemperMKTrans::from_params("tempMK2", p); },
                &TemperMKTrans::param_specs);
        TR::reg("laggedTempering",
                [](const Params& p){ return LaggedTempering::from_params(p); },
                &LaggedTempering::param_specs);
        return 0;
    }();
    (void)once;
}

}  // namespace

// ── Public factory entry points (pure registry lookups) ─────────────────

std::unique_ptr<Generator> create_generator(
    const std::string& family, const Params& params, int L)
{
    register_all_generators();
    return GR::lookup(family).from_params(params, L);
}

std::vector<ParamSpec> get_gen_param_specs(const std::string& family)
{
    register_all_generators();
    return GR::lookup(family).param_specs();
}

bool family_is_enumerable(const std::string& family)
{
    register_all_generators();
    auto* info = GR::find(family);
    return info && static_cast<bool>(info->make_enumerator);
}

std::unique_ptr<GenEnumerator> make_gen_enumerator(
    const std::string& family, const Params& resolved, int L)
{
    register_all_generators();
    auto* info = GR::find(family);
    if (!info || !info->make_enumerator) return nullptr;
    return info->make_enumerator(resolved, L);
}

std::unique_ptr<Transformation> create_transformation(
    const std::string& type, const Params& params)
{
    register_all_transformations();
    return TR::lookup(type).from_params(params);
}

std::vector<ParamSpec> get_trans_param_specs(const std::string& type)
{
    register_all_transformations();
    return TR::lookup(type).param_specs();
}

}  // namespace regpoly::core
