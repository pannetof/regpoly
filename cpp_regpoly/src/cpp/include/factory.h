#pragma once
#include "generateur.h"
#include "transformation.h"
#include "params.h"
#include "param_spec.h"
#include "gen_enumerator.h"
#include <memory>
#include <string>
#include <vector>

// Forward declaration for pybind11 module type
namespace pybind11 { class module_; }

std::unique_ptr<Generateur> create_generator(
    const std::string& family, const Params& params, int L);

std::unique_ptr<Transformation> create_transformation(
    const std::string& type, const Params& params);

std::vector<ParamSpec> get_gen_param_specs(const std::string& family);
std::vector<ParamSpec> get_trans_param_specs(const std::string& type);

// Register all generator subclass types with pybind11
// (defined in factory.cpp where all headers are available)
void register_generator_types(pybind11::module_& m);
