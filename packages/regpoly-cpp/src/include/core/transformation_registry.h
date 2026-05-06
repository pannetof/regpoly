// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "param_spec.h"
#include "params.h"
#include "transformation.h"

#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Mirror of GeneratorRegistry for Transformation subclasses.
//
// Adding a new transformation type means:
//   1. Write bar.h / bar.cpp with the Transformation subclass.
//   2. Add ONE line in factory.cpp's register_all_transformations()
//      block: TR::reg("bar", &BarTrans::from_params, &BarTrans::param_specs);
//   3. (No pybind11 binding needed — Transformations are exposed via
//      a generic create_transformation()/get_trans_param_specs()
//      surface, no per-class py::class_.)
//
// The factory entry points create_transformation() and
// get_trans_param_specs() are pure registry lookups.
class TransformationRegistry {
public:
    using FromParamsFn = std::function<
        std::unique_ptr<Transformation>(const Params& params)>;
    using ParamSpecsFn = std::function<std::vector<ParamSpec>()>;

    struct Info {
        std::string canonical_name;
        FromParamsFn from_params;
        ParamSpecsFn param_specs;
    };

    static int reg(const std::string& canonical_name,
                   FromParamsFn from_params,
                   ParamSpecsFn param_specs);

    static const Info& lookup(const std::string& name);
    static const Info* find(const std::string& name);
};
