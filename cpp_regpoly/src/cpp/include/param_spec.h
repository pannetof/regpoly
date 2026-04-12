#pragma once
#include <string>
#include <vector>
#include <cstdint>

struct ParamSpec {
    std::string name;
    std::string type;       // "int", "bool", "int_vec", "uint_vec"
    bool structural;        // true → defines k, user must provide, never random
    bool has_default;
    int64_t default_val;    // meaningful only when has_default && type=="int"|"bool"
    std::string rand_type;  // "bitmask", "range", "poly_exponents",
                            // "bitmask_vec", "none", ""
    std::string rand_args;  // depends on rand_type:
                            //   bitmask:         param name for bit count (e.g. "w")
                            //   range:           "min,max-expr" (e.g. "1,r-1")
                            //   poly_exponents:  param name for degree (e.g. "k")
                            //   bitmask_vec:     "bits_param,length_param" (e.g. "w,nocoeff")
                            //   none/"":         not randomizable
};
