#pragma once

// Phase 3.3: parsers for the old text-based ".dat" parameter files.
//
// Mirrors packages/regpoly/src/regpoly/io/legacy_reader.py 1:1.
// The file's first whitespace-separated token is the dispatch tag
// (e.g. "polylcg", "taus", "tgfsr", "MT", "genf2w", "carry",
// "marsaxorshift", "matsumoto"). The transformation reader is
// tag-less and is dispatched by the per-line type ("permut",
// "tempMK", "tempMK2", "tempMKopt", "tempMK2opt").
//
// Two flavours of API are provided:
//   * read_*_specs       — return data-only descriptors (family +
//                          Params), suitable for the pybind11 layer
//                          that wants to round-trip via the Python
//                          Generator.create() / Transformation.create()
//                          factories (which preserve the Python-side
//                          _params dict for serialisation).
//   * read_generators /
//     read_transformations — build the C++ objects directly via
//                          create_generator() / create_transformation()
//                          from <regpoly/factory.h>.

#include "generator.h"
#include "transformation.h"
#include "params.h"

#include <memory>
#include <string>
#include <vector>

namespace regpoly_legacy {

struct LegacyGeneratorSpec {
    std::string family;   ///< canonical C++ family name (e.g. "PolyLCG",
                          ///< "F2wLFSRGen", "WELLGen"); see translation
                          ///< table at the top of legacy_reader.cpp.
    Params params;        ///< already-built Params object (coeff/nocoeff
                          ///< flattening for genf2w/carry already applied).
};

struct LegacyTransformationSpec {
    std::string trans_type;   ///< "permut", "tempMK", "tempMK2"
    Params params;
};

struct LegacyTransformationsResult {
    std::vector<LegacyTransformationSpec> specs;
    bool mk_opt = false;      ///< true iff any tempMKopt / tempMK2opt was seen.
};

struct LegacyTransformationsBuilt {
    std::vector<std::unique_ptr<Transformation>> transformations;
    bool mk_opt = false;
};

// ── Generator dispatch ────────────────────────────────────────────────

// Parse a legacy generator file and return one spec per generator.
// The file's first token names the family ("polylcg", "taus[2]",
// "tgfsr", "MT", "genf2w", "carry", "marsaxorshift", "matsumoto").
// Throws std::runtime_error on parse / I/O errors (the message
// includes the filename).
std::vector<LegacyGeneratorSpec>
read_generator_specs(const std::string& filename, int L);

// Same as read_generator_specs but eagerly builds the generators via
// create_generator() (so a C++-only consumer never touches Params).
std::vector<std::unique_ptr<Generator>>
read_generators(const std::string& filename, int L);

// ── Transformation dispatch ───────────────────────────────────────────

// Parse a legacy transformations file. The file has no top-level tag —
// the first non-blank line is the count, then each subsequent line
// starts with the per-transformation type word.
LegacyTransformationsResult
read_transformation_specs(const std::string& filename);

LegacyTransformationsBuilt
read_transformations(const std::string& filename);

}  // namespace regpoly_legacy
