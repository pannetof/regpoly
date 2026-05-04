// regpoly-cli — standalone C++ command-line driver.
//
// Phase 4.1 ships two foundational subcommands that exercise the
// already-ported C++ machinery without requiring a YAML search-config
// loader:
//
//   regpoly-cli catalog list [--library DIR]
//                            print every paper id + display title.
//   regpoly-cli catalog show PAPER_ID [--library DIR]
//                            print one paper's full record (authors,
//                            citation, generator list).
//   regpoly-cli catalog gen  GEN_ID  [--library DIR]
//                            print one generator's record + paper.
//   regpoly-cli legacy-info  FILE.dat [-L N]
//                            parse a legacy generator file via the C++
//                            legacy reader and print one line per
//                            parsed generator (family + params).
//   regpoly-cli legacy-trans FILE.dat
//                            same for a legacy transformations file.
//
// Phase 4.2 adds `search FILE.yaml` — load a seek-style YAML config
// and run the equidistribution search loop via the existing C++
// drivers. Phase 4.3 adds `show <result.yaml>` — display a
// tested-generator file (single- or multi-component shape) including
// its tempering chain and equidistribution / collision-free /
// tuplets results. The `publish` subcommand (catalog write) lands
// later.

#include "catalog.h"
#include "combination.h"
#include "legacy_reader.h"
#include "params.h"
#include "seek_config.h"
#include "seek_search.h"

#include <yaml-cpp/yaml.h>

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kVersion = "regpoly-cli 2.0.0 (Phase 4.3)";
constexpr const char* kUsage =
    "Usage: regpoly-cli <command> [options]\n"
    "\n"
    "Commands:\n"
    "  catalog list [--library DIR]\n"
    "      List every paper id + display title.\n"
    "  catalog show PAPER_ID [--library DIR]\n"
    "      Print one paper's full record.\n"
    "  catalog gen GEN_ID [--library DIR]\n"
    "      Print one generator's record + its paper.\n"
    "  legacy-info FILE.dat [-L N]\n"
    "      Parse a legacy .dat generator file (default L=32).\n"
    "  legacy-trans FILE.dat\n"
    "      Parse a legacy .dat transformations file.\n"
    "  search FILE.yaml\n"
    "      Load a seek-style YAML search config and run the\n"
    "      equidistribution search loop. Output format does not\n"
    "      mirror `uv run regpoly`; use that command for the\n"
    "      Python-side display.\n"
    "  show FILE.yaml\n"
    "      Display a tested-generator YAML (single- or multi-\n"
    "      component): components + tempering chain + results.\n"
    "\n"
    "Options:\n"
    "  -h, --help        show this message and exit\n"
    "  -V, --version     print version and exit\n"
    "\n"
    "Default catalog dir is the docs/library/ next to this binary,\n"
    "found by walking up from the current directory.\n";

// Locate the docs/library directory by walking up from cwd looking
// for a sibling docs/. Returns empty string if not found within 8
// levels.
std::string find_default_library_dir() {
    fs::path here = fs::current_path();
    for (int up = 0; up < 8; ++up) {
        fs::path cand = here / "docs" / "library";
        if (fs::is_directory(cand)) return cand.string();
        if (here == here.parent_path()) break;
        here = here.parent_path();
    }
    return {};
}

// Parse `--library DIR` or `-l DIR` if present in args; returns the
// directory and removes the args from the vector.
std::string consume_library_flag(std::vector<std::string>& args) {
    for (size_t i = 0; i < args.size(); ++i) {
        if ((args[i] == "--library" || args[i] == "-l")
            && i + 1 < args.size()) {
            std::string dir = args[i + 1];
            args.erase(args.begin() + i, args.begin() + i + 2);
            return dir;
        }
    }
    return find_default_library_dir();
}

// Parse `-L N` if present; returns the int (default 32) and removes
// from the args.
int consume_L_flag(std::vector<std::string>& args, int def = 32) {
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "-L" && i + 1 < args.size()) {
            int v = std::atoi(args[i + 1].c_str());
            args.erase(args.begin() + i, args.begin() + i + 2);
            return v;
        }
    }
    return def;
}

std::string params_to_one_line(const Params& p) {
    std::ostringstream os;
    bool first = true;
    auto sep = [&]() { if (!first) os << ", "; first = false; };
    for (const auto& kv : p.ints())     { sep(); os << kv.first << "=" << kv.second; }
    for (const auto& kv : p.bools())    { sep(); os << kv.first << "=" << (kv.second ? "true" : "false"); }
    for (const auto& kv : p.strings())  { sep(); os << kv.first << "='" << kv.second << "'"; }
    for (const auto& kv : p.int_vecs()) {
        sep(); os << kv.first << "=[";
        for (size_t i = 0; i < kv.second.size(); ++i) {
            if (i) os << ",";
            os << kv.second[i];
        }
        os << "]";
    }
    for (const auto& kv : p.uint_vecs()) {
        sep(); os << kv.first << "=[";
        for (size_t i = 0; i < kv.second.size(); ++i) {
            if (i) os << ",";
            os << "0x" << std::hex << kv.second[i] << std::dec;
        }
        os << "]";
    }
    return os.str();
}

int cmd_catalog(std::vector<std::string> args) {
    if (args.empty()) {
        std::cerr << "regpoly-cli: catalog requires a sub-action "
                  << "(list | show | gen)\n";
        return 2;
    }
    std::string action = args[0];
    args.erase(args.begin());

    std::string library_dir = consume_library_flag(args);
    if (library_dir.empty()) {
        std::cerr << "regpoly-cli: could not locate docs/library/. "
                  << "Pass --library DIR.\n";
        return 2;
    }

    regpoly_catalog::Catalog cat(library_dir);
    cat.load();

    if (action == "list") {
        regpoly_catalog::Catalog::PapersFilter f;
        f.include_invalid = true;
        auto papers = cat.papers(f);
        for (const auto& p : papers) {
            std::cout << p.id << " — " << p.display();
            if (!p.valid()) std::cout << " (INVALID: " << p.errors.size()
                                       << " error" << (p.errors.size() == 1 ? "" : "s")
                                       << ")";
            std::cout << "\n";
        }
        return 0;
    }

    if (action == "show") {
        if (args.empty()) {
            std::cerr << "regpoly-cli: catalog show requires PAPER_ID\n";
            return 2;
        }
        auto p = cat.paper(args[0]);
        if (!p.has_value()) {
            std::cerr << "regpoly-cli: no such paper: " << args[0] << "\n";
            return 1;
        }
        std::cout << "id:        " << p->id << "\n";
        std::cout << "display:   " << p->display() << "\n";
        std::cout << "year:      " << p->year << "\n";
        std::cout << "title:     " << p->title << "\n";
        std::cout << "venue:     " << p->venue << "\n";
        if (!p->doi.empty())   std::cout << "doi:       " << p->doi << "\n";
        if (!p->bibkey.empty())std::cout << "bibkey:    " << p->bibkey << "\n";
        std::cout << "starred:   " << (p->starred ? "yes" : "no") << "\n";
        std::cout << "valid:     " << (p->valid() ? "yes" : "no") << "\n";
        for (const auto& e : p->errors) std::cout << "  ! " << e << "\n";
        std::cout << "citation:  " << p->acmtrans_citation() << "\n";
        std::cout << "generators (" << p->generators.size() << "):\n";
        for (const auto& g : p->generators) {
            std::cout << "  - " << g.id << " [" << g.family << "]"
                      << "  L=" << g.Lmax;
            if (g.starred) std::cout << "  *";
            if (!g.valid()) std::cout << "  INVALID";
            std::cout << "\n";
        }
        return 0;
    }

    if (action == "gen") {
        if (args.empty()) {
            std::cerr << "regpoly-cli: catalog gen requires GEN_ID\n";
            return 2;
        }
        auto loc = cat.generator(args[0]);
        if (!loc.has_value()) {
            std::cerr << "regpoly-cli: no such generator: " << args[0] << "\n";
            return 1;
        }
        const auto& [paper, gen] = *loc;
        std::cout << "generator: " << gen.id << "\n";
        std::cout << "display:   " << gen.display << "\n";
        std::cout << "family:    " << gen.family << "\n";
        std::cout << "Lmax:      " << gen.Lmax << "\n";
        std::cout << "target:    " << gen.target << "\n";
        std::cout << "combined:  " << (gen.combined ? "yes" : "no") << "\n";
        std::cout << "components: " << gen.components.size() << "\n";
        std::cout << "paper:     " << paper.id << " — " << paper.display() << "\n";
        return 0;
    }

    std::cerr << "regpoly-cli: unknown catalog sub-action '"
              << action << "'\n";
    return 2;
}

int cmd_legacy_info(std::vector<std::string> args) {
    if (args.empty()) {
        std::cerr << "regpoly-cli: legacy-info requires FILE.dat\n";
        return 2;
    }
    int L = consume_L_flag(args, 32);
    if (args.empty()) {
        std::cerr << "regpoly-cli: legacy-info requires FILE.dat\n";
        return 2;
    }
    std::string path = args[0];
    try {
        auto specs = regpoly_legacy::read_generator_specs(path, L);
        std::cout << "file:      " << path << "\n";
        std::cout << "L:         " << L << "\n";
        std::cout << "generators: " << specs.size() << "\n";
        for (size_t i = 0; i < specs.size(); ++i) {
            const auto& s = specs[i];
            std::cout << "  [" << i << "] " << s.family << "  "
                      << params_to_one_line(s.params) << "\n";
        }
    } catch (const std::exception& exc) {
        std::cerr << "regpoly-cli: " << exc.what() << "\n";
        return 1;
    }
    return 0;
}

const char* test_kind_name(SeekTestKind k) {
    switch (k) {
        case SeekTestKind::EquidistributionMatricial:        return "equidist[matricial]";
        case SeekTestKind::EquidistributionLattice:          return "equidist[lattice]";
        case SeekTestKind::EquidistributionHarase:           return "equidist[harase]";
        case SeekTestKind::EquidistributionNotPrimitive:     return "equidist[notprimitive]";
        case SeekTestKind::EquidistributionSimdNotPrimitive: return "equidist[simd_notprimitive]";
        case SeekTestKind::EquidistributionNothing:          return "equidist[nothing]";
        case SeekTestKind::CollisionFree:                    return "collision_free";
        case SeekTestKind::Tuplets:                          return "tuplets";
    }
    return "?";
}

int cmd_search(std::vector<std::string> args) {
    if (args.empty()) {
        std::cerr << "regpoly-cli: search requires FILE.yaml\n";
        return 2;
    }
    const std::string yaml_path = args[0];

    regpoly_yaml_config::SeekConfig cfg;
    try {
        cfg = regpoly_yaml_config::load_seek_config(yaml_path);
    } catch (const std::exception& exc) {
        std::cerr << "regpoly-cli: " << exc.what() << "\n";
        return 1;
    }

    regpoly_yaml_config::BuiltSearch built;
    try {
        built = regpoly_yaml_config::build_search(cfg);
    } catch (const std::exception& exc) {
        std::cerr << "regpoly-cli: " << exc.what() << "\n";
        return 1;
    }

    // Header — brief, machine-readable-ish.
    std::cout << "regpoly-cli search\n";
    std::cout << "  config:    " << yaml_path << "\n";
    std::cout << "  seed:      [" << cfg.seed1 << ", " << cfg.seed2 << "]\n";
    std::cout << "  Lmax:      " << cfg.Lmax << "\n";
    std::cout << "  J:         " << cfg.components.size() << "\n";
    std::cout << "  nbtries:   " << cfg.nbtries << "\n";
    std::cout << "  tests:     ";
    for (size_t i = 0; i < cfg.tests.size(); ++i) {
        if (i) std::cout << ", ";
        std::cout << test_kind_name(cfg.tests[i].kind);
    }
    std::cout << "\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout.flush();

    int64_t selection_count = 0;
    auto on_iter = [&](Combination& comb, const SeekIterResult& r) {
        ++selection_count;
        std::cout << "  [" << std::setw(6) << selection_count << "] ";
        std::cout << "k_g=" << comb.k_g() << " L=" << comb.L();
        if (r.me_ran) {
            std::cout << "  me_se=" << r.me_se;
            std::cout << " (verified=" << (r.me_verified ? "yes" : "no")
                      << (r.me_is_me ? ", ME" : "") << ")";
        }
        if (r.cf_ran) {
            std::cout << "  cf_secf=" << r.cf_secf
                      << " (verified=" << (r.cf_verified ? "yes" : "no") << ")";
        }
        if (r.tup_ran) {
            std::cout << "  tup_first_max=" << r.tup_firstpart_max
                      << " tup_first_sum=" << r.tup_firstpart_sum;
        }
        std::cout << "\n";
        std::cout.flush();
    };

    auto on_progress = [&](const SeekProgress& p) {
        std::cout << "  ... progress: combos=" << p.nbgen
                  << " selected=" << p.nb_select
                  << " ME=" << p.nb_me
                  << " elapsed=" << std::fixed << std::setprecision(2)
                  << p.elapsed_seconds << "s\n";
        std::cout.flush();
    };

    SeekResult result;
    try {
        result = run_seek_search(*built.combination, cfg.tests,
                                 cfg.nbtries, /*progress_interval=*/1000,
                                 /*on_prep=*/nullptr,
                                 on_iter, on_progress);
    } catch (const std::exception& exc) {
        std::cerr << "regpoly-cli: search failed: " << exc.what() << "\n";
        return 1;
    }

    std::cout << std::string(60, '=') << "\n";
    std::cout << "Summary:\n";
    std::cout << "  combos:    " << result.nbgen << "\n";
    std::cout << "  selected:  " << result.nb_select << "\n";
    std::cout << "  ME (full): " << result.nb_me << "\n";
    std::cout << "  elapsed:   " << std::fixed << std::setprecision(3)
              << result.elapsed_seconds << "s\n";
    std::cout.flush();
    return 0;
}

// Render one YAML scalar value compactly. Used by cmd_show; the
// tested-generator schema's leaf values are scalars (ints, hex
// strings, bools, floats) or short int sequences.
std::string scalar_to_str(const YAML::Node& n) {
    if (!n.IsScalar()) return "<non-scalar>";
    return n.as<std::string>();
}

void print_kv_block(const YAML::Node& m, const std::string& prefix) {
    if (!m || !m.IsMap()) return;
    for (auto kv : m) {
        auto k = kv.first.as<std::string>();
        const auto& v = kv.second;
        if (v.IsScalar()) {
            std::cout << prefix << k << ": " << scalar_to_str(v) << "\n";
        } else if (v.IsSequence()) {
            std::cout << prefix << k << ": [";
            for (size_t i = 0; i < v.size(); ++i) {
                if (i) std::cout << ", ";
                if (v[i].IsScalar()) std::cout << scalar_to_str(v[i]);
                else std::cout << "?";
            }
            std::cout << "]\n";
        } else if (v.IsMap()) {
            std::cout << prefix << k << ":\n";
            print_kv_block(v, prefix + "  ");
        }
    }
}

void print_component(const YAML::Node& gen, const YAML::Node& tempering,
                     int idx) {
    std::cout << "  component " << idx << ":\n";
    if (gen && gen.IsMap()) {
        std::cout << "    generator:\n";
        print_kv_block(gen, "      ");
    }
    if (tempering && tempering.IsSequence() && tempering.size() > 0) {
        std::cout << "    tempering:\n";
        for (size_t i = 0; i < tempering.size(); ++i) {
            std::cout << "      - ";
            const auto& step = tempering[i];
            if (step["type"]) {
                std::cout << "type: " << step["type"].as<std::string>() << "\n";
            } else {
                std::cout << "(no type)\n";
            }
            for (auto kv : step) {
                auto k = kv.first.as<std::string>();
                if (k == "type") continue;
                std::cout << "        " << k << ": "
                          << scalar_to_str(kv.second) << "\n";
            }
        }
    }
}

void print_results(const YAML::Node& results) {
    if (!results || !results.IsMap()) return;
    std::cout << "  results:\n";
    for (auto rkv : results) {
        std::cout << "    " << rkv.first.as<std::string>() << ":\n";
        print_kv_block(rkv.second, "      ");
    }
}

int cmd_show(std::vector<std::string> args) {
    if (args.empty()) {
        std::cerr << "regpoly-cli: show requires FILE.yaml\n";
        return 2;
    }
    std::string path = args[0];
    YAML::Node doc;
    try {
        doc = YAML::LoadFile(path);
    } catch (const std::exception& exc) {
        std::cerr << "regpoly-cli: failed to load " << path
                  << ": " << exc.what() << "\n";
        return 1;
    }
    if (!doc || !doc.IsMap()) {
        std::cerr << "regpoly-cli: " << path
                  << ": top-level YAML must be a mapping\n";
        return 1;
    }

    std::cout << "tested generator: " << path << "\n";

    // Two top-level shapes: single-component (`generator` + `tempering`)
    // or multi-component (`components: [{generator, tempering}, ...]`).
    if (doc["components"] && doc["components"].IsSequence()) {
        std::cout << "  J: " << doc["components"].size() << "\n";
        int idx = 0;
        for (auto comp : doc["components"]) {
            print_component(comp["generator"], comp["tempering"], idx++);
        }
    } else if (doc["generator"]) {
        std::cout << "  J: 1\n";
        print_component(doc["generator"], doc["tempering"], 0);
    } else {
        std::cerr << "regpoly-cli: " << path
                  << ": missing `generator` or `components` key\n";
        return 1;
    }

    if (doc["results"]) print_results(doc["results"]);
    return 0;
}

int cmd_legacy_trans(std::vector<std::string> args) {
    if (args.empty()) {
        std::cerr << "regpoly-cli: legacy-trans requires FILE.dat\n";
        return 2;
    }
    std::string path = args[0];
    try {
        auto r = regpoly_legacy::read_transformation_specs(path);
        std::cout << "file:      " << path << "\n";
        std::cout << "transformations: " << r.specs.size() << "\n";
        std::cout << "mk_opt:    " << (r.mk_opt ? "yes" : "no") << "\n";
        for (size_t i = 0; i < r.specs.size(); ++i) {
            const auto& s = r.specs[i];
            std::cout << "  [" << i << "] " << s.trans_type << "  "
                      << params_to_one_line(s.params) << "\n";
        }
    } catch (const std::exception& exc) {
        std::cerr << "regpoly-cli: " << exc.what() << "\n";
        return 1;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << kUsage;
        return 0;
    }
    std::string cmd = argv[1];
    if (cmd == "-h" || cmd == "--help")    { std::cout << kUsage;            return 0; }
    if (cmd == "-V" || cmd == "--version") { std::cout << kVersion << "\n";  return 0; }

    std::vector<std::string> rest;
    for (int i = 2; i < argc; ++i) rest.emplace_back(argv[i]);

    if (cmd == "catalog")      return cmd_catalog(std::move(rest));
    if (cmd == "legacy-info")  return cmd_legacy_info(std::move(rest));
    if (cmd == "legacy-trans") return cmd_legacy_trans(std::move(rest));
    if (cmd == "search")       return cmd_search(std::move(rest));
    if (cmd == "show")         return cmd_show(std::move(rest));

    std::cerr << "regpoly-cli: unknown command '" << cmd << "'.\n"
              << "Run `regpoly-cli --help`.\n";
    return 2;
}
