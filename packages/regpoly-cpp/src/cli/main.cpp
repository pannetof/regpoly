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
// Phase 4.2+ will add `search`, `show <result.yaml>`, and `publish`.

#include "catalog.h"
#include "legacy_reader.h"
#include "params.h"

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kVersion = "regpoly-cli 2.0.0 (Phase 4.1)";
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

    std::cerr << "regpoly-cli: unknown command '" << cmd << "'.\n"
              << "Run `regpoly-cli --help`.\n";
    return 2;
}
