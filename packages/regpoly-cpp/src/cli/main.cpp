// regpoly-cli — standalone C++ command-line driver.
//
// Phase 1 ships an empty stub that prints version + usage and exits.
// Phase 4 fills in the real subcommands (search, display tested
// generator, publish to library) so a C++-only user has full feature
// parity with the Python `regpoly` CLI.

#include <cstring>
#include <iostream>

namespace {

constexpr const char* kVersion = "regpoly-cli 2.0.0 (Phase 1 stub)";
constexpr const char* kUsage =
    "Usage: regpoly-cli <command> [options]\n"
    "\n"
    "Commands (Phase 4):\n"
    "  search   <config.yaml>   run an equidistribution / full-period search\n"
    "  show     <result.yaml>   display a tested-generator result\n"
    "  publish  <result.yaml>   publish a tested generator to the catalog\n"
    "\n"
    "Options:\n"
    "  -h, --help        show this message and exit\n"
    "  -V, --version     print version and exit\n"
    "\n"
    "Phase 1 status: command dispatch is not yet implemented. The C++\n"
    "core (regpoly_core static library + public API headers under\n"
    "<regpoly/regpoly.h>) is available and the SearchDriver family will\n"
    "land in Phase 2.\n";

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << kUsage;
        return 0;
    }
    const char* arg = argv[1];
    if (std::strcmp(arg, "-h") == 0 || std::strcmp(arg, "--help") == 0) {
        std::cout << kUsage;
        return 0;
    }
    if (std::strcmp(arg, "-V") == 0 || std::strcmp(arg, "--version") == 0) {
        std::cout << kVersion << "\n";
        return 0;
    }
    std::cerr << "regpoly-cli: command '" << arg
              << "' is not yet implemented (Phase 4).\n";
    std::cerr << "Run `regpoly-cli --help` for the available commands.\n";
    return 2;
}
