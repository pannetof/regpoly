// Minimal standalone consumer of the installed regpoly C++ library.
//
// Proves, with no Python present:
//   - <regpoly/regpoly.h> (the umbrella) compiles,
//   - regpoly::regpoly links (so NTL + yaml-cpp resolve via find_dependency),
//   - real out-of-line symbols (Catalog::load → yaml-cpp-backed code) link.

#include <regpoly/regpoly.h>   // umbrella public header
#include <regpoly/catalog.h>   // Catalog API (catalog.cpp references yaml-cpp)

#include <iostream>

int main() {
    regpoly::core::BitVect v(32);
    v.set_bit(0, 1);

    // Touch yaml-cpp-backed code; the directory is intentionally absent, so
    // load() throws — we only need the symbol to link.
    regpoly::library::Catalog cat("/nonexistent-regpoly-catalog");
    try {
        cat.load();
    } catch (const std::exception&) {
        // expected: missing catalog dir
    }

    std::cout << "regpoly find_package consumer OK (bitvect nbits="
              << v.nbits() << ")\n";
    return 0;
}
