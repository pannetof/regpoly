#!/usr/bin/env python3
"""
Regenerate packages/regpoly-cpp/src/algebra/primitive_factors_data.cpp
from the canonical JSON factor table at:

    packages/regpoly/src/regpoly/data/primitive_factors.json

Run from anywhere — paths are computed relative to this script — or as

    python -m regpoly.tools.gen_primitive_factors_data

The factor table holds the prime factors of Phi_n(2) for n up to a
useful bound (covers all generator state sizes we care about). It is
the sync point between the JSON (read by Python tests + hand-tuned by
humans, produced upstream by ``regpoly.tools.generate_factors``) and
the C++ data table that compiles into ``_regpoly_cpp.so``. The
generated .cpp file is committed to the repo so a fresh build does
not depend on running this script.
"""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    # Path layout: packages/regpoly/src/regpoly/tools/THIS_FILE.py
    #   parents[0] = tools
    #   parents[1] = regpoly (inner)
    #   parents[2] = src
    #   parents[3] = regpoly (outer, the package dir)
    #   parents[4] = packages
    #   parents[5] = repo root
    repo = Path(__file__).resolve().parents[5]
    src = repo / "packages/regpoly/src/regpoly/data/primitive_factors.json"
    dst = repo / "packages/regpoly-cpp/src/algebra/primitive_factors_data.cpp"

    with open(src) as f:
        data = json.load(f)

    out: list[str] = []
    out.append("// AUTO-GENERATED from "
               "packages/regpoly/src/regpoly/data/primitive_factors.json")
    out.append("// Regenerate via "
               "`python -m regpoly.tools.gen_primitive_factors_data`")
    out.append("// Embeds prime factors of Phi_n(2) for primitivity testing.")
    out.append("")
    out.append('#include "primitivity.h"')
    out.append("")
    out.append("#include <unordered_map>")
    out.append("#include <vector>")
    out.append("#include <string>")
    out.append("")
    out.append("namespace regpoly_internal {")
    out.append("")
    out.append("struct FactorEntry {")
    out.append("    bool complete;")
    out.append("    std::vector<std::string> factors;")
    out.append("};")
    out.append("")
    out.append("static const std::unordered_map<int, FactorEntry>& "
               "factors_table() {")
    out.append("    static const std::unordered_map<int, FactorEntry> t = {")

    for k in sorted(int(k) for k in data):
        entry = data[str(k)]
        facs = ", ".join('"' + str(p) + '"' for p in entry["factors"])
        comp = "true" if entry["complete"] else "false"
        out.append(f"        {{{k}, {{{comp}, {{{facs}}}}}}},")

    out.append("    };")
    out.append("    return t;")
    out.append("}")
    out.append("")
    out.append("const std::vector<std::string>* lookup_factors("
               "int k, bool& complete) {")
    out.append("    const auto& t = factors_table();")
    out.append("    auto it = t.find(k);")
    out.append("    if (it == t.end()) return nullptr;")
    out.append("    complete = it->second.complete;")
    out.append("    return &it->second.factors;")
    out.append("}")
    out.append("")
    out.append("}  // namespace regpoly_internal")
    out.append("")

    with open(dst, "w") as f:
        f.write("\n".join(out))

    print(f"Wrote {dst} ({len(data)} entries)")


if __name__ == "__main__":
    main()
