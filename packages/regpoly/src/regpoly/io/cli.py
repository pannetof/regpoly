# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
regpoly.io.cli — Command-line entry point.

Usage:
    regpoly equidist/config.yaml                                (equidistribution search)
    regpoly fullperiodsearch/Family.desc.yaml                   (full-period search)
    regpoly testedgenerators/Family/params.equidist.0001.yaml   (re-test a saved generator)

The YAML format is auto-detected:
  - "search.family" key → full-period search
  - "generator" or "components" + "results" → tested generator (display results)
  - "components" + "tests" → equidistribution search (Seek)

The pre-v2 positional-argument form (``regpoly nb_comp test_file
gen_file1 [...]``) has moved to the optional ``regpoly-legacy`` add-on
— install it and run ``regpoly-legacy seek ...`` instead.
"""

import sys


def main() -> None:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <config.yaml>")
        sys.exit(1)

    arg = sys.argv[1]

    if not (arg.endswith(".yaml") or arg.endswith(".yml")):
        print(
            f"Error: expected a .yaml / .yml config file, got {arg!r}.\n"
            f"The legacy positional invocation "
            f"(nb_comp test_file gen_file1 ...) moved to the regpoly-legacy "
            f"add-on; install it and use `regpoly-legacy seek ...` instead.",
            file=sys.stderr,
        )
        sys.exit(1)

    import yaml
    with open(arg) as f:
        config = yaml.safe_load(f)

    search = config.get("search", {})

    if "family" in search:
        # Full-period search
        from regpoly.search.search_primitive import PrimitiveSearch
        PrimitiveSearch.from_yaml(arg).run()

    elif "generator" in config or (
        "components" in config and "results" in config
    ):
        # Tested generator file — display its contents
        from regpoly.io.tested_generator import load_tested_generator
        comb, results = load_tested_generator(arg)
        print(f"Loaded tested generator from {arg}")
        print(f"  k_g = {comb.k_g}, L = {comb.L}, J = {comb.J}")
        if results:
            for test_name, res in results.items():
                print(f"  {test_name}: {res}")

    else:
        # Equidistribution search
        from regpoly.search.seek import Seek
        Seek.from_yaml(arg).run()


if __name__ == "__main__":
    main()
