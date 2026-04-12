"""
regpoly.cli — Command-line entry point.

Usage:
    regpoly search.config.yaml                          (equidistribution search)
    regpoly search_primitive.yaml                       (full-period search)
    regpoly nb_comp test_file gen_file1 [gen_file2 ...]  (legacy format)

The YAML format is auto-detected:
  - If the file has a top-level "search.family" key → primitive search
  - If the file has a "components" key → equidistribution search (Seek)
"""

import sys


def main() -> None:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <config.yaml>")
        print(f"       {sys.argv[0]} nb_comp test_file gen_file1 [gen_file2 ...]")
        sys.exit(1)

    arg = sys.argv[1]

    if arg.endswith(".yaml") or arg.endswith(".yml"):
        import yaml
        with open(arg) as f:
            config = yaml.safe_load(f)

        search = config.get("search", {})
        if "family" in search:
            from regpoly.search_primitive import PrimitiveSearch
            PrimitiveSearch.from_yaml(arg).run()
        else:
            from regpoly.seek import Seek
            Seek.from_yaml(arg).run()
    else:
        from regpoly.seek import Seek
        if len(sys.argv) < 4:
            print(f"Usage: {sys.argv[0]} nb_comp test_file gen_file1 [gen_file2 ...]")
            sys.exit(1)
        nb_comp = int(sys.argv[1])
        test_file = sys.argv[2]
        gen_data_files = sys.argv[3:]
        if len(gen_data_files) != nb_comp:
            print(f"Expected {nb_comp} generator file(s), got {len(gen_data_files)}")
            sys.exit(1)
        Seek.from_legacy(nb_comp, test_file, gen_data_files).run()


if __name__ == "__main__":
    main()
