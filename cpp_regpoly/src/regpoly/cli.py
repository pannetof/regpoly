"""
regpoly.cli — Command-line entry point for the combined-generator search.

Usage:
    regpoly search.config.yaml                          (YAML format)
    regpoly nb_comp test_file gen_file1 [gen_file2 ...]  (legacy format)
"""

import sys

from regpoly.seek import Seek


def main() -> None:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <search.yaml>")
        print(f"       {sys.argv[0]} nb_comp test_file gen_file1 [gen_file2 ...]")
        sys.exit(1)

    # Detect format: if first arg is a .yaml file, use YAML path
    if sys.argv[1].endswith(".yaml") or sys.argv[1].endswith(".yml"):
        Seek.from_yaml(sys.argv[1]).run()
    else:
        # Legacy format
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
