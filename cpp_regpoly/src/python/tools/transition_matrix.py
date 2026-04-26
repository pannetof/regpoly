"""
transition_matrix.py — Display the transition matrix of every generator
in a data file.

Usage
-----
    python -m regpoly.tools.transition_matrix <generator_file>

The file type is detected from its first token (e.g. "polylcg", "taus",
"taus2", "tgfsr").  Each generator is loaded and its transition matrix
is computed and printed.
"""

import sys

from regpoly.io.legacy_reader import LegacyReader
from regpoly.core.component import Component


def main() -> None:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <generator_file>")
        sys.exit(1)

    path = sys.argv[1]

    # Detect generator type for display
    with open(path) as f:
        gen_type = f.readline().split()[0]

    # Map tag to class name (matching C test_transition_matrix_generator)
    _class_names = {
        "polylcg": "PolyLCG",
        "taus": "Tausworthe", "taus2": "Tausworthe",
        "tgfsr": "TGFSR",
        "MT": "MersenneTwister",
        "genf2w": "GenF2w",
    }
    class_name = _class_names.get(gen_type, gen_type)

    L = 4096
    gen_list = LegacyReader.read_generators(path, L)

    comp = Component()
    for gen in gen_list:
        comp.add_gen(gen)

    print(f"File   : {path}")
    print(f"Type   : {gen_type}  ({class_name})")
    print(f"Generators: {len(comp)}")
    print()

    for i, gen in enumerate(comp.gens):
        gen.display()
        try:
            A = gen.transition_matrix()
            A.display()
        except Exception as exc:
            print(f"  [{exc}]")
        print()


if __name__ == "__main__":
    main()
