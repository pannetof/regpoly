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

from regpoly.combinaison import Combinaison
from regpoly.component import Component


def main() -> None:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <generator_file>")
        sys.exit(1)

    path = sys.argv[1]

    # Detect generator type and resolve class
    with open(path) as f:
        gen_type = f.readline().split()[0]

    dispatch = Combinaison._get_dispatch()
    gen_cls = dispatch.get(gen_type)
    if gen_cls is None:
        print(f"Unknown generator type '{gen_type}' in {path}")
        sys.exit(1)

    # Use a large L so every generator fits
    L = 4096
    gen_list = gen_cls.CreateListFromFile(path, L)

    comp = Component()
    for gen in gen_list:
        comp.add_gen(gen)

    print(f"File   : {path}")
    print(f"Type   : {gen_type}  ({gen_cls.__name__})")
    print(f"Generators: {len(comp)}")
    print()

    for i, gen in enumerate(comp.gens):
        gen.display()
        try:
            At = gen.transition_matrix(transposed=True)
            At.display_transposed()
        except Exception as exc:
            print(f"  [{exc}]")
        print()


if __name__ == "__main__":
    main()
