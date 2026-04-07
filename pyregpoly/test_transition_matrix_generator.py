"""
test_transition_matrix_generator.py — Display the transition matrix of every
generator in a data file.

Usage
-----
    python test_transition_matrix_generator.py <generator_file>

The file type is detected from its first token (e.g. "polylcg", "taus",
"taus2").  Each generator is loaded into a Component and its transition
matrix is computed and printed in turn.
"""

import sys

from bitvect import BitVect
from combinaison import Combinaison
from component import Component


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
            # Temporarily set L >= k so transition matrix can be computed
            saved_L = gen.L
            if gen.L < gen.k:
                gen.L = gen.k
                gen.gen_state = BitVect.zeros(gen.L)
            A = gen.transition_matrix()
            gen.L = saved_L
            gen.gen_state = BitVect.zeros(saved_L)
            A.display()
        except ValueError as exc:
            print(f"  [{exc}]")
        print()


if __name__ == "__main__":
    main()
