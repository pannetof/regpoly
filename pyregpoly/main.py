"""
main.py — Entry point for the combined-generator search.

Usage:
    python main.py nb_comp test_file gen_file1 [gen_file2 ...]

Mirrors the C binary interface: ./POL nb_comp test_file gen_file1 [gen_file2 ...]
"""

import sys
from seek import Seek

if len(sys.argv) < 4:
    print(f"Usage: {sys.argv[0]} nb_comp test_file gen_file1 [gen_file2 ...]")
    sys.exit(1)

nb_comp        = int(sys.argv[1])
test_file      = sys.argv[2]
gen_data_files = sys.argv[3:]

if len(gen_data_files) != nb_comp:
    print(f"Expected {nb_comp} generator file(s), got {len(gen_data_files)}")
    sys.exit(1)

Seek(nb_comp, test_file, gen_data_files).run()
