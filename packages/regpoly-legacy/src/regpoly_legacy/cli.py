# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""regpoly-legacy console script.

Subcommands:
    info  <file.dat>                    parse + print generator specs
    trans <file.dat>                    parse + print transformation specs
    seek  <nb_comp> <test_file> <gen_file1> [<gen_file2> ...]
                                        run a search driven by legacy
                                        positional inputs (replaces the
                                        old ``regpoly nb_comp ...`` mode)
"""

from __future__ import annotations

import argparse
import sys


def _cmd_info(args: argparse.Namespace) -> int:
    from regpoly_legacy.reader import parse_generator_specs
    specs = parse_generator_specs(args.path, args.L)
    print(f"File: {args.path}")
    print(f"Generators: {len(specs)}")
    for i, (family, params) in enumerate(specs):
        bits = ", ".join(f"{k}={v}" for k, v in params.items())
        print(f"  [{i}] family={family}  {bits}")
    return 0


def _cmd_trans(args: argparse.Namespace) -> int:
    from regpoly_legacy.reader import parse_transformation_specs
    specs, mk_opt = parse_transformation_specs(args.path)
    print(f"File: {args.path}")
    print(f"Transformations: {len(specs)}  mk_opt={mk_opt}")
    for i, (trans_type, params) in enumerate(specs):
        bits = ", ".join(f"{k}={v}" for k, v in params.items())
        print(f"  [{i}] type={trans_type}  {bits}")
    return 0


def _cmd_seek(args: argparse.Namespace) -> int:
    from regpoly_legacy.seek_factory import seek_from_legacy
    if len(args.gen_data_files) != args.nb_comp:
        print(
            f"Expected {args.nb_comp} generator file(s), got "
            f"{len(args.gen_data_files)}",
            file=sys.stderr,
        )
        return 1
    seek_from_legacy(args.nb_comp, args.test_file, args.gen_data_files).run()
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="regpoly-legacy",
        description="Read pre-v2 .dat parameter files into regpoly objects.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    info = sub.add_parser("info", help="print generator specs from a .dat")
    info.add_argument("path")
    info.add_argument("--L", type=int, default=32,
                      help="output resolution in bits (default 32)")
    info.set_defaults(func=_cmd_info)

    trans = sub.add_parser("trans", help="print transformation specs from a .dat")
    trans.add_argument("path")
    trans.set_defaults(func=_cmd_trans)

    seek = sub.add_parser("seek", help="run Seek from legacy positional inputs")
    seek.add_argument("nb_comp", type=int)
    seek.add_argument("test_file")
    seek.add_argument("gen_data_files", nargs="+")
    seek.set_defaults(func=_cmd_seek)

    args = p.parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
