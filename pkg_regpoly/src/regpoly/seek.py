"""
seek.py — Seek class: reads the same files as rechercheguide.c and runs
the combined-generator search, using EquidistributionTest/Results and
CollisionFreeTest/Results.
"""

from __future__ import annotations

import sys
import time

from regpoly.tests.equidistribution_test import (
    EquidistributionTest,
    METHOD_MATRICIAL,
    METHOD_DUALLATTICE,
    METHOD_NOTHING,
)
from regpoly.tests.collision_free_test import CollisionFreeTest
from regpoly.tests.tuplets_test import TupletsTest
from regpoly.tests.tuplets_results import _MAX_TYPE
from regpoly.combinaison import Combinaison
from regpoly.transformation import Transformation

_SEP      = "\n\n" + "+" * 104
_EQ66     = "=" * 66
_EQ66dash = "-" * 66


class Seek:
    """
    Search for combined generators with good equidistribution properties.

    Reads the same file formats as rechercheguide.c:
      - test_file     : search-parameter file (ReadSearch format)
      - gen_data_files: one generator-data file per component (ReadGenDataFiles format)

    Call run() to execute the search and print results to stdout.

    Parameters
    ----------
    nb_comp        : int       — number of components J
    test_file      : str       — path to the search-parameter file
    gen_data_files : list[str] — paths to the J generator-data files
    """

    def __init__(
        self,
        nb_comp: int,
        test_file: str,
        gen_data_files: list,
    ) -> None:
        self.nb_comp        = nb_comp
        self.test_file      = test_file
        self.gen_data_files = list(gen_data_files)

    # -- Public entry point -----------------------------------------------

    def run(self) -> None:
        """Execute the search, printing results to stdout."""
        comb, me_test, tuplets_test, nbtries, mk_opt, Lmax = self._read_search()

        if not comb.reset():
            print("First combined generator not found or invalid Combinaison")
            return

        if mk_opt:
            import warnings
            warnings.warn(
                "MK tempering optimisation (TestMETemperMK) is not implemented; "
                "falling back to standard TestEquid without optimisation.",
                stacklevel=2,
            )

        nbgen = nbsel = nbME = nbCF = 0
        no_try  = 1
        t_start = time.time()

        while True:
            comb.update_trans()
            nbgen += 1

            me_results = me_test.run(comb)

            if me_results.is_presque_me():
                tup_results = tuplets_test.run(comb)
                if tup_results.is_ok():
                    self._display_current_comb(comb)
                    if not me_results.is_me():
                        print("\n  Dimension gaps for every resolution", end="")
                        me_results.display_table(comb, 'l')
                    else:
                        me_results.display()
                        nbME += 1
                    tup_results.display()
                    nbsel += 1
                    print(_SEP)

            sys.stdout.flush()

            if no_try < nbtries:
                no_try += 1
            else:
                try:
                    next(comb)
                except StopIteration:
                    break
                no_try = 1

        elapsed = time.time() - t_start
        self._display_summary(nbgen, nbME, nbCF, nbsel, elapsed)

    # -- File reading -----------------------------------------------------

    def _read_search(self):
        """
        ReadSearch: parse the test file, build all test objects, and return:
          (comb, me_test, tuplets_test, nbtries, mk_opt, Lmax)

        File format (comments after #):
          seed1 seed2          (or -1 for a time-based seed)
          Lmax
          deftrans1 [trans_file1]   (1 = apply tempering to component j)
          ...
          defTransJ [trans_fileJ]
          nbtries
          mse [method]         (mse == -1 means no ME test)
          nb_lines_delta
          jmin jmax ec         (repeated nb_lines_delta times)
          T [d t1 ... td mDD]  (T == 0 means no tuplets test)
        """
        with open(self.test_file, encoding="latin-1") as f:
            lines = [ln.split('#')[0].strip() for ln in f]
        lines = [ln for ln in lines if ln]
        i = 0

        # Seeds
        parts = lines[i].split(); i += 1
        seed1 = int(parts[0])
        seed2 = int(parts[1]) if seed1 != -1 else 0
        seeds = [seed1, seed2]

        # Lmax
        Lmax = int(lines[i]); i += 1

        # Tempering per component
        temperings = []
        mk_opt     = False
        has_tempering = False
        for _ in range(self.nb_comp):
            parts    = lines[i].split(); i += 1
            deftrans = int(parts[0])
            if deftrans:
                trans_list, mk_j = Transformation.CreateFromFiles(parts[1])
                temperings.append(trans_list)
                mk_opt        |= mk_j
                has_tempering  = True
            else:
                temperings.append([])

        # nbtries
        nbtries = int(lines[i]); i += 1
        if not has_tempering:
            nbtries = 1

        # mse and method
        parts   = lines[i].split(); i += 1
        mse_val = int(parts[0])
        if mse_val == -1:
            meverif = False
            mse     = sys.maxsize
            method  = METHOD_NOTHING
        else:
            meverif    = True
            mse        = mse_val
            method_str = parts[1].lower()
            if method_str == "matrix":
                method = METHOD_MATRICIAL
            elif method_str == "lattice":
                method = METHOD_DUALLATTICE
            else:
                raise ValueError(f"Unknown method '{method_str}' in {self.test_file}")

        # Delta lines
        delta       = [10_000_000] * (Lmax + 1)
        delta_lines = []
        nb_lines    = int(lines[i]); i += 1
        for _ in range(nb_lines):
            parts = lines[i].split(); i += 1
            jmin, jmax, ec = int(parts[0]), int(parts[1]), int(parts[2])
            delta_lines.append((jmin, jmax, ec))
            for j in range(jmin, jmax + 1):
                delta[j] = ec

        # Tuplets
        parts    = lines[i].split(); i += 1
        tupverif = int(parts[0])
        if tupverif:
            d   = int(parts[1])
            s   = [0] + [int(parts[2 + j]) for j in range(d)]   # 1-indexed
            mDD = float(parts[2 + d])
            tuplets_test = TupletsTest(tupletsverif=True, d=d, s=s, mDD=mDD, testtype=_MAX_TYPE)
        else:
            tuplets_test = TupletsTest(tupletsverif=False)

        # Build generator lists — when path is "same", reuse previous list object
        # so that product_no_repeat_ordered enforces C(n,k) automatically.
        dispatch     = Combinaison._get_dispatch()
        gen_lists    = []
        old_gen_list = None
        for path in self.gen_data_files:
            if path.lower() == "same":
                gen_lists.append(old_gen_list)
            else:
                with open(path) as f:
                    gen_type = f.readline().split()[0]
                gen_cls = dispatch.get(gen_type)
                if gen_cls is None:
                    raise NotImplementedError(
                        f"Generator type '{gen_type}' is not yet implemented in Python"
                    )
                old_gen_list = gen_cls.CreateListFromFile(path, Lmax)
                gen_lists.append(old_gen_list)

        # Build Combinaison — prints header, tempering names, component summary
        comb = Combinaison.CreateFromFiles(
            gen_lists, Lmax, temperings, seeds
        )

        # Print remaining search parameters
        if has_tempering:
            print(f"Number of tries per combined generator : {nbtries}")
        if meverif:
            print(f"Upperbound for the sum of dimension gaps for resolutions in psi_12 : {mse}")
        if nb_lines > 0:
            print("delta for particular bits:")
            for jmin, jmax, ec in delta_lines:
                print(f"   >>Bits from {jmin} to {jmax} delta_i={ec}")
        if tupverif:
            print("Verification of DELTA( " +
                  ", ".join(str(s[j]) for j in range(1, d)) +
                  f"{s[d]})")
        print("=" * 68)
        for j, gen_list in enumerate(gen_lists):
            print(f"- Component {j + 1}: {gen_list[0].name()} ")
        sys.stdout.flush()

        me_test = EquidistributionTest(
            L=Lmax, delta=delta, mse=mse, meverif=meverif, method=method
        )
        return comb, me_test, tuplets_test, nbtries, mk_opt, Lmax

    # -- Display helpers --------------------------------------------------

    @staticmethod
    def _display_current_comb(comb: Combinaison) -> None:
        """
        DispCurrentComb: print a summary of the current combined generator.

        Matches the C output format exactly: 66-char separators, generator
        name + parameters, transformation chain for each component.
        """
        print(_EQ66)
        print(f"Number of points   : 2^({comb.k_g:3d})")
        if comb.smax != sys.maxsize:
            print(f"Maximum dimension : {comb.smax:3d}")
        print()

        poly_bv = comb[0].char_poly()
        hw      = bin(poly_bv._val).count('1') + 1   # +1 for the implicit leading x^k term
        print(f"hammingweigth = {hw}")
        sys.stdout.flush()

        for j, comp in enumerate(comb.components):
            gen = comb[j]
            dn = getattr(gen, 'display_name', gen.name)
            print(f"{dn()}:")
            gen.display()
            comp.display()
            if j < comb.J - 1:
                print(_EQ66dash)
            else:
                print(_EQ66)

    @staticmethod
    def _display_summary(
        nbgen: int, nbME: int, nbCF: int, nbsel: int, elapsed: float
    ) -> None:
        """Print the final search summary."""
        print("\n===========================")
        print(f"   Total   =  {nbgen:10d}  ")
        print()
        print(f"     ME    =  {nbME:10d}  ")
        print(f"   CF-ME   =  {nbCF:10d}  ")
        print(f"  retained =  {nbsel:10d}  ")
        print("---------------------------")
        print(f" CPU (sec) =   {elapsed:5.2f}       ")
        print("===========================")
