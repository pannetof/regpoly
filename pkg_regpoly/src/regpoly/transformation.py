"""
transformation.py — Abstract base class for a linear tempering transformation.

A Transformation is callable: t(state) applies the forward transformation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from regpoly.bitvect import BitVect


class Transformation(ABC):
    """
    Abstract base for a linear tempering transformation.

    A transformation is callable: t(state) applies it, t.inverse(state)
    inverts it.

    Attributes
    ----------
    w          : int — current bit width used by this transformation
    w_original : int — original bit width (-1 means inherit generator L at runtime)
    """

    def __init__(self) -> None:
        self.w: int = 0
        self.w_original: int = 0

    @property
    @abstractmethod
    def display_name(self) -> str:
        """Human-readable name used in the search-parameter summary."""

    @abstractmethod
    def display(self) -> None:
        """Print this transformation to stdout."""

    @abstractmethod
    def __call__(self, A: BitVect) -> BitVect:
        """Apply the transformation to A and return the result."""

    @abstractmethod
    def inverse(self, A: BitVect) -> BitVect:
        """Apply the inverse transformation to A and return the result."""

    @abstractmethod
    def update_params(self) -> None:
        """Recompute internal parameters after w changes."""

    @property
    @abstractmethod
    def type_id(self) -> int:
        """Unique integer identifier for this transformation type."""

    @abstractmethod
    def copy(self) -> "Transformation":
        """Return a deep copy of self."""

    # -- File I/O ---------------------------------------------------------

    @classmethod
    def CreateFromFiles(cls, filename: str) -> tuple:
        """
        Read a tempering-transformation file and return all transformations.

        File format:
            <n>
            <type> <param1> ...    (repeated n times, possibly with extra lines)

        Prints the transformation summary as it reads (matching C ReadTempering
        output).

        Returns
        -------
        (transformations, mk_opt)
            transformations : list[Transformation]
            mk_opt          : True if any 'tempMKopt' / 'tempMK2opt' entry found
        """
        from regpoly.transformations.permutation import Permutation
        from regpoly.transformations.temper_mk import TemperMK

        _dispatch = {
            "permut"    : Permutation,
            "tempMK"    : TemperMK,
            "tempMKopt" : TemperMK,
            "tempMK2"   : TemperMK,
            "tempMK2opt": TemperMK,
        }
        _mk_opt_types = {"tempMKopt", "tempMK2opt"}
        _mk2_types    = {"tempMK2",   "tempMK2opt"}

        mk_opt = False
        transformations = []
        with open(filename) as f:
            # Skip blank lines at the start (C's ReadLn skips them)
            line = f.readline()
            while line and line.strip() == '':
                line = f.readline()
            n = int(line)
            for _ in range(n):
                tokens    = f.readline().split()
                type_str  = tokens[0]
                trans_cls = _dispatch.get(type_str, cls)
                full_tokens = tokens + trans_cls._extra_tokens(tokens, f)
                t = trans_cls.read_params(full_tokens)
                transformations.append(t)
                if type_str in _mk_opt_types:
                    mk_opt = True
        return transformations, mk_opt

    @classmethod
    def _extra_tokens(cls, header_tokens: list[str], f) -> list[str]:
        """Read any extra file lines needed beyond the header line. Default: none."""
        return []

    @classmethod
    def read_params(cls, tokens: list[str]) -> "Transformation":
        """
        Parse transformation parameters from the full token list and return
        a new instance.  tokens[0] is the type string; all data needed by
        this transformation (including any multi-line file content) has
        already been flattened into tokens by CreateFromFiles.

        The default implementation raises NotImplementedError so that
        unrecognised transformation types produce a clear error message.
        """
        params = {f"param{i}": v for i, v in enumerate(tokens[1:])}
        raise NotImplementedError(
            f"Transformation type '{tokens[0]}' is not implemented. "
            f"Parsed token dict: {params}"
        )
