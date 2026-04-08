"""
matrix.py — BitMatrix class for binary matrices over GF(2).

Internal storage: a list of Python ints, each representing a packed row.
Bit 0 (MSB) is at the highest position of the int.
"""

from regpoly.bitvect import BitVect


_BLOCK = [" ", "▄", "▀", "█"]  # index = top*2 + bot


class BitMatrix:
    __slots__ = ("nblignes", "l", "_rows")

    def __init__(self, nblignes: int, l: int):
        """
        Allocate a zero-filled binary matrix.

        nblignes : number of rows
        l        : number of columns (bits per row)
        """
        self.nblignes: int = nblignes
        self.l: int = l
        self._rows: list[int] = [0] * nblignes

    def copy(self) -> "BitMatrix":
        m = BitMatrix(self.nblignes, self.l)
        m._rows = self._rows[:]
        return m

    def get_bit(self, row: int, col: int) -> int:
        """Get bit at (row, col). col 0 = MSB."""
        return (self._rows[row] >> (self.l - 1 - col)) & 1

    def set_bit(self, row: int, col: int, val: int) -> None:
        """Set bit at (row, col). col 0 = MSB."""
        pos = self.l - 1 - col
        if val:
            self._rows[row] |= (1 << pos)
        else:
            self._rows[row] &= ~(1 << pos)

    def get_bitvect(self, row: int) -> BitVect:
        """Return row as a BitVect."""
        return BitVect(self.l, self._rows[row])

    def set_bitvect(self, row: int, bv: BitVect) -> None:
        """Set row from a BitVect."""
        self._rows[row] = bv._val & ((1 << self.l) - 1)

    # ── Display ──────────────────────────────────────────────────────────

    def display(self, fmt: str = "raw") -> str:
        """
        Render the matrix as a string.

        fmt:
            "raw"     — binary 0/1 with row numbers
            "hex4"    — hexadecimal, 4-bit groups
            "hex8"    — hexadecimal, 8-bit groups
            "block"   — Unicode block characters (1 col x 2 rows)
            "png"     — save as PNG bitmap, returns the file path
        """
        if fmt == "raw":
            return self._fmt_raw()
        elif fmt == "hex4":
            return self._fmt_hex(4)
        elif fmt == "hex8":
            return self._fmt_hex(8)
        elif fmt == "block":
            return self._fmt_block()
        else:
            raise ValueError(f"Unknown format: {fmt!r}")

    def to_image(self, scale: int = 8):
        """
        Render the matrix as a PIL Image (black-and-white).

        scale : pixels per bit (1 for huge matrices, 8-16 for small ones).
        Returns a PIL.Image.Image object.
        """
        try:
            from PIL import Image
        except ImportError:
            raise ImportError("Pillow is required for image output: pip install Pillow")

        import numpy as np
        M = self._to_numpy()
        img_data = ((1 - M) * 255).astype(np.uint8)
        img = Image.fromarray(img_data, mode="L")
        if scale > 1:
            w, h = img.size
            img = img.resize((w * scale, h * scale), Image.NEAREST)
        return img

    # ── Internal format methods ──────────────────────────────────────────

    def _bit(self, row: int, col: int) -> int:
        """Fast bit access without bounds check."""
        if col >= self.l:
            return 0
        return (self._rows[row] >> (self.l - 1 - col)) & 1

    def _fmt_raw(self) -> str:
        lines = []
        for r in range(self.nblignes):
            bits = "".join(
                str((self._rows[r] >> (self.l - 1 - c)) & 1)
                for c in range(self.l)
            )
            lines.append(f"{r:3d}| {bits}")
        return "\n".join(lines)

    def _fmt_hex(self, bits_per_group: int) -> str:
        digits = (bits_per_group + 3) // 4
        lines = []
        for r in range(self.nblignes):
            groups = []
            for c in range(0, self.l, bits_per_group):
                val = 0
                for b in range(bits_per_group):
                    if c + b < self.l:
                        val = (val << 1) | self._bit(r, c + b)
                    else:
                        val <<= 1
                groups.append(f"{val:0{digits}X}")
            lines.append(f"{r:3d}| {' '.join(groups)}")
        return "\n".join(lines)


    def _fmt_block(self) -> str:
        nrows = self.nblignes
        ncols = self.l
        lines = []
        for r in range(0, nrows, 2):
            chars = []
            for c in range(ncols):
                top = self._bit(r, c)
                bot = self._bit(r + 1, c) if r + 1 < nrows else 0
                chars.append(_BLOCK[top * 2 + bot])
            lines.append("".join(chars))
        return "\n".join(lines)

    def _to_numpy(self):
        import numpy as np
        M = np.zeros((self.nblignes, self.l), dtype=np.uint8)
        for r in range(self.nblignes):
            for c in range(self.l):
                M[r, c] = self._bit(r, c)
        return M

    def __repr__(self) -> str:
        return f"BitMatrix(nblignes={self.nblignes}, l={self.l})"
