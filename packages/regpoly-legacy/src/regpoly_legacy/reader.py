# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Pure-Python parser for the pre-v2 ``.dat`` parameter file format.

Mirrors the now-deleted C++ ``legacy_reader.cpp`` byte-for-byte: same tag
dispatch, same per-tag field layouts, same error wording. Produces
``(family_canonical_name, params_dict)`` tuples that feed straight into
``regpoly.Generator.create(family, L, **params)``.

Tag → family mapping (from C++ ``read_generator_specs`` dispatch):

    polylcg                 → PolyLCG
    taus / taus2            → Tausworthe
    tgfsr                   → TGFSR
    MT                      → MersenneTwister
    genf2w (gen_type==1)    → GenF2wLFSR
    genf2w (gen_type!=1)    → GenF2wPolyLCG
    carry                   → WELLGen        (legacy "Carry2Gen" file alias)
    marsaxorshift           → MarsaXorshiftGen

All numeric parsing uses base-10 (``parse_int``) or base-16 (``parse_hex``)
with explicit prefix stripping — the C++ ``std::stoull`` accepts ``0x``
in hex mode, so we replicate that with a leading-prefix strip before
``int(s, 16)``. Unbounded Python ints are masked to 32 / 64 bits at the
same points the C++ relied on ``uint32_t`` / ``uint64_t`` overflow.
"""

from __future__ import annotations

from typing import Any

from regpoly_legacy.well_legacy_decode import decode_d_s_mask, decode_test_mask

_M32 = 0xFFFFFFFF
_M64 = 0xFFFFFFFFFFFFFFFF


# ── Numeric parsing ─────────────────────────────────────────────────────


def _parse_int(s: str, filename: str) -> int:
    """Mirror C++ ``parse_int``: base-10, signed, rejects trailing junk.

    Python ``int(s, 10)`` accepts optional leading sign and rejects
    trailing junk by default. Matches C++ ``std::stoll(..., 10)`` with
    a ``pos != s.size()`` check.
    """
    try:
        return int(s, 10)
    except ValueError as exc:
        raise RuntimeError(
            f"legacy_reader: not an integer: '{s}' (in {filename})"
        ) from exc


def _parse_hex(s: str, filename: str) -> int:
    """Mirror C++ ``parse_hex``: base-16, unsigned, ``0x``/``0X`` prefix
    is accepted (the C++ ``std::stoull`` accepts it; Python ``int(s, 16)``
    does NOT unless we strip first).
    """
    raw = s
    if s.startswith(("0x", "0X")):
        s = s[2:]
    try:
        v = int(s, 16)
    except ValueError as exc:
        raise RuntimeError(
            f"legacy_reader: not a hex literal: '{raw}' (in {filename})"
        ) from exc
    if v < 0:
        raise RuntimeError(
            f"legacy_reader: not a hex literal: '{raw}' (in {filename})"
        )
    return v


# ── Tokenisation helpers ────────────────────────────────────────────────


def _tokenise_rest(lines: list[str], start: int) -> list[str]:
    """Whitespace-split the remainder of ``lines`` starting at ``start``
    into a flat token list (mirrors C++ ``tokenise_rest``)."""
    out: list[str] = []
    for line in lines[start:]:
        out.extend(line.split())
    return out


def _next_nonempty_line(lines: list[str], idx: int, filename: str) -> tuple[list[str], int]:
    """Return tokens of the next non-empty line and updated index.
    Raises if EOF.
    """
    while idx < len(lines):
        toks = lines[idx].split()
        idx += 1
        if toks:
            return toks, idx
    raise RuntimeError(f"legacy_reader: unexpected EOF in {filename}")


# ── File loading ────────────────────────────────────────────────────────


def _read_lines(filename: str) -> list[str]:
    # latin-1 mirrors the C++ side, which reads raw bytes via std::getline
    # without any encoding awareness. Some fixtures embed non-UTF-8 bytes
    # in their human-readable comments (example2 has French text from the
    # original C codebase).
    try:
        with open(filename, "r", encoding="latin-1") as fh:
            return fh.readlines()
    except OSError as exc:
        raise RuntimeError(
            f"legacy_reader: cannot open file: {filename}"
        ) from exc


# ── Cursor over flat token list (Python "tokens" loops) ─────────────────


class _Cursor:
    """Bounded token-stream cursor; mirrors C++ ``Cursor`` struct."""

    __slots__ = ("_toks", "_i", "_filename")

    def __init__(self, toks: list[str], filename: str):
        self._toks = toks
        self._i = 0
        self._filename = filename

    def take(self) -> str:
        if self._i >= len(self._toks):
            raise RuntimeError(
                f"legacy_reader: unexpected end of file in {self._filename}"
            )
        v = self._toks[self._i]
        self._i += 1
        return v

    def take_int(self) -> int:
        return _parse_int(self.take(), self._filename)

    def take_hex(self) -> int:
        return _parse_hex(self.take(), self._filename)


# ── tempMK helper: assemble b/c from hex words ──────────────────────────


def _read_hex_words(tokens: list[str], nb_words: int, w: int, filename: str) -> int:
    """Mirror C++ ``read_hex_words``: pack ``nb_words`` × 32-bit hex values
    into one integer, then right-shift to fit ``w`` bits.

    The C++ enforces ``32 * nb_words <= 64`` (so the result fits in a
    ``uint64_t``). We replicate the same guard so behavior matches: in
    Python the int could be wider, but raising at the same boundary
    preserves identical error semantics.
    """
    assert nb_words >= 1
    if 32 * nb_words > 64:
        raise RuntimeError(
            f"legacy_reader: hex words too wide (32*nb_words="
            f"{32 * nb_words} > 64) in {filename}"
        )
    if len(tokens) < nb_words:
        raise RuntimeError(
            f"legacy_reader: expected {nb_words} hex word(s), got "
            f"{len(tokens)} in {filename}"
        )
    val = 0
    for i in range(nb_words):
        val = ((val << 32) | _parse_hex(tokens[i], filename)) & _M64
    shift = 32 * nb_words - w
    if shift > 0:
        val >>= shift
    return val


# ── Per-tag generator readers ───────────────────────────────────────────


def _read_polylcg(lines: list[str], filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    # Skip tag line (line 0).
    if len(lines) < 2:
        raise RuntimeError(f"legacy_reader: missing nbgen line in {filename}")
    ntoks = lines[1].split()
    if not ntoks:
        raise RuntimeError(f"legacy_reader: empty nbgen line in {filename}")
    n = _parse_int(ntoks[0], filename)

    toks = _tokenise_rest(lines, 2)
    cur = _Cursor(toks, filename)

    out: list[tuple[str, dict[str, Any]]] = []
    for _ in range(n):
        k = cur.take_int()
        exponents: list[int] = []
        while True:
            e = cur.take_int()
            exponents.append(e)
            if e == 0:
                break
        out.append(("PolyLCG", {"k": k, "poly": exponents}))
    return out


def _read_tausworthe(lines: list[str], filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    # Skip tag line; read header on next non-empty line.
    idx = 1
    header, idx = _next_nonempty_line(lines, idx, filename)
    if len(header) < 2:
        raise RuntimeError(
            f"legacy_reader: tausworthe header needs >= 2 tokens in {filename}"
        )
    n = _parse_int(header[0], filename)
    quicktaus = _parse_int(header[1], filename) != 0
    file_smax = 0
    if not quicktaus:
        if len(header) < 3:
            raise RuntimeError(
                f"legacy_reader: tausworthe header needs smax token "
                f"when not quicktaus, in {filename}"
            )
        file_smax = _parse_int(header[2], filename)

    toks = _tokenise_rest(lines, idx)
    cur = _Cursor(toks, filename)

    out: list[tuple[str, dict[str, Any]]] = []
    for _ in range(n):
        R: list[int] = []
        while True:
            e = cur.take_int()
            R.append(e)
            if e == 0:
                break
        Q = list(reversed(R))
        nq = len(Q)
        k = Q[nq - 1]
        if not quicktaus:
            s = min(file_smax, L - k)
            if s < 1:
                s = 1
        else:
            s = k - Q[nq - 2]
        out.append((
            "Tausworthe",
            {"poly": Q, "s": s, "quicktaus": quicktaus},
        ))
    return out


def _read_tgfsr(lines: list[str], filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    idx = 1
    header, idx = _next_nonempty_line(lines, idx, filename)
    if len(header) < 2:
        raise RuntimeError(
            f"legacy_reader: tgfsr header needs (w, r) in {filename}"
        )
    w = _parse_int(header[0], filename)
    r = _parse_int(header[1], filename)
    if w > 32:
        raise RuntimeError(
            f"legacy_reader: tgfsr w must be <= 32, got {w} in {filename}"
        )
    count_line, idx = _next_nonempty_line(lines, idx, filename)
    nbgen = _parse_int(count_line[0], filename)

    out: list[tuple[str, dict[str, Any]]] = []
    for _ in range(nbgen):
        toks, idx = _next_nonempty_line(lines, idx, filename)
        if len(toks) < 2:
            raise RuntimeError(
                f"legacy_reader: tgfsr row needs (a_hex, m) in {filename}"
            )
        a_val = _parse_hex(toks[0], filename)
        m_val = _parse_int(toks[1], filename)
        out.append((
            "TGFSR",
            {"w": w, "r": r, "m": m_val, "a": a_val},
        ))
    return out


def _read_mt(lines: list[str], filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    idx = 1
    count_line, idx = _next_nonempty_line(lines, idx, filename)
    nbgen = _parse_int(count_line[0], filename)

    out: list[tuple[str, dict[str, Any]]] = []
    for _ in range(nbgen):
        toks, idx = _next_nonempty_line(lines, idx, filename)
        if len(toks) < 5:
            raise RuntimeError(
                f"legacy_reader: MT row needs (w, r, m, p, a_hex) in {filename}"
            )
        w = _parse_int(toks[0], filename)
        r = _parse_int(toks[1], filename)
        m = _parse_int(toks[2], filename)
        p_v = _parse_int(toks[3], filename)
        a = _parse_hex(toks[4], filename)
        if w > 32:
            raise RuntimeError(
                f"legacy_reader: MT w must be <= 32, got {w} in {filename}"
            )
        out.append((
            "MersenneTwister",
            {"w": w, "r": r, "m": m, "p": p_v, "a": a},
        ))
    return out


def _read_genf2w(lines: list[str], filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    idx = 1
    gen_type_line, idx = _next_nonempty_line(lines, idx, filename)
    gen_type = _parse_int(gen_type_line[0], filename)
    nb_line, idx = _next_nonempty_line(lines, idx, filename)
    nb_token = _parse_int(nb_line[0], filename)
    normal_basis = nb_token != 0
    count_line, idx = _next_nonempty_line(lines, idx, filename)
    nbgen = _parse_int(count_line[0], filename)

    toks = _tokenise_rest(lines, idx)
    cur = _Cursor(toks, filename)

    fam = "GenF2wLFSR" if gen_type == 1 else "GenF2wPolyLCG"

    out: list[tuple[str, dict[str, Any]]] = []
    for _ in range(nbgen):
        w = cur.take_int()
        r = cur.take_int()
        modM = cur.take_hex()
        step = cur.take_int()
        nbcoeff = cur.take_int()

        coeff_vals: list[int] = []
        coeff_pos: list[int] = []
        for _ in range(nbcoeff):
            v = cur.take_hex()
            pos = cur.take_int()
            coeff_vals.append(v)
            coeff_pos.append(pos)

        out.append((
            fam,
            {
                "w": w,
                "r": r,
                "modM": modM,
                "normal_basis": normal_basis,
                "step": step,
                "type": "lfsr" if gen_type == 1 else "polylcg",
                "coeff": coeff_vals,
                "nocoeff": coeff_pos,
            },
        ))
    return out


# Legacy 0..7 type indices remap to paper Mi (0..6). Old type 6 (multi-shift
# XOR) has no paper equivalent — rejected at parse time.
_OLD_TO_NEW_MI = (3, 1, 4, 2, 5, 6, -1, 0)


def _read_carry(lines: list[str], filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    idx = 1
    header, idx = _next_nonempty_line(lines, idx, filename)
    if len(header) < 3:
        raise RuntimeError(
            f"legacy_reader: carry header needs (w, r, p) in {filename}"
        )
    w = _parse_int(header[0], filename)
    r = _parse_int(header[1], filename)
    p = _parse_int(header[2], filename)
    count_line, idx = _next_nonempty_line(lines, idx, filename)
    nbgen = _parse_int(count_line[0], filename)

    out: list[tuple[str, dict[str, Any]]] = []
    for _ in range(nbgen):
        toks, idx = _next_nonempty_line(lines, idx, filename)
        # Per Python: 3 m's + 8 × (1 type + 3 pi + 3 pu) = 3 + 8*7 = 59.
        if len(toks) < 59:
            raise RuntimeError(
                f"legacy_reader: carry row needs >=59 tokens, got "
                f"{len(toks)} in {filename}"
            )
        m1 = _parse_int(toks[0], filename)
        m2 = _parse_int(toks[1], filename)
        m3 = _parse_int(toks[2], filename)

        matrices: dict[str, dict[str, Any]] = {}
        cursor = 3
        for k in range(8):
            raw_type = _parse_int(toks[cursor], filename); cursor += 1
            pi0 = _parse_int(toks[cursor], filename); cursor += 1
            _pi1 = _parse_int(toks[cursor], filename); cursor += 1  # noqa: F841
            _pi2 = _parse_int(toks[cursor], filename); cursor += 1  # noqa: F841
            pu0 = _parse_hex(toks[cursor], filename); cursor += 1
            pu1 = _parse_hex(toks[cursor], filename); cursor += 1
            pu2 = _parse_hex(toks[cursor], filename); cursor += 1

            if raw_type < 0 or raw_type >= 8 or _OLD_TO_NEW_MI[raw_type] < 0:
                raise RuntimeError(
                    f"legacy_reader: carry .dat in {filename} uses "
                    f"obsolete WELL transformation type {raw_type} "
                    f"(no paper Mi equivalent; was the multi-shift extension)"
                )
            Mi = _OLD_TO_NEW_MI[raw_type]
            slot = f"T{k}"
            ctx = f"{filename} (slot {slot})"
            entry: dict[str, Any] = {"M": Mi}
            if Mi in (0, 1):
                pass  # no args
            elif Mi in (2, 3):
                entry["t"] = pi0
            elif Mi == 4:
                entry["a"] = pu0 & _M32
            elif Mi == 5:
                entry["t"] = pi0
                entry["b"] = pu0 & _M32
            elif Mi == 6:
                s = decode_d_s_mask(pu1, ctx)
                t_bit = decode_test_mask(pu2, ctx)
                entry["q"] = pi0
                entry["t"] = t_bit
                entry["s"] = s
                entry["a"] = pu0 & _M32
            matrices[slot] = entry

        out.append((
            "WELLGen",
            {
                "w": w, "r": r, "p": p,
                "m1": m1, "m2": m2, "m3": m3,
                "matrices": matrices,
            },
        ))
    return out


def _read_marsaxorshift(lines: list[str], filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    idx = 1
    header, idx = _next_nonempty_line(lines, idx, filename)
    if len(header) < 3:
        raise RuntimeError(
            f"legacy_reader: marsaxorshift header needs (nbgen, r, w) in {filename}"
        )
    nbgen = _parse_int(header[0], filename)
    r = _parse_int(header[1], filename)
    w = _parse_int(header[2], filename)
    _skip, idx = _next_nonempty_line(lines, idx, filename)  # numberGenToAllocate

    def push_t1(dst: list[tuple[str, dict[str, Any]]], a: int, b: int, c: int) -> None:
        dst.append((
            "MarsaXorshiftGen",
            {"type": 1, "w": w, "r": 1, "shifts": [a, b, c]},
        ))

    def push_t2x(
        dst: list[tuple[str, dict[str, Any]]],
        m_val: int,
        pv: list[int],
        qv: list[int],
    ) -> None:
        dst.append((
            "MarsaXorshiftGen",
            {"type": 2, "w": w, "r": r, "m": m_val, "p": pv, "q": qv},
        ))

    out: list[tuple[str, dict[str, Any]]] = []
    for _ in range(nbgen):
        toks, idx = _next_nonempty_line(lines, idx, filename)
        i = 0

        def take() -> int:
            nonlocal i
            if i >= len(toks):
                raise RuntimeError(
                    f"legacy_reader: marsaxorshift row truncated in {filename}"
                )
            v = _parse_int(toks[i], filename)
            i += 1
            return v

        gen_type = take()

        if gen_type == 1:
            a = take(); b = take(); c = take()
            push_t1(out, a, b, c)
            push_t1(out, c, b, a)
            push_t1(out, -a, -b, -c)
            push_t1(out, a, -c, -b)

        elif 21 <= gen_type <= 25:
            m_val = take()
            op = [take(), take(), take()]
            oq = [take(), take(), take()]
            variants: list[tuple[list[int], list[int]]] = []
            if gen_type == 21:
                variants.append((op, oq))
                variants.append(([-op[0], -op[1],  op[2]], [-oq[0],  oq[1],  oq[2]]))
                variants.append(([-op[1], -op[0],  op[2]], [-oq[0],  oq[1],  oq[2]]))
            elif gen_type == 22:
                variants.append((op, oq))
                variants.append(([-op[0],  op[1],  op[2]], [-oq[0], -oq[1],  oq[2]]))
                variants.append(([-op[0],  op[1],  op[2]], [-oq[1], -oq[0],  oq[2]]))
            elif gen_type == 23:
                variants.append((op, oq))
                variants.append(([-op[0], -op[1],  op[2]], [-oq[0],  oq[1],  oq[2]]))
            elif gen_type == 24:
                variants.append((op, oq))
                variants.append(([-op[0],  op[1],  op[2]], [-oq[0], -oq[1],  oq[2]]))
            else:  # 25
                variants.append((op, oq))
                variants.append(([ op[2],  op[1],  op[0]], oq))
                variants.append(([-op[0], -op[1], -op[2]], oq))
                variants.append(([ op[0],  op[2],  op[1]], oq))
            for pv, qv in variants:
                push_t2x(out, m_val, pv, qv)

        elif gen_type == 3:
            n3 = take()
            tap_pos: list[int] = []
            tap_shifts: list[int] = []
            for _ in range(n3):
                tap_pos.append(take())
                tap_shifts.append(take())
            out.append((
                "MarsaXorshiftGen",
                {"type": 3, "w": w, "r": r,
                 "tap_positions": tap_pos, "tap_shifts": tap_shifts},
            ))

        elif gen_type == 4:
            r_val = take(); m_val = take()
            a = take(); b = take(); c = take(); d = take()
            out.append((
                "MarsaXorshiftGen",
                {"type": 4, "w": w, "r": r_val, "m": m_val,
                 "p": [a, b], "q": [c, d]},
            ))

        elif gen_type == 100:
            nbmi = take()
            _nbxorshift = take()  # accepted but unused — Python ignores it
            mi_pos: list[int] = []
            mi_counts: list[int] = []
            mi_shifts: list[int] = []
            for _ in range(nbmi):
                mi_pos.append(take())
                cnt = take()
                mi_counts.append(cnt)
                for _ in range(cnt):
                    mi_shifts.append(take())
            out.append((
                "MarsaXorshiftGen",
                {"type": 100, "w": w, "r": r,
                 "mi_positions": mi_pos, "mi_counts": mi_counts,
                 "mi_shifts": mi_shifts},
            ))

        else:
            raise RuntimeError(
                f"legacy_reader: unknown marsaxorshift type {gen_type} in {filename}"
            )
    return out


# ── Public dispatch entry points ────────────────────────────────────────


_GENERATOR_DISPATCH = {
    "polylcg":       _read_polylcg,
    "taus":          _read_tausworthe,
    "taus2":         _read_tausworthe,
    "tgfsr":         _read_tgfsr,
    "MT":            _read_mt,
    "genf2w":        _read_genf2w,
    "carry":         _read_carry,
    "marsaxorshift": _read_marsaxorshift,
}


def parse_generator_specs(filename: str, L: int) -> list[tuple[str, dict[str, Any]]]:
    """Parse a legacy ``.dat`` file into ``[(family, params_dict), ...]``.

    The first non-empty line's first token selects the per-tag parser.
    Each ``params_dict`` matches the schema accepted by
    ``regpoly.Generator.create(family, L, **params_dict)``.
    """
    lines = _read_lines(filename)
    if not lines:
        raise RuntimeError(f"legacy_reader: empty file: {filename}")
    first = lines[0].split()
    if not first:
        raise RuntimeError(f"legacy_reader: missing tag on line 1 of {filename}")
    tag = first[0]
    reader = _GENERATOR_DISPATCH.get(tag)
    if reader is None:
        raise RuntimeError(
            f"legacy_reader: unknown generator tag '{tag}' in {filename}"
        )
    return reader(lines, filename, L)


def parse_transformation_specs(filename: str) -> tuple[list[tuple[str, dict[str, Any]]], bool]:
    """Parse a legacy transformation ``.dat`` file.

    Returns ``(specs, mk_opt)`` where ``specs`` is ``[(trans_type, params), ...]``
    and ``mk_opt`` is True iff at least one row used a ``*opt`` variant.
    """
    lines = _read_lines(filename)
    if not lines:
        raise RuntimeError(f"legacy_reader: empty file: {filename}")

    idx = 0
    count_line, idx = _next_nonempty_line(lines, idx, filename)
    n = _parse_int(count_line[0], filename)

    mk_opt = False
    specs: list[tuple[str, dict[str, Any]]] = []
    for _ in range(n):
        toks, idx = _next_nonempty_line(lines, idx, filename)
        if not toks:
            raise RuntimeError(
                f"legacy_reader: empty transformation row in {filename}"
            )
        type_str = toks[0]
        params: dict[str, Any] = {}

        if type_str == "permut":
            if len(toks) < 4:
                raise RuntimeError(
                    f"legacy_reader: permut row needs (w, p, q) in {filename}"
                )
            w = _parse_int(toks[1], filename)
            p_v = _parse_int(toks[2], filename)
            q_v = _parse_int(toks[3], filename)
            # -1 sentinels → 0 (the C++ PermutationTrans factory will randomise
            # from 0 the same way Python does via its random_specs path).
            if p_v == -1:
                p_v = 0
            if q_v == -1:
                q_v = 0
            params = {"w": w, "p": p_v, "q": q_v}
            trans_type = "permut"

        elif type_str in ("tempMK", "tempMKopt", "tempMK2", "tempMK2opt"):
            if len(toks) < 4:
                raise RuntimeError(
                    f"legacy_reader: tempMK[2][opt] row needs at least "
                    f"(w, eta, mu) in {filename}"
                )
            w = _parse_int(toks[1], filename)
            eta = _parse_int(toks[2], filename)
            mu = _parse_int(toks[3], filename)
            cursor = 4
            is_mk2 = type_str in ("tempMK2", "tempMK2opt")
            is_opt = type_str in ("tempMKopt", "tempMK2opt")
            u = 0
            l = 0
            if is_mk2:
                if cursor + 1 >= len(toks):
                    raise RuntimeError(
                        f"legacy_reader: tempMK2 row missing (u, l) in {filename}"
                    )
                u = _parse_int(toks[cursor], filename); cursor += 1
                l = _parse_int(toks[cursor], filename); cursor += 1
                trans_type = "tempMK2"
            else:
                trans_type = "tempMK"

            if cursor >= len(toks):
                raise RuntimeError(
                    f"legacy_reader: tempMK row missing nb_words in {filename}"
                )
            nb_words = _parse_int(toks[cursor], filename); cursor += 1

            if is_opt:
                # disp_progress + limit_v (consumed but unused)
                if cursor + 1 >= len(toks):
                    raise RuntimeError(
                        f"legacy_reader: tempMK*opt row missing "
                        f"(disp_progress, limit_v) in {filename}"
                    )
                cursor += 2

            params = {"w": w, "eta": eta, "mu": mu, "u": u, "l": l}
            if nb_words != -1:
                b_tokens, idx = _next_nonempty_line(lines, idx, filename)
                c_tokens, idx = _next_nonempty_line(lines, idx, filename)
                b = _read_hex_words(b_tokens, nb_words, w, filename)
                c = _read_hex_words(c_tokens, nb_words, w, filename)
                params["b"] = b
                params["c"] = c
            # else: b/c omitted → randomized later by fill_params.

            if type_str in ("tempMKopt", "tempMK2opt"):
                mk_opt = True

        else:
            raise RuntimeError(
                f"legacy_reader: unknown transformation type '{type_str}' "
                f"in {filename}"
            )

        specs.append((trans_type, params))

    return specs, mk_opt


# ── Public class API ────────────────────────────────────────────────────


class LegacyReader:
    """Pre-v2 `.dat` parser façade.

    Drop-in replacement for the historical
    `regpoly.io.legacy_reader.LegacyReader` (deleted when the legacy
    code moved out of `regpoly` into this add-on). Constructs
    :class:`regpoly.core.generator.Generator` /
    :class:`regpoly.core.transformation.Transformation`
    objects via the existing C++ factory — no behavioural drift versus
    the deleted C++ parser.

    The lower-level
    :func:`regpoly_legacy.reader.parse_generator_specs`
    and
    :func:`regpoly_legacy.reader.parse_transformation_specs`
    functions return raw `(family, params_dict)` tuples for callers
    that want to inspect parsed specs before construction.
    """

    @staticmethod
    def read_generators(filename: str, L: int) -> list:
        """Parse a `.dat` generator file and build `Generator` instances.

        Parameters
        ----------
        filename
            Path to a legacy `.dat` generator file.
        L
            Output resolution in bits, applied to every generator
            in the file.

        Returns
        -------
        list of Generator
            List of :class:`regpoly.core.generator.Generator`
            instances, one per row in the input file.

        Raises
        ------
        RuntimeError
            For any parse failure — see the per-tag
            error messages in
            :func:`regpoly_legacy.reader.parse_generator_specs`.
        """
        from regpoly.core.generator import Generator
        specs = parse_generator_specs(filename, L)
        return [Generator.create(family, L, **params) for family, params in specs]

    @staticmethod
    def read_transformations(filename: str) -> tuple:
        """Parse a `.dat` transformation file and build `Transformation` instances.

        Parameters
        ----------
        filename
            Path to a legacy `.dat` transformation file.

        Returns
        -------
        tuple
            A tuple ``(transformations, mk_opt)``:

            - `transformations`: list of
              :class:`regpoly.core.transformation.Transformation`
              instances.
            - `mk_opt` (`bool`): True iff at least one row used a
              `tempMKopt` / `tempMK2opt` variant.
        """
        from regpoly.core.transformation import Transformation
        specs, mk_opt = parse_transformation_specs(filename)
        return (
            [Transformation.create(t_type, **params) for t_type, params in specs],
            mk_opt,
        )
