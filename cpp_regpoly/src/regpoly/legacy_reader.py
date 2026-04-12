"""
legacy_reader.py — Parsers for the old text-based file formats.

All creation goes through Generateur.create() and Transformation().
Legacy file tags (e.g. "polylcg", "taus") are accepted by
Generateur.create() which resolves them to C++ class names.
"""

from __future__ import annotations

import re
import sys



class LegacyReader:
    """Parsers for legacy text-based generator and transformation files."""

    # -- Generator auto-detect and dispatch --------------------------------

    _GEN_DISPATCH: dict = {}

    @classmethod
    def read_generators(cls, filename: str, L: int) -> list:
        if not cls._GEN_DISPATCH:
            cls._GEN_DISPATCH = {
                "polylcg": cls._read_polylcg,
                "taus"   : cls._read_tausworthe,
                "taus2"  : cls._read_tausworthe,
                "tgfsr"  : cls._read_tgfsr,
                "MT"     : cls._read_mt,
                "genf2w" : cls._read_genf2w,
                "carry"       : cls._read_carry,
                "marsaxorshift": cls._read_marsaxorshift,
                "matsumoto"   : cls._read_matsumoto,
            }
        with open(filename) as f:
            tag = f.readline().split()[0]
        reader = cls._GEN_DISPATCH.get(tag)
        if reader is None:
            raise ValueError(f"Unknown generator type '{tag}' in {filename}")
        return reader(filename, L)

    # -- PolyLCG ----------------------------------------------------------

    @staticmethod
    def _read_polylcg(filename: str, L: int) -> list:
        from regpoly.generateur import Generateur

        with open(filename) as f:
            f.readline()
            n = int(f.readline())
            tokens = []
            for line in f:
                tokens.extend(line.split())

        generators = []
        idx = 0
        for _ in range(n):
            k = int(tokens[idx]); idx += 1
            exponents = []
            while True:
                e = int(tokens[idx]); idx += 1
                exponents.append(e)
                if e == 0:
                    break
            params = {"k": k, "poly": exponents}
            generators.append(Generateur.create("PolyLCG", params, L))
        return generators

    # -- Tausworthe -------------------------------------------------------

    @staticmethod
    def _read_tausworthe(filename: str, L: int) -> list:
        from regpoly.generateur import Generateur

        with open(filename) as f:
            f.readline()
            header = f.readline().split()
            n = int(header[0])
            quicktaus = bool(int(header[1]))
            file_smax = int(header[2]) if not quicktaus else 0

            tokens = []
            for line in f:
                tokens.extend(line.split())

        generators = []
        idx = 0
        for _ in range(n):
            R = []
            while True:
                e = int(tokens[idx]); idx += 1
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
                kmq = k - Q[nq - 2]
                s = kmq

            params = {"poly": Q, "s": s, "quicktaus": quicktaus}
            generators.append(Generateur.create("Tausworthe", params, L))
        return generators

    # -- TGFSR ------------------------------------------------------------

    @staticmethod
    def _read_tgfsr(filename: str, L: int) -> list:
        from regpoly.generateur import Generateur

        with open(filename) as f:
            f.readline()
            parts = f.readline().split()
            w = int(parts[0])
            r = int(parts[1])
            if w > 32:
                raise ValueError(f"Error: w must be <= 32 bits, got {w}")
            nbgen = int(f.readline().split()[0])

            generators = []
            for _ in range(nbgen):
                tokens = f.readline().split()
                a_val = int(tokens[0], 16)
                m_val = int(tokens[1])
                params = {"w": w, "r": r, "m": m_val, "a": a_val}
                generators.append(Generateur.create("TGFSRGen", params, L))
        return generators

    # -- Mersenne Twister -------------------------------------------------

    @staticmethod
    def _read_mt(filename: str, L: int) -> list:
        from regpoly.generateur import Generateur

        with open(filename) as f:
            f.readline()
            nbgen = int(f.readline().split()[0])

            generators = []
            for _ in range(nbgen):
                tokens = f.readline().split()
                w = int(tokens[0])
                r = int(tokens[1])
                m = int(tokens[2])
                p = int(tokens[3])
                a = int(tokens[4], 16)
                if w > 32:
                    raise ValueError(f"w must be <= 32 bits, got {w}")
                params = {"w": w, "r": r, "m": m, "p": p, "a": a}
                generators.append(Generateur.create("MersenneTwister", params, L))
        return generators

    # -- GenF2w -----------------------------------------------------------

    @staticmethod
    def _read_genf2w(filename: str, L: int) -> list:
        from regpoly.generateur import Generateur

        _token_re = re.compile(r'\S+')

        with open(filename) as f:
            f.readline()
            gen_type = int(f.readline().split()[0])
            normal_basis = bool(int(f.readline().split()[0]))
            nbgen = int(f.readline().split()[0])

            def _tokens():
                for line in f:
                    for m in _token_re.finditer(line):
                        yield m.group()

            it = _tokens()
            generators = []
            for _ in range(nbgen):
                w = int(next(it))
                r = int(next(it))
                modM = int(next(it), 16)
                step_count = int(next(it))
                nbcoeff = int(next(it))

                coeffs = []
                for _ in range(nbcoeff):
                    val = int(next(it), 16)
                    pos = int(next(it))
                    coeffs.append({"value": val, "position": pos})

                params = {
                    "w": w, "r": r, "modM": modM,
                    "normal_basis": normal_basis,
                    "step": step_count,
                    "type": "lfsr" if gen_type == 1 else "polylcg",
                    "coeffs": coeffs,
                }
                generators.append(Generateur.create("genf2w", params, L))
        return generators

    # -- Carry ------------------------------------------------------------

    @staticmethod
    def _read_carry(filename: str, L: int) -> list:
        from regpoly.generateur import Generateur

        with open(filename) as f:
            f.readline()  # skip tag
            parts = f.readline().split()
            w = int(parts[0])
            r = int(parts[1])
            p = int(parts[2])
            nbgen = int(f.readline().split()[0])

            generators = []
            for _ in range(nbgen):
                tokens = f.readline().split()
                m1 = int(tokens[0])
                m2 = int(tokens[1])
                m3 = int(tokens[2])

                mat_types = []
                mat_pi = []
                mat_pu = []
                idx = 3
                for _ in range(8):
                    mat_types.append(int(tokens[idx])); idx += 1
                    mat_pi.append(int(tokens[idx])); idx += 1
                    mat_pi.append(int(tokens[idx])); idx += 1
                    mat_pi.append(int(tokens[idx])); idx += 1
                    mat_pu.append(int(tokens[idx], 16)); idx += 1
                    mat_pu.append(int(tokens[idx], 16)); idx += 1
                    mat_pu.append(int(tokens[idx], 16)); idx += 1

                params = {
                    "w": w, "r": r, "p": p,
                    "m1": m1, "m2": m2, "m3": m3,
                    "mat_types": mat_types,
                    "mat_pi": mat_pi,
                    "mat_pu": mat_pu,
                }
                generators.append(Generateur.create("Carry2Gen", params, L))
        return generators

    # -- Matsumoto --------------------------------------------------------

    @staticmethod
    def _read_matsumoto(filename: str, L: int) -> list:
        """
        File format:
            matsumoto
            nbgens
            type n m nbparamsint p0 p1 ... nbparamsunsigned u0 u1 ...
        """
        from regpoly.generateur import Generateur

        with open(filename) as f:
            f.readline()  # skip tag
            nbgens = int(f.readline().split()[0])

            generators = []
            for _ in range(nbgens):
                tokens = f.readline().split()
                idx = 0
                gen_type = int(tokens[idx]); idx += 1
                n = int(tokens[idx]); idx += 1
                m = int(tokens[idx]); idx += 1
                nbpi = int(tokens[idx]); idx += 1
                paramsint = []
                for _ in range(nbpi):
                    paramsint.append(int(tokens[idx])); idx += 1
                nbpu = int(tokens[idx]); idx += 1
                paramsunsigned = []
                for _ in range(nbpu):
                    paramsunsigned.append(int(tokens[idx], 16)); idx += 1

                params = {
                    "type": gen_type, "n": n, "m": m,
                    "paramsint": paramsint,
                    "paramsunsigned": paramsunsigned,
                }
                generators.append(Generateur.create("MatsumotoGen", params, L))
        return generators

    # -- MarsaXorshift ----------------------------------------------------

    @staticmethod
    def _read_marsaxorshift(filename: str, L: int) -> list:
        """
        File format:
            marsaxorshift
            nbgen r w
            numberGenToAllocate
            Per config line (type-dependent):
              TYPE1:  type a b c → creates 4 variants
              TYPE21-25: type m p0 p1 p2 q0 q1 q2 → creates 2-4 variants
              TYPE3:  type N [pos shift]×N
              TYPE4:  type r m a b c d
              TYPE100: type nbmi nbxorshift [mi mi_nb [shifts]]×nbmi
        """
        from regpoly.generateur import Generateur

        TYPE1 = 1
        TYPE21, TYPE22, TYPE23, TYPE24, TYPE25 = 21, 22, 23, 24, 25
        TYPE3 = 3
        TYPE4 = 4
        TYPEGENERAL = 100

        with open(filename) as f:
            f.readline()  # skip tag
            parts = f.readline().split()
            nbgen = int(parts[0])
            r = int(parts[1])
            w = int(parts[2])
            f.readline()  # skip numberGenToAllocate

            generators = []
            for _ in range(nbgen):
                tokens = f.readline().split()
                idx = 0
                gen_type = int(tokens[idx]); idx += 1

                if gen_type == TYPE1:
                    a = int(tokens[idx]); idx += 1
                    b = int(tokens[idx]); idx += 1
                    c = int(tokens[idx]); idx += 1
                    # Create 4 variants: (a,b,c), (c,b,a), (-a,-b,-c), (a,-c,-b)
                    for shifts in [(a,b,c), (c,b,a), (-a,-b,-c), (a,-c,-b)]:
                        params = {"type": 1, "w": w, "r": 1, "shifts": list(shifts)}
                        generators.append(Generateur.create("MarsaXorshiftGen", params, L))

                elif gen_type >= TYPE21 and gen_type <= TYPE25:
                    m_val = int(tokens[idx]); idx += 1
                    op = [int(tokens[idx]) for _ in range(3)]; idx += 3
                    oq = [int(tokens[idx]) for _ in range(3)]; idx += 3

                    # Build variants based on type
                    variants = []
                    if gen_type == TYPE21:
                        variants.append((list(op), list(oq)))
                        variants.append(([-op[0], -op[1], op[2]], [-oq[0], oq[1], oq[2]]))
                        variants.append(([-op[1], -op[0], op[2]], [-oq[0], oq[1], oq[2]]))
                    elif gen_type == TYPE22:
                        variants.append((list(op), list(oq)))
                        variants.append(([-op[0], op[1], op[2]], [-oq[0], -oq[1], oq[2]]))
                        variants.append(([-op[0], op[1], op[2]], [-oq[1], -oq[0], oq[2]]))
                    elif gen_type == TYPE23:
                        variants.append((list(op), list(oq)))
                        variants.append(([-op[0], -op[1], op[2]], [-oq[0], oq[1], oq[2]]))
                    elif gen_type == TYPE24:
                        variants.append((list(op), list(oq)))
                        variants.append(([-op[0], op[1], op[2]], [-oq[0], -oq[1], oq[2]]))
                    elif gen_type == TYPE25:
                        variants.append((list(op), list(oq)))
                        variants.append(([op[2], op[1], op[0]], list(oq)))
                        variants.append(([-op[0], -op[1], -op[2]], list(oq)))
                        variants.append(([op[0], op[2], op[1]], list(oq)))

                    for p_v, q_v in variants:
                        params = {
                            "type": gen_type, "w": w, "r": r, "m": m_val,
                            "p": p_v, "q": q_v,
                        }
                        generators.append(Generateur.create("MarsaXorshiftGen", params, L))

                elif gen_type == TYPE3:
                    n3 = int(tokens[idx]); idx += 1
                    tap_pos = []
                    tap_shifts = []
                    for _ in range(n3):
                        tap_pos.append(int(tokens[idx])); idx += 1
                        tap_shifts.append(int(tokens[idx])); idx += 1
                    params = {
                        "type": 3, "w": w, "r": r,
                        "tap_positions": tap_pos, "tap_shifts": tap_shifts,
                    }
                    generators.append(Generateur.create("MarsaXorshiftGen", params, L))

                elif gen_type == TYPE4:
                    r_val = int(tokens[idx]); idx += 1
                    m_val = int(tokens[idx]); idx += 1
                    a = int(tokens[idx]); idx += 1
                    b = int(tokens[idx]); idx += 1
                    c = int(tokens[idx]); idx += 1
                    d = int(tokens[idx]); idx += 1
                    params = {
                        "type": 4, "w": w, "r": r_val, "m": m_val,
                        "p": [a, b], "q": [c, d],
                    }
                    generators.append(Generateur.create("MarsaXorshiftGen", params, L))

                elif gen_type == TYPEGENERAL:
                    nbmi = int(tokens[idx]); idx += 1
                    nbxorshift = int(tokens[idx]); idx += 1
                    mi_pos = []
                    mi_counts = []
                    mi_shifts = []
                    for _ in range(nbmi):
                        mi_pos.append(int(tokens[idx])); idx += 1
                        cnt = int(tokens[idx]); idx += 1
                        mi_counts.append(cnt)
                        for _ in range(cnt):
                            mi_shifts.append(int(tokens[idx])); idx += 1
                    params = {
                        "type": 100, "w": w, "r": r,
                        "mi_positions": mi_pos,
                        "mi_counts": mi_counts,
                        "mi_shifts": mi_shifts,
                    }
                    generators.append(Generateur.create("MarsaXorshiftGen", params, L))

                else:
                    raise ValueError(f"Unknown marsaxorshift type: {gen_type}")

        return generators

    # -- Transformations --------------------------------------------------

    @staticmethod
    def read_transformations(filename: str) -> tuple:
        """
        Read transformations from the legacy text format.

        Returns (transformations, mk_opt).
        """
        from regpoly.transformation import Transformation

        _mk_opt_types = {"tempMKopt", "tempMK2opt"}
        _mk2_types = {"tempMK2", "tempMK2opt"}

        mk_opt = False
        transformations = []
        with open(filename) as f:
            line = f.readline()
            while line and line.strip() == '':
                line = f.readline()
            n = int(line)
            for _ in range(n):
                tokens = f.readline().split()
                type_str = tokens[0]

                if type_str == "permut":
                    w = int(tokens[1])
                    p = int(tokens[2])
                    q = int(tokens[3])
                    params = {"w": w, "p": p, "q": q}
                    random_specs = {}
                    if p == -1:
                        random_specs["p"] = {"random": "coprime", "mod": "w"}
                        params["p"] = 0
                    if q == -1:
                        random_specs["q"] = {"random": "int", "min": 0, "max": "w"}
                        params["q"] = 0
                    trans_type = "permut"

                elif type_str in ("tempMK", "tempMKopt", "tempMK2", "tempMK2opt"):
                    w = int(tokens[1])
                    eta = int(tokens[2])
                    mu = int(tokens[3])
                    idx = 4

                    if type_str in _mk2_types:
                        trans_type = "tempMK2"
                        u = int(tokens[idx]); idx += 1
                        l = int(tokens[idx]); idx += 1
                    else:
                        trans_type = "tempMK"
                        u = l = 0

                    nb_words = int(tokens[idx]); idx += 1

                    if "opt" in type_str:
                        idx += 1  # skip disp_progress
                        idx += 1  # skip limit_v

                    random_specs = {}
                    if nb_words != -1:
                        b_tokens = f.readline().split()
                        c_tokens = f.readline().split()
                        b_val = LegacyReader._read_hex_words(b_tokens, nb_words, w)
                        c_val = LegacyReader._read_hex_words(c_tokens, nb_words, w)
                    else:
                        b_val = 0
                        c_val = 0
                        random_specs = {
                            "b": {"random": "bitvect", "bits": "w"},
                            "c": {"random": "bitvect", "bits": "w"},
                        }

                    params = {
                        "w": w, "eta": eta, "mu": mu,
                        "u": u, "l": l,
                        "b": b_val, "c": c_val,
                    }
                else:
                    raise ValueError(f"Unknown transformation type '{type_str}'")

                t = Transformation(trans_type, params, w_original=w)
                t._random_specs = random_specs
                transformations.append(t)
                if type_str in _mk_opt_types:
                    mk_opt = True

        return transformations, mk_opt

    @staticmethod
    def _read_hex_words(hex_tokens: list, nb_words: int, w: int) -> int:
        val = 0
        for token in hex_tokens[:nb_words]:
            val = (val << 32) | int(token, 16)
        shift = 32 * nb_words - w
        val >>= shift
        return val
