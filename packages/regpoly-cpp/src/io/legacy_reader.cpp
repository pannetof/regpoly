// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 3.3: legacy ".dat" parameter file reader (C++ port of
// packages/regpoly/src/regpoly/io/legacy_reader.py). See the header
// for the public API contract.

#include "legacy_reader.h"

#include "factory.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace regpoly_legacy {

namespace {

// ── Tokenisation helpers ──────────────────────────────────────────────

// Whitespace-split the entire remainder of a stream into tokens.
// Mirrors the Python idiom `for line in f: tokens.extend(line.split())`.
std::vector<std::string> tokenise_rest(std::istream& in) {
    std::vector<std::string> out;
    std::string tok;
    while (in >> tok) out.push_back(std::move(tok));
    return out;
}

// Tokenise the next non-empty line. Throws on EOF.
std::vector<std::string>
tokenise_line(std::istream& in, const std::string& filename) {
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> toks;
        std::istringstream ls(line);
        std::string tok;
        while (ls >> tok) toks.push_back(std::move(tok));
        if (!toks.empty()) return toks;
    }
    throw std::runtime_error("legacy_reader: unexpected EOF in " + filename);
}

// Tokenise the next line, returning empty vector on EOF (no throw).
std::vector<std::string> tokenise_line_or_empty(std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> toks;
        std::istringstream ls(line);
        std::string tok;
        while (ls >> tok) toks.push_back(std::move(tok));
        if (!toks.empty()) return toks;
    }
    return {};
}

// ── Numeric parsing ───────────────────────────────────────────────────

int64_t parse_dec(const std::string& s, const std::string& filename) {
    try {
        size_t pos = 0;
        int64_t v = std::stoll(s, &pos, 10);
        if (pos != s.size()) {
            throw std::runtime_error(
                "legacy_reader: trailing junk in integer '"
                + s + "' (in " + filename + ")");
        }
        return v;
    } catch (const std::invalid_argument&) {
        throw std::runtime_error(
            "legacy_reader: not an integer: '" + s
            + "' (in " + filename + ")");
    } catch (const std::out_of_range&) {
        throw std::runtime_error(
            "legacy_reader: integer out of range: '" + s
            + "' (in " + filename + ")");
    }
}

int parse_int(const std::string& s, const std::string& filename) {
    int64_t v = parse_dec(s, filename);
    return static_cast<int>(v);
}

// Hex parse via stoull — handles values up to UINT64_MAX. The Python
// `int(x, 16)` returns an unbounded int; the legacy file format only
// uses hex for values that fit in 64 bits (typically 32 bits — modM,
// matrix entries, etc.).
uint64_t parse_hex(const std::string& s, const std::string& filename) {
    try {
        size_t pos = 0;
        uint64_t v = std::stoull(s, &pos, 16);
        if (pos != s.size()) {
            throw std::runtime_error(
                "legacy_reader: trailing junk in hex '"
                + s + "' (in " + filename + ")");
        }
        return v;
    } catch (const std::invalid_argument&) {
        throw std::runtime_error(
            "legacy_reader: not a hex literal: '" + s
            + "' (in " + filename + ")");
    } catch (const std::out_of_range&) {
        throw std::runtime_error(
            "legacy_reader: hex literal out of range: '" + s
            + "' (in " + filename + ")");
    }
}

std::ifstream open_or_throw(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error(
            "legacy_reader: cannot open file: " + filename);
    }
    return in;
}

// ── Bounded token-stream cursor (for the Python "tokens" loops) ──────

struct Cursor {
    const std::vector<std::string>* toks;
    size_t i;
    const std::string* filename;

    const std::string& take() {
        if (i >= toks->size()) {
            throw std::runtime_error(
                "legacy_reader: unexpected end of file in " + *filename);
        }
        return (*toks)[i++];
    }
    int      take_int() { return parse_int(take(), *filename); }
    int64_t  take_dec() { return parse_dec(take(), *filename); }
    uint64_t take_hex() { return parse_hex(take(), *filename); }
};

// ── _read_hex_words: tempMK b/c reconstruction ────────────────────────

uint64_t read_hex_words(const std::vector<std::string>& tokens,
                        int nb_words, int w,
                        const std::string& filename) {
    assert(nb_words >= 1);
    // The Python helper does an arbitrary-precision shift; in C++ we
    // require the assembled value to fit in 64 bits. tempMK / tempMK2
    // in the legacy format always satisfy 32 * nb_words <= 64 because
    // w <= 32 ⇒ nb_words == 1 (single 32-bit word).
    if (32 * nb_words > 64) {
        throw std::runtime_error(
            "legacy_reader: hex words too wide (32*nb_words="
            + std::to_string(32 * nb_words)
            + " > 64) in " + filename);
    }
    if ((int)tokens.size() < nb_words) {
        throw std::runtime_error(
            "legacy_reader: expected " + std::to_string(nb_words)
            + " hex word(s), got " + std::to_string(tokens.size())
            + " in " + filename);
    }
    uint64_t val = 0;
    for (int i = 0; i < nb_words; ++i) {
        val = (val << 32) | parse_hex(tokens[i], filename);
    }
    int shift = 32 * nb_words - w;
    if (shift > 0) val >>= shift;
    return val;
}

// ── Per-tag generator readers ─────────────────────────────────────────
//
// Each returns a list of LegacyGeneratorSpec objects. Family names
// here are the canonical C++ class names accepted by
// create_generator() in factory.cpp:
//
//   polylcg                 → PolyLCG
//   taus / taus2            → Tausworthe
//   tgfsr                   → TGFSR
//   MT                      → MersenneTwister
//   genf2w (gen_type==1)    → GenF2wLFSR
//   genf2w (gen_type!=1)    → GenF2wPolyLCG
//   carry                   → WELLGen   (the old "Carry2Gen" file alias)
//   matsumoto               → MatsumotoGen
//   marsaxorshift           → MarsaXorshiftGen

std::vector<LegacyGeneratorSpec>
read_polylcg(std::ifstream& f, int L, const std::string& filename) {
    std::string tag_line;
    std::getline(f, tag_line);   // skip tag

    std::string line;
    if (!std::getline(f, line)) {
        throw std::runtime_error(
            "legacy_reader: missing nbgen line in " + filename);
    }
    int n;
    {
        std::istringstream ls(line);
        std::string ntok;
        if (!(ls >> ntok)) {
            throw std::runtime_error(
                "legacy_reader: empty nbgen line in " + filename);
        }
        n = parse_int(ntok, filename);
    }

    auto toks = tokenise_rest(f);
    Cursor cur{&toks, 0, &filename};

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(n);
    for (int gi = 0; gi < n; ++gi) {
        int k = cur.take_int();
        std::vector<int> exponents;
        while (true) {
            int e = cur.take_int();
            exponents.push_back(e);
            if (e == 0) break;
        }
        Params p;
        p.set_int("k", k);
        p.set_int_vec("poly", exponents);
        out.push_back({"PolyLCG", std::move(p)});
    }
    return out;
}

std::vector<LegacyGeneratorSpec>
read_tausworthe(std::ifstream& f, int L, const std::string& filename) {
    std::string line;
    std::getline(f, line);   // skip tag
    auto header = tokenise_line(f, filename);
    if (header.size() < 2) {
        throw std::runtime_error(
            "legacy_reader: tausworthe header needs >= 2 tokens in "
            + filename);
    }
    int n = parse_int(header[0], filename);
    bool quicktaus = parse_int(header[1], filename) != 0;
    int file_smax = 0;
    if (!quicktaus) {
        if (header.size() < 3) {
            throw std::runtime_error(
                "legacy_reader: tausworthe header needs smax token "
                "when not quicktaus, in " + filename);
        }
        file_smax = parse_int(header[2], filename);
    }
    auto toks = tokenise_rest(f);
    Cursor cur{&toks, 0, &filename};

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(n);
    for (int gi = 0; gi < n; ++gi) {
        std::vector<int> R;
        while (true) {
            int e = cur.take_int();
            R.push_back(e);
            if (e == 0) break;
        }
        std::vector<int> Q(R.rbegin(), R.rend());
        int nq = (int)Q.size();
        int k = Q[nq - 1];
        int s;
        if (!quicktaus) {
            s = std::min(file_smax, L - k);
            if (s < 1) s = 1;
        } else {
            s = k - Q[nq - 2];
        }
        Params p;
        p.set_int_vec("poly", Q);
        p.set_int("s", s);
        p.set_bool("quicktaus", quicktaus);
        out.push_back({"Tausworthe", std::move(p)});
    }
    return out;
}

std::vector<LegacyGeneratorSpec>
read_tgfsr(std::ifstream& f, int L, const std::string& filename) {
    std::string line;
    std::getline(f, line);   // skip tag
    auto header = tokenise_line(f, filename);
    if (header.size() < 2) {
        throw std::runtime_error(
            "legacy_reader: tgfsr header needs (w, r) in " + filename);
    }
    int w = parse_int(header[0], filename);
    int r = parse_int(header[1], filename);
    if (w > 32) {
        throw std::runtime_error(
            "legacy_reader: tgfsr w must be <= 32, got "
            + std::to_string(w) + " in " + filename);
    }
    auto count_line = tokenise_line(f, filename);
    int nbgen = parse_int(count_line[0], filename);

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(nbgen);
    for (int gi = 0; gi < nbgen; ++gi) {
        auto toks = tokenise_line(f, filename);
        if (toks.size() < 2) {
            throw std::runtime_error(
                "legacy_reader: tgfsr row needs (a_hex, m) in "
                + filename);
        }
        uint64_t a_val = parse_hex(toks[0], filename);
        int m_val = parse_int(toks[1], filename);
        Params p;
        p.set_int("w", w);
        p.set_int("r", r);
        p.set_int("m", m_val);
        p.set_int("a", static_cast<int64_t>(a_val));
        out.push_back({"TGFSR", std::move(p)});
    }
    return out;
}

std::vector<LegacyGeneratorSpec>
read_mt(std::ifstream& f, int L, const std::string& filename) {
    std::string line;
    std::getline(f, line);   // skip tag
    auto count_line = tokenise_line(f, filename);
    int nbgen = parse_int(count_line[0], filename);

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(nbgen);
    for (int gi = 0; gi < nbgen; ++gi) {
        auto toks = tokenise_line(f, filename);
        if (toks.size() < 5) {
            throw std::runtime_error(
                "legacy_reader: MT row needs (w, r, m, p, a_hex) in "
                + filename);
        }
        int w = parse_int(toks[0], filename);
        int r = parse_int(toks[1], filename);
        int m = parse_int(toks[2], filename);
        int p_v = parse_int(toks[3], filename);
        uint64_t a = parse_hex(toks[4], filename);
        if (w > 32) {
            throw std::runtime_error(
                "legacy_reader: MT w must be <= 32, got "
                + std::to_string(w) + " in " + filename);
        }
        Params pp;
        pp.set_int("w", w);
        pp.set_int("r", r);
        pp.set_int("m", m);
        pp.set_int("p", p_v);
        pp.set_int("a", static_cast<int64_t>(a));
        out.push_back({"MersenneTwister", std::move(pp)});
    }
    return out;
}

std::vector<LegacyGeneratorSpec>
read_genf2w(std::ifstream& f, int L, const std::string& filename) {
    std::string line;
    std::getline(f, line);   // skip tag
    int gen_type    = parse_int(tokenise_line(f, filename)[0], filename);
    int nb_token    = parse_int(tokenise_line(f, filename)[0], filename);
    bool normal_basis = (nb_token != 0);
    int nbgen       = parse_int(tokenise_line(f, filename)[0], filename);

    auto toks = tokenise_rest(f);
    Cursor cur{&toks, 0, &filename};

    const std::string fam = (gen_type == 1) ? "GenF2wLFSR" : "GenF2wPolyLCG";

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(nbgen);
    for (int gi = 0; gi < nbgen; ++gi) {
        int w        = cur.take_int();
        int r        = cur.take_int();
        uint64_t modM = cur.take_hex();
        int step     = cur.take_int();
        int nbcoeff  = cur.take_int();

        std::vector<uint64_t> coeff_vals;
        std::vector<int>      coeff_pos;
        coeff_vals.reserve(nbcoeff);
        coeff_pos.reserve(nbcoeff);
        for (int j = 0; j < nbcoeff; ++j) {
            uint64_t v = cur.take_hex();
            int pos = cur.take_int();
            coeff_vals.push_back(v);
            coeff_pos.push_back(pos);
        }

        Params p;
        p.set_int("w", w);
        p.set_int("r", r);
        p.set_int("modM", static_cast<int64_t>(modM));
        p.set_bool("normal_basis", normal_basis);
        p.set_int("step", step);
        // type kept for parity with the Python dict, although the
        // family name already encodes lfsr-vs-polylcg.
        p.set_string("type", (gen_type == 1) ? "lfsr" : "polylcg");
        // Flatten the Python list-of-{value,position}-dicts into the
        // two parallel vectors that F2wLFSRGen / F2wPolyLCGGen expect.
        p.set_uint_vec("coeff", coeff_vals);
        p.set_int_vec("nocoeff", coeff_pos);
        out.push_back({fam, std::move(p)});
    }
    return out;
}

std::vector<LegacyGeneratorSpec>
read_carry(std::ifstream& f, int L, const std::string& filename) {
    std::string line;
    std::getline(f, line);   // skip tag
    auto header = tokenise_line(f, filename);
    if (header.size() < 3) {
        throw std::runtime_error(
            "legacy_reader: carry header needs (w, r, p) in " + filename);
    }
    int w = parse_int(header[0], filename);
    int r = parse_int(header[1], filename);
    int p = parse_int(header[2], filename);
    int nbgen = parse_int(tokenise_line(f, filename)[0], filename);

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(nbgen);
    for (int gi = 0; gi < nbgen; ++gi) {
        auto toks = tokenise_line(f, filename);
        // Per Python: 3 m's + 8 × (1 type + 3 pi + 3 pu) = 3 + 8*7 = 59.
        if (toks.size() < 59) {
            throw std::runtime_error(
                "legacy_reader: carry row needs >=59 tokens, got "
                + std::to_string(toks.size()) + " in " + filename);
        }
        int m1 = parse_int(toks[0], filename);
        int m2 = parse_int(toks[1], filename);
        int m3 = parse_int(toks[2], filename);

        std::vector<int> mat_types;
        std::vector<int> mat_pi;
        std::vector<uint64_t> mat_pu;
        mat_types.reserve(8);
        mat_pi.reserve(24);
        mat_pu.reserve(24);
        size_t idx = 3;
        for (int k = 0; k < 8; ++k) {
            mat_types.push_back(parse_int(toks[idx++], filename));
            mat_pi.push_back(parse_int(toks[idx++], filename));
            mat_pi.push_back(parse_int(toks[idx++], filename));
            mat_pi.push_back(parse_int(toks[idx++], filename));
            mat_pu.push_back(parse_hex(toks[idx++], filename));
            mat_pu.push_back(parse_hex(toks[idx++], filename));
            mat_pu.push_back(parse_hex(toks[idx++], filename));
        }

        Params pp;
        pp.set_int("w", w);
        pp.set_int("r", r);
        pp.set_int("p", p);
        pp.set_int("m1", m1);
        pp.set_int("m2", m2);
        pp.set_int("m3", m3);
        pp.set_int_vec("mat_types", mat_types);
        pp.set_int_vec("mat_pi", mat_pi);
        pp.set_uint_vec("mat_pu", mat_pu);
        // Python uses Generator.create("Carry2Gen", ...); the Python
        // alias table maps "Carry2Gen" → "WELLGen" before dispatch.
        // Skip the indirection here.
        out.push_back({"WELLGen", std::move(pp)});
    }
    return out;
}

std::vector<LegacyGeneratorSpec>
read_matsumoto(std::ifstream& f, int L, const std::string& filename) {
    std::string line;
    std::getline(f, line);   // skip tag
    int nbgens = parse_int(tokenise_line(f, filename)[0], filename);

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(nbgens);
    for (int gi = 0; gi < nbgens; ++gi) {
        auto toks = tokenise_line(f, filename);
        size_t i = 0;
        auto need = [&](const char* what) {
            if (i >= toks.size()) {
                throw std::runtime_error(
                    std::string("legacy_reader: matsumoto row missing ")
                    + what + " in " + filename);
            }
        };
        need("type");      int gen_type = parse_int(toks[i++], filename);
        need("n");         int n = parse_int(toks[i++], filename);
        need("m");         int m = parse_int(toks[i++], filename);
        need("nbpi");      int nbpi = parse_int(toks[i++], filename);
        std::vector<int> paramsint;
        paramsint.reserve(nbpi);
        for (int j = 0; j < nbpi; ++j) {
            need("paramsint entry");
            paramsint.push_back(parse_int(toks[i++], filename));
        }
        need("nbpu");      int nbpu = parse_int(toks[i++], filename);
        std::vector<uint64_t> paramsunsigned;
        paramsunsigned.reserve(nbpu);
        for (int j = 0; j < nbpu; ++j) {
            need("paramsunsigned entry");
            paramsunsigned.push_back(parse_hex(toks[i++], filename));
        }

        Params p;
        p.set_int("type", gen_type);
        p.set_int("n", n);
        p.set_int("m", m);
        p.set_int_vec("paramsint", paramsint);
        p.set_uint_vec("paramsunsigned", paramsunsigned);
        out.push_back({"MatsumotoGen", std::move(p)});
    }
    return out;
}

std::vector<LegacyGeneratorSpec>
read_marsaxorshift(std::ifstream& f, int L, const std::string& filename) {
    std::string line;
    std::getline(f, line);   // skip tag
    auto header = tokenise_line(f, filename);
    if (header.size() < 3) {
        throw std::runtime_error(
            "legacy_reader: marsaxorshift header needs (nbgen, r, w) in "
            + filename);
    }
    int nbgen = parse_int(header[0], filename);
    int r = parse_int(header[1], filename);
    int w = parse_int(header[2], filename);
    (void)tokenise_line(f, filename);   // skip numberGenToAllocate

    auto push_t1 = [&](std::vector<LegacyGeneratorSpec>& dst,
                       int a, int b, int c) {
        Params p;
        p.set_int("type", 1);
        p.set_int("w", w);
        p.set_int("r", 1);
        p.set_int_vec("shifts", {a, b, c});
        dst.push_back({"MarsaXorshiftGen", std::move(p)});
    };
    auto push_t2x = [&](std::vector<LegacyGeneratorSpec>& dst,
                        int gen_type, int m_val,
                        const std::vector<int>& pv,
                        const std::vector<int>& qv) {
        Params p;
        p.set_int("type", gen_type);
        p.set_int("w", w);
        p.set_int("r", r);
        p.set_int("m", m_val);
        p.set_int_vec("p", pv);
        p.set_int_vec("q", qv);
        dst.push_back({"MarsaXorshiftGen", std::move(p)});
    };

    std::vector<LegacyGeneratorSpec> out;
    out.reserve(nbgen * 4);   // crude upper bound
    for (int gi = 0; gi < nbgen; ++gi) {
        auto toks = tokenise_line(f, filename);
        size_t i = 0;
        auto take = [&]() -> int {
            if (i >= toks.size()) {
                throw std::runtime_error(
                    "legacy_reader: marsaxorshift row truncated in "
                    + filename);
            }
            return parse_int(toks[i++], filename);
        };
        int gen_type = take();

        if (gen_type == 1) {
            int a = take();
            int b = take();
            int c = take();
            push_t1(out, a, b, c);
            push_t1(out, c, b, a);
            push_t1(out, -a, -b, -c);
            push_t1(out, a, -c, -b);

        } else if (gen_type >= 21 && gen_type <= 25) {
            int m_val = take();
            std::vector<int> op(3), oq(3);
            for (int j = 0; j < 3; ++j) op[j] = take();
            for (int j = 0; j < 3; ++j) oq[j] = take();

            std::vector<std::pair<std::vector<int>, std::vector<int>>> variants;
            if (gen_type == 21) {
                variants.push_back({op, oq});
                variants.push_back({{-op[0], -op[1],  op[2]},
                                    {-oq[0],  oq[1],  oq[2]}});
                variants.push_back({{-op[1], -op[0],  op[2]},
                                    {-oq[0],  oq[1],  oq[2]}});
            } else if (gen_type == 22) {
                variants.push_back({op, oq});
                variants.push_back({{-op[0],  op[1],  op[2]},
                                    {-oq[0], -oq[1],  oq[2]}});
                variants.push_back({{-op[0],  op[1],  op[2]},
                                    {-oq[1], -oq[0],  oq[2]}});
            } else if (gen_type == 23) {
                variants.push_back({op, oq});
                variants.push_back({{-op[0], -op[1],  op[2]},
                                    {-oq[0],  oq[1],  oq[2]}});
            } else if (gen_type == 24) {
                variants.push_back({op, oq});
                variants.push_back({{-op[0],  op[1],  op[2]},
                                    {-oq[0], -oq[1],  oq[2]}});
            } else { // 25
                variants.push_back({op, oq});
                variants.push_back({{ op[2],  op[1],  op[0]}, oq});
                variants.push_back({{-op[0], -op[1], -op[2]}, oq});
                variants.push_back({{ op[0],  op[2],  op[1]}, oq});
            }
            for (auto& [pv, qv] : variants) {
                push_t2x(out, gen_type, m_val, pv, qv);
            }

        } else if (gen_type == 3) {
            int n3 = take();
            std::vector<int> tap_pos, tap_shifts;
            tap_pos.reserve(n3);
            tap_shifts.reserve(n3);
            for (int j = 0; j < n3; ++j) {
                tap_pos.push_back(take());
                tap_shifts.push_back(take());
            }
            Params p;
            p.set_int("type", 3);
            p.set_int("w", w);
            p.set_int("r", r);
            p.set_int_vec("tap_positions", tap_pos);
            p.set_int_vec("tap_shifts", tap_shifts);
            out.push_back({"MarsaXorshiftGen", std::move(p)});

        } else if (gen_type == 4) {
            int r_val = take();
            int m_val = take();
            int a = take();
            int b = take();
            int c = take();
            int d = take();
            Params p;
            p.set_int("type", 4);
            p.set_int("w", w);
            p.set_int("r", r_val);
            p.set_int("m", m_val);
            p.set_int_vec("p", {a, b});
            p.set_int_vec("q", {c, d});
            out.push_back({"MarsaXorshiftGen", std::move(p)});

        } else if (gen_type == 100) {
            int nbmi = take();
            int nbxorshift = take();
            (void)nbxorshift;   // accepted but unused — Python ignores it
            std::vector<int> mi_pos, mi_counts, mi_shifts;
            for (int j = 0; j < nbmi; ++j) {
                mi_pos.push_back(take());
                int cnt = take();
                mi_counts.push_back(cnt);
                for (int s = 0; s < cnt; ++s) mi_shifts.push_back(take());
            }
            Params p;
            p.set_int("type", 100);
            p.set_int("w", w);
            p.set_int("r", r);
            p.set_int_vec("mi_positions", mi_pos);
            p.set_int_vec("mi_counts", mi_counts);
            p.set_int_vec("mi_shifts", mi_shifts);
            out.push_back({"MarsaXorshiftGen", std::move(p)});

        } else {
            throw std::runtime_error(
                "legacy_reader: unknown marsaxorshift type "
                + std::to_string(gen_type) + " in " + filename);
        }
    }
    return out;
}

}  // namespace

// ── Public dispatch entry points ──────────────────────────────────────

std::vector<LegacyGeneratorSpec>
read_generator_specs(const std::string& filename, int L) {
    auto in = open_or_throw(filename);

    // Peek the tag without consuming the first line — each per-tag
    // reader expects to skip its own header line.
    auto pos = in.tellg();
    std::string first_line;
    if (!std::getline(in, first_line)) {
        throw std::runtime_error(
            "legacy_reader: empty file: " + filename);
    }
    std::istringstream ls(first_line);
    std::string tag;
    if (!(ls >> tag)) {
        throw std::runtime_error(
            "legacy_reader: missing tag on line 1 of " + filename);
    }
    in.clear();
    in.seekg(pos);

    if (tag == "polylcg")        return read_polylcg(in, L, filename);
    if (tag == "taus")           return read_tausworthe(in, L, filename);
    if (tag == "taus2")          return read_tausworthe(in, L, filename);
    if (tag == "tgfsr")          return read_tgfsr(in, L, filename);
    if (tag == "MT")             return read_mt(in, L, filename);
    if (tag == "genf2w")         return read_genf2w(in, L, filename);
    if (tag == "carry")          return read_carry(in, L, filename);
    if (tag == "marsaxorshift")  return read_marsaxorshift(in, L, filename);
    if (tag == "matsumoto")      return read_matsumoto(in, L, filename);

    throw std::runtime_error(
        "legacy_reader: unknown generator tag '" + tag
        + "' in " + filename);
}

std::vector<std::unique_ptr<Generator>>
read_generators(const std::string& filename, int L) {
    auto specs = read_generator_specs(filename, L);
    std::vector<std::unique_ptr<Generator>> out;
    out.reserve(specs.size());
    for (auto& s : specs) {
        out.push_back(create_generator(s.family, s.params, L));
    }
    return out;
}

LegacyTransformationsResult
read_transformation_specs(const std::string& filename) {
    auto in = open_or_throw(filename);

    LegacyTransformationsResult result;
    auto count_line = tokenise_line(in, filename);
    int n = parse_int(count_line[0], filename);

    for (int gi = 0; gi < n; ++gi) {
        auto toks = tokenise_line(in, filename);
        if (toks.empty()) {
            throw std::runtime_error(
                "legacy_reader: empty transformation row in " + filename);
        }
        const std::string& type_str = toks[0];

        Params params;
        std::string trans_type;

        if (type_str == "permut") {
            if (toks.size() < 4) {
                throw std::runtime_error(
                    "legacy_reader: permut row needs (w, p, q) in "
                    + filename);
            }
            int w = parse_int(toks[1], filename);
            int p = parse_int(toks[2], filename);
            int q = parse_int(toks[3], filename);
            // Per the spec: convert -1 sentinels into 0 (the C++
            // PermutationTrans factory will randomise from 0 the same
            // way Python does via its random_specs path).
            if (p == -1) p = 0;
            if (q == -1) q = 0;
            params.set_int("w", w);
            params.set_int("p", p);
            params.set_int("q", q);
            trans_type = "permut";

        } else if (type_str == "tempMK"  || type_str == "tempMKopt"
                || type_str == "tempMK2" || type_str == "tempMK2opt") {
            if (toks.size() < 4) {
                throw std::runtime_error(
                    "legacy_reader: tempMK[2][opt] row needs at least "
                    "(w, eta, mu) in " + filename);
            }
            int w   = parse_int(toks[1], filename);
            int eta = parse_int(toks[2], filename);
            int mu  = parse_int(toks[3], filename);
            size_t idx = 4;

            const bool is_mk2 =
                (type_str == "tempMK2" || type_str == "tempMK2opt");
            const bool is_opt =
                (type_str == "tempMKopt" || type_str == "tempMK2opt");
            int u = 0, l = 0;
            if (is_mk2) {
                if (idx + 1 >= toks.size()) {
                    throw std::runtime_error(
                        "legacy_reader: tempMK2 row missing (u, l) in "
                        + filename);
                }
                u = parse_int(toks[idx++], filename);
                l = parse_int(toks[idx++], filename);
                trans_type = "tempMK2";
            } else {
                trans_type = "tempMK";
            }

            if (idx >= toks.size()) {
                throw std::runtime_error(
                    "legacy_reader: tempMK row missing nb_words in "
                    + filename);
            }
            int nb_words = parse_int(toks[idx++], filename);

            if (is_opt) {
                // skip disp_progress + limit_v (consumed but unused)
                if (idx + 1 >= toks.size()) {
                    throw std::runtime_error(
                        "legacy_reader: tempMK*opt row missing "
                        "(disp_progress, limit_v) in " + filename);
                }
                idx += 2;
            }

            params.set_int("w", w);
            params.set_int("eta", eta);
            params.set_int("mu", mu);
            params.set_int("u", u);
            params.set_int("l", l);

            if (nb_words != -1) {
                auto b_tokens = tokenise_line(in, filename);
                auto c_tokens = tokenise_line(in, filename);
                uint64_t b = read_hex_words(b_tokens, nb_words, w, filename);
                uint64_t c = read_hex_words(c_tokens, nb_words, w, filename);
                params.set_int("b", static_cast<int64_t>(b));
                params.set_int("c", static_cast<int64_t>(c));
            }
            // else: b/c omitted → randomized later by fill_params.

            if (type_str == "tempMKopt" || type_str == "tempMK2opt") {
                result.mk_opt = true;
            }

        } else {
            throw std::runtime_error(
                "legacy_reader: unknown transformation type '"
                + type_str + "' in " + filename);
        }

        result.specs.push_back({std::move(trans_type), std::move(params)});
    }
    return result;
}

LegacyTransformationsBuilt
read_transformations(const std::string& filename) {
    auto specs_result = read_transformation_specs(filename);
    LegacyTransformationsBuilt out;
    out.mk_opt = specs_result.mk_opt;
    out.transformations.reserve(specs_result.specs.size());
    for (auto& s : specs_result.specs) {
        out.transformations.push_back(
            create_transformation(s.trans_type, s.params));
    }
    return out;
}

}  // namespace regpoly_legacy
