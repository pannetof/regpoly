"""
generate_factors.py — Generate the primitive_factors.json data file.

Computes the prime factorization of Φ_k(2) (the k-th cyclotomic
polynomial evaluated at 2) for k = 2..N using:
  - Cunningham tables (2^n ± 1 factorizations)
  - FactorDB (for 2^a * m cases with a >= 2)
  - SymPy (for small numbers)

Output: src/regpoly/data/primitive_factors.json

Usage:
    python -m regpoly.tools.generate_factors [max_k]
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.request
import urllib.parse

from sympy import divisors, isprime
from sympy.ntheory import factorint as sym_factorint
from sympy.functions.combinatorial.numbers import mobius

URL_MINUS = "https://homes.cerias.purdue.edu/~ssw/cun/pmain126.txt"


# ═══════════════════════════════════════════════════════════════════════════
# Cunningham table parsing
# ═══════════════════════════════════════════════════════════════════════════

def fetch(url):
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode("ascii", errors="replace")


def parse_block(text, header):
    start = text.find(header)
    if start == -1:
        return {}
    next_table = text.find("Table", start + len(header))
    chunk = text[start: next_table if next_table != -1 else None]

    result = {}
    entry_re = re.compile(r'^\s{0,6}(\d+)\s+(.*)')
    cont_re = re.compile(r'^\s{20,}(.*)')
    current_n, current_tail = None, ""

    def flush(n, tail):
        if n is None:
            return
        tail2 = re.sub(r'\([^)]*\)', '', tail).strip()
        tokens = [t.strip().rstrip('*') for t in tail2.split('.') if t.strip()]
        primes, complete = [], True
        for tok in tokens:
            if re.fullmatch(r'\d+', tok):
                primes.append(int(tok))
            elif re.fullmatch(r'[Pp]\d+', tok):
                complete = False
        result[n] = {"primes": primes, "complete": complete}

    for line in chunk.splitlines():
        m = entry_re.match(line)
        if m:
            flush(current_n, current_tail)
            current_n, current_tail = int(m.group(1)), m.group(2)
        else:
            c = cont_re.match(line)
            if c and current_n is not None:
                current_tail += c.group(1)
    flush(current_n, current_tail)
    return result


def load_tables():
    print("Fetching Cunningham tables...", end=" ", flush=True)
    text = fetch(URL_MINUS)
    minus = parse_block(text, "Table 2-")
    plus = parse_block(text, "Table 2+")
    print(f"done ({len(minus)} entries in 2-, {len(plus)} entries in 2+)")
    return minus, plus


# ═══════════════════════════════════════════════════════════════════════════
# Cyclotomic polynomial Φ_k(2)
# ═══════════════════════════════════════════════════════════════════════════

def phi2(n):
    num, den = 1, 1
    for d in divisors(n):
        mu = mobius(n // d)
        if mu == 1:
            num *= (2 ** d - 1)
        elif mu == -1:
            den *= (2 ** d - 1)
    return num // den


# ═══════════════════════════════════════════════════════════════════════════
# FactorDB
# ═══════════════════════════════════════════════════════════════════════════

def factordb_factors(n):
    try:
        url = f"http://factordb.com/api?query={n}"
        with urllib.request.urlopen(url, timeout=10) as r:
            data = json.load(r)
        status = data.get("status", "U")
        if status in ("FF", "P", "Prp"):
            primes = [int(p) for p, _ in data["factors"]]
            return primes, True
        elif status == "CF":
            primes = [int(p) for p, _ in data["factors"]
                      if isprime(int(p))]
            return primes, False
    except Exception:
        pass
    return [], False


# ═══════════════════════════════════════════════════════════════════════════
# Factor Φ_k(2)
# ═══════════════════════════════════════════════════════════════════════════

def primitive_prime_factors(k, minus_table, plus_table):
    a, m = 0, k
    while m % 2 == 0:
        a += 1
        m //= 2

    if a == 0:
        if k not in minus_table:
            return None, False
        entry = minus_table[k]
        return entry["primes"], entry["complete"]

    elif a == 1:
        if m not in plus_table:
            return None, False
        entry = plus_table[m]
        return entry["primes"], entry["complete"]

    else:
        val = phi2(k)
        if val == 1:
            return [], True
        if isprime(val):
            return [int(val)], True
        primes, complete = factordb_factors(val)
        if complete:
            return primes, True
        if val.bit_length() <= 100:
            fdict = sym_factorint(int(val))
            return sorted(fdict.keys()), True
        return primes, complete


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def generate(max_k: int = 1499):
    minus, plus = load_tables()

    data = {}
    incomplete = 0
    missing = 0

    for k in range(2, max_k + 1):
        primes, complete = primitive_prime_factors(k, minus, plus)

        if primes is None:
            data[str(k)] = {"factors": [], "complete": False}
            missing += 1
        else:
            data[str(k)] = {
                "factors": [int(p) for p in primes],
                "complete": complete,
            }
            if not complete:
                incomplete += 1

        if k % 100 == 0:
            print(f"  k={k:5d}  ({incomplete} incomplete, {missing} missing)",
                  flush=True)

    # Rate-limit FactorDB calls (be polite)
    time.sleep(0.1)

    print(f"\nDone: k=2..{max_k}")
    print(f"  Complete: {max_k - 1 - incomplete - missing}")
    print(f"  Incomplete: {incomplete}")
    print(f"  Missing: {missing}")

    return data


def main():
    max_k = int(sys.argv[1]) if len(sys.argv) > 1 else 1499

    data = generate(max_k)

    out_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "primitive_factors.json")

    with open(out_path, "w") as f:
        json.dump(data, f, separators=(",", ":"))

    size_kb = os.path.getsize(out_path) / 1024
    print(f"\nWritten to {out_path} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
