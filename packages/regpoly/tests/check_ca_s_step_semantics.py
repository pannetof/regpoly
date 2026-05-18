"""Verify that CellularAutomataGen's time-spacing parameter `s` is
implemented correctly: each call to next() should advance the state by
exactly `s` single-step CA updates.

The implementation at cellular_automata.cpp:57-60 is:

    void next() {
        for (int i = 0; i < s_; ++i)
            next_one_step_();
    }

We verify this by direct simulation:

  (1) State-equivalence: for several s values and several step counts N,
      a CA with parameter s after N next() calls must equal the same CA
      with parameter 1 after s*N next() calls (same initial state, same
      k, same rule150_positions).

  (2) Output-sequence consistency: the bit-i output sampled every step
      from a CA at s=1 (then decimated by s) equals the bit-i output
      sampled every step from a CA at s.

  (3) Char-poly relationship: with full-period T at s=1 (k=31 primitive),
      the s-stepped CA has T^s as its transition matrix.  T^s remains
      primitive iff gcd(s, 2^k-1) = 1.  The recovered char poly degree
      must be k for all such s.
"""

from __future__ import annotations

import math
from regpoly_cpp._regpoly_cpp import BitVect
from regpoly.core.generator import Generator


def state_to_int(g):
    """Read the current state as an integer (bit i of int = cell i)."""
    st = g._cpp_gen.state()
    k = g.k
    n = 0
    for i in range(k):
        if st.get_bit(i):
            n |= 1 << (k - 1 - i)
    return n


def make_ca(k, positions, s):
    return Generator.create("CellularAutomataGen",
                            L=min(k, 64), k=k, rule150_positions=positions, s=s)


def init_with_bit0(g):
    init = BitVect(g.k)
    init.set_bit(0, 1)
    g._cpp_gen.init(init)


def check_state_equivalence():
    print("=" * 72)
    print("Test (1): state equivalence — CA(s) after N steps == CA(1) after s*N steps")
    print("=" * 72)
    cases = [
        (31, [10],      [1, 2, 3, 5, 7, 8, 10]),
        (32, [0, 14],   [1, 2, 4, 7, 8, 10]),
        (65, [0],       [1, 5, 10]),                # Adak-Das CA(90'), k=65
        (35, list(range(1, 35)), [1, 7, 10]),       # Adak-Das CA(150'), k=35
    ]
    N = 20  # number of next() calls under the s-stepped CA

    all_ok = True
    for k, positions, s_values in cases:
        for s in s_values:
            g_s   = make_ca(k, positions, s=s)
            g_1   = make_ca(k, positions, s=1)
            init_with_bit0(g_s)
            init_with_bit0(g_1)
            for _ in range(N):
                g_s._cpp_gen.next()
            for _ in range(s * N):
                g_1._cpp_gen.next()
            a = state_to_int(g_s)
            b = state_to_int(g_1)
            ok = (a == b)
            all_ok &= ok
            tag = "OK " if ok else "FAIL"
            print(f"  [{tag}] k={k:3d} positions={str(positions)[:25]:25s} s={s:2d}  "
                  f"after N={N} steps:  CA(s) state == CA(1)^(s·N) state ? {ok}")
    return all_ok


def check_output_decimation():
    print()
    print("=" * 72)
    print("Test (2): output decimation — bit i of CA(s) outputs == every s-th bit i of CA(1)")
    print("=" * 72)
    cases = [
        (31, [10],    7),   # Table 9
        (32, [0, 14], 8),   # Table 10
        (32, [0, 14], 4),   # Table 8
        (31, [10],    2),
    ]
    N = 100  # number of CA(s) samples

    all_ok = True
    for k, positions, s in cases:
        g_s = make_ca(k, positions, s=s)
        g_1 = make_ca(k, positions, s=1)
        init_with_bit0(g_s)
        init_with_bit0(g_1)
        seq_s = []
        seq_1 = []
        for _ in range(N):
            g_s._cpp_gen.next()
            seq_s.append(g_s._cpp_gen.get_output().get_bit(0))
        for _ in range(s * N):
            g_1._cpp_gen.next()
            seq_1.append(g_1._cpp_gen.get_output().get_bit(0))
        # Decimate seq_1 to every s-th value: indices s-1, 2s-1, 3s-1, ...
        # (because after s next() calls on g_1, its state equals g_s after 1 next()).
        decimated = seq_1[s - 1::s]
        ok = (seq_s == decimated)
        all_ok &= ok
        tag = "OK " if ok else "FAIL"
        print(f"  [{tag}] k={k} positions={positions} s={s}: "
              f"seq_s[0:{N}] == seq_1[s-1::s][:{N}] ? {ok}")
        if not ok:
            # Find first mismatch.
            for i, (a, b) in enumerate(zip(seq_s, decimated)):
                if a != b:
                    print(f"        first mismatch at i={i}: seq_s={a} decimated={b}")
                    break
    return all_ok


def check_primitivity_under_s():
    print()
    print("=" * 72)
    print("Test (3): primitivity under T^s — gcd(s, 2^k-1) = 1 ⇒ CA(s) primitive iff CA(1) primitive")
    print("=" * 72)
    cases = [
        (31, [10]),
        (32, [0, 14]),
    ]
    s_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    all_ok = True
    for k, positions in cases:
        g1 = make_ca(k, positions, s=1)
        prim_1 = g1._cpp_gen.is_full_period()
        print(f"  k={k} positions={positions}: CA(s=1) primitive = {prim_1}")
        rho = (1 << k) - 1
        for s in s_values:
            g_s = make_ca(k, positions, s=s)
            prim_s = g_s._cpp_gen.is_full_period()
            gcd = math.gcd(s, rho)
            expected = (gcd == 1) and prim_1
            ok = (prim_s == expected)
            all_ok &= ok
            tag = "OK " if ok else "FAIL"
            print(f"    [{tag}] s={s:2d}  gcd(s, 2^k-1)={gcd:>6d}  "
                  f"primitive(T^s) = {prim_s}  (expected {expected})")
    return all_ok


if __name__ == "__main__":
    r1 = check_state_equivalence()
    r2 = check_output_decimation()
    r3 = check_primitivity_under_s()
    print()
    print("=" * 72)
    print(f"Overall: state_equiv={r1}, output_decimation={r2}, primitivity_under_s={r3}")
    print(f"All pass: {r1 and r2 and r3}")
    print("=" * 72)
