# Cellular Automata (Bhuvaneswari–Bhattacharjee)

**C++ class:** `CellularAutomataGen`
**Name in code:** `"CellularAutomataGen"`
**Legacy aliases:** `CA`

A k-cell linear maximal-length Elementary Cellular Automaton (ECA) over
the null boundary, using only rules 90 and 150.  The family encodes the
single-component CAs used by Bhuvaneswari & Bhattacharjee (2026) to
build their 2-component combined PRNGs with time-spacing.  Two such
components are combined via the existing `CombinedGenerator` to form
the full paper-style generator (no new combined-family class needed).

## Mathematical recurrence

The state is a packed $k$-bit `BitVect` with cell 0 stored at the MSB and
cell $k-1$ at the LSB.  Each cell updates under one of two rules:

- **Rule 90** on cell $i$: $S_i' = S_{i-1} \oplus S_{i+1}$
- **Rule 150** on cell $i$: $S_i' = S_{i-1} \oplus S_i \oplus S_{i+1}$

Outside-boundary cells ($i = -1$, $i = k$) are pinned to 0 (null boundary).

Stacking the per-cell update over all $k$ cells, one CA step collapses
to three bitvector operations:

```
state ← (state << 1) XOR (state >> 1) XOR (rules AND state)
        ^ brings cell i+1      ^ brings cell i-1   ^ middle bit only
        ^ to cell i            ^ to cell i         ^ at rule-150 cells
```

where `rules` is the $k$-bit vector with bit $i = 1$ iff cell $i$ uses rule 150.
`next()` applies this $s$ times in a row, where $s$ is the time-spacing.
The matrix $T^s$ is never materialized — the equidistribution machinery
extracts it lazily via the base `Generator::transition_matrix()` helper.

### Characteristic polynomial

Recovered at runtime by the inherited Berlekamp–Massey routine
(`Generator::char_poly()`).  For an isolated CA at $s = 1$ the resulting
polynomial is primitive iff the rule vector is one of the published
maximal-length CA constructions — see Cattell & Zhang 1995 (Table 2)
or Adak & Das 2021 (CA(90′) / CA(150′) families).  At time-spacing $s$,
the recovered polynomial is the characteristic polynomial of $T^s$,
which is primitive iff $T$ is primitive AND $\gcd(s, 2^k - 1) = 1$.

## Parameters

| Name                | Type      | Role       | Rand type | Description |
|---------------------|-----------|------------|-----------|-------------|
| `k`                 | `int`     | structural | --        | Number of cells. Defines the state size. |
| `rule150_positions` | `int_vec` | structural | --        | 0-indexed cell positions using rule 150; all others use rule 90. |
| `s`                 | `int`     | structural | --        | Time-spacing (default 1). Each `next()` advances by `s` CA steps. |

There are no non-structural (search) parameters — the rule vector and
time-spacing are pinned by the paper-published constructions.

## State size

- $k = $ number of cells (paper uses $31 \le k \le 128$ except for R2 where $k = 1409$)
- period $= 2^k - 1$ when $T$ is primitive AND $\gcd(s, 2^k - 1) = 1$
- period $= (2^k - 1) / \gcd(s, 2^k - 1)$ otherwise

## Catalog entries

The full set of paper-published constructions is at
[`library/bhuvaneswari-bhattacharjee-2026-cellular-automata`](../papers/bhuvaneswari-bhattacharjee-2026-cellular-automata.md).
That entry catalogs ~116 generators:

- **R1** (k=32 single CA, weak — Table 1).
- **(31, 32) at s ∈ {1, 2, 4, 7, 8}** — Tables 3, 7, 8, 9, 10.
- **49 combined Cattell–Zhang CA-PRNGs** — Table 12.
- **~60 Adak–Das CA(90′)/CA(150′) combinations** — Table 11.
- **R3, R4, R2 weak generators** — Section 3.1 / Table 1 baseline.

Rule vectors for the Cattell–Zhang component CAs come from
`packages/regpoly/src/regpoly/data/cellular_automata.json` (transcribed
from paper Table 2; convention locked by Example 1).  Adak–Das
CA(90′)/CA(150′) constructions are derived programmatically.

## Implementation notes

- **Bit numbering.** `BitVect` stores bit 0 at the MSB; the recurrence
  shifts maintain this convention.  Paper Table 2 lists 1-indexed cell
  positions; we store 0-indexed positions throughout (subtract 1).
- **No matrix materialization.** The shift–XOR recurrence is `O(s · ⌈k/64⌉)`
  word operations per `next()`.  No `T^s` matrix is stored.
- **Combined two-component PRNGs** are built via `CombinedGenerator`
  with two `CellularAutomataGen` components.  Each component carries
  its own time-spacing s; the combined output is the bitwise XOR.
- **Output convention is MSB-aligned by default** — cell 0 of each
  component lives at the most significant output bit; if k < L, the
  low-order bits are zero-padded.  The paper hedges between "left or
  right padding"; both conventions give similar rank deficiency
  patterns (see paper-disagreement note below).

## Two conventions for the equidistribution matrix $B$

The paper's per-$t$ rank values in Tables 1, 3, 7-10 are only
reproducible under a **non-standard** matrix construction.  The
standard ME definition asks "$\mathrm{rank}(B) = t \cdot \ell$ where $B$ has $t \cdot \ell$ rows
from $t$ advances of the recurrence" — that is what our C++
matricial-DE kernel implements.  The paper uses:

> **Paper convention.** $B$ has $(t+1) \cdot \ell$ rows: block 0 is the
> *initial-state* output (no advance), blocks $1, \ldots, t$ come from $t$
> advances.  $(t, \ell)$-equidistributed iff $\mathrm{rank}(B) \ge t \cdot \ell$.

This convention is identified empirically.  For R1 (single 32-bit
CA, paper Table 1) the paper-convention rank reproduces 8 of 9
published ranks exactly; the remaining row at $t = 16$ is a likely
paper typo.  For combined $(31, 32, s \in \{1, 2, 4, 7, 8\})$ (paper
Tables 3, 7, 8, 9, 10) the paper-convention binary verdict matches
the paper at **every per-row cell**.

## Reproduction of paper Tables under paper convention

| Paper data | Reproduction |
|---|---|
| Table 1 — R1 single CA  | 8/9 per-row ranks match; $t = 16$ likely typo |
| Table 2 — Cattell-Zhang positions, $k = 29\ldots128$ | 100/100 primitive ✓ |
| Table 3 — combined $(31, 32)$, $s = 1$, NOT ME       | 13/13 binary verdicts match ✓ |
| Table 7 — combined $(31, 32)$, $s = 2$, NOT ME       | matches at the ME-overall level ✓ |
| Table 8 — combined $(31, 32)$, $s = 4$, NOT ME       | matches at the ME-overall level ✓ |
| Table 9 — combined $(31, 32)$, $s = 7$, ME           | 13/13 binary verdicts match (all equi) ✓ |
| Table 10 — combined $(31, 32)$, $s = 8$, ME          | 13/13 binary verdicts match (all equi) ✓ |
| Table 11 — Adak-Das CA(90′)/CA(150′) combos  | primitivity verified; ME claims tracked |
| Table 12 — 49 combined CA-PRNGs              | 35/49 match paper-conv ME; 14/49 disagree |
| Section 3.1 — R1, R2, R3, R4 (weak)          | all primitive, all NOT ME ✓ |

## 14 Table 12 entries that disagree with the paper

Under paper convention, these entries have specific $(t, \ell^\ast)$ cells
where $\mathrm{rank}(B) < t \cdot \ell$, contradicting the paper's "ME" claim:

| Entry | Deficit row(s) — `(t, l*, our_rank, t·l)` |
|---|---|
| (47, 56, s=8)  | (3, 32, 93, 96), (4, 25, 98, 100) |
| (31, 32, s=5)  | (3, 21, 62, 63) |
| (33, 40, s=10) | (36, 2, 70, 72) |
| (35, 64, s=10) | (3, 32, 95, 96), (33, 3, 98, 99) |
| (39, 56, s=9)  | (47, 2, 93, 94) |
| (41, 64, s=10) | (3, 32, 95, 96) |
| (45, 64, s=9)  | (3, 32, 94, 96), (4, 27, 107, 108) |
| (45, 64, s=10) | (54, 2, 107, 108) |
| (47, 56, s=9)  | (51, 2, 100, 102) |
| (47, 72, s=10) | (4, 29, 115, 116) |
| (59, 64, s=10) | (41, 3, 122, 123) |
| **(63, 64, s=10)** | (4, 31, 115, 124) — large deficit |
| **(63, 80, s=10)** | (4, 32, 115, 128), (5, 28, 128, 140) — large deficit |
| **(67, 72, s=10)** | (4, 32, 122, 128) — large deficit |

11 of these have small deficits (1–3 ranks short).  The remaining 3
have substantial deficits (6, 9, 12+).  A block-count study shows
that $(t+2) \cdot \ell$ blocks would close 12/14 and $(t+3) \cdot \ell$ all 14 — but
those block counts over-shoot the R1 Table 1 ranks, so they cannot
be the universal convention.  Most plausible reading: the paper
made arithmetic errors or applied "almost ME" labelling for these
specific entries.

## Kernel convention vs paper convention

Our C++ matricial-DE kernel implements the **standard** convention
($t \cdot \ell$ rows, equi iff $\mathrm{rank} = t \cdot \ell$).  For combined generators where
the conventions diverge (e.g., $(31, 32, s = 7)$), the kernel reports
NOT ME even though the paper convention (and the paper) say ME.  The
period and primitivity machinery is unaffected — those are reproduced
byte-for-byte.

See
`packages/regpoly/tests/test_cellular_automata_paper_disagreement.py`
for the ground-truth GF(2) rank computations under both conventions,
the R1 Table 1 reproduction, and the 14 enumerated Table 12
disagreement deficits.

---

## References

- Bhuvaneswari A & Bhattacharjee K, 2026.
  *Cellular Automata based Resource Efficient Maximally Equidistributed
  Pseudo-Random Number Generators.*
  arXiv preprint 2603.19656.
  [DOI](https://doi.org/10.48550/arXiv.2603.19656).
  Published as
  [`bhuvaneswari-bhattacharjee-2026-cellular-automata`](../papers/bhuvaneswari-bhattacharjee-2026-cellular-automata.md)
  in the library.
- Cattell K & Zhang S, 1995.
  *Minimal cost one-dimensional linear hybrid cellular automata of
  degree through 500.*
  Journal of Electronic Testing 6(2), 255–258.
  Source of the maximal-length CA rule vectors used by
  `CellularAutomataGen` for k = 29..128.
- Adak S & Das S, 2021.
  *(Imperfect) strategies to generate primitive polynomials over GF(2).*
  Theoretical Computer Science 872, 79–96.
  Source of the CA(90′) / CA(150′) constructions used at Sophie-Germain
  prime sizes (k ∈ {26, 29, 35, 39, 65, 69, 105, 113, 119}).
