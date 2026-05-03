# Tempering Optimization Algorithm

Optimization of tempering bitmask parameters to minimize the total
dimension defect of an F2-linear generator.  The algorithm uses
**safe masks** with **random multi-bit perturbation** and an
**incremental dual-lattice basis** (StackBase), mirroring the C
code's `OptimizeTemper()` from `temperMK.c`.

## 1. Background

### Why tempering?

An F2-linear generator produces output by extracting $w$ bits from a
$k$-bit state.  Without tempering, the raw output often has poor
equidistribution: certain bit patterns are overrepresented or
underrepresented in tuples of consecutive outputs.

A **tempering transformation** $T$ is a bijective GF(2)-linear map
applied to the output word before it is returned to the user.
Because $T$ is bijective, it does not change the generator's period
or its state-update structure.  It only rearranges the output bits.
The goal of tempering optimization is to choose $T$'s free parameters
(bitmasks $b$, $c$, etc.) so that the equidistribution is as good
as possible.

### Equidistribution

The **dimension of equidistribution with $v$-bit accuracy** is $t_v$:
the largest $t$ such that the $t$-tuples $(u_i, \ldots, u_{i+t-1})$,
each truncated to their $v$ most significant bits, cover all $2^{tv}$
patterns equally over one full period.  The upper bound is
$t_v^* = \lfloor k / v \rfloor$.

The **dimension gap** at resolution $v$ is:

$$d(v) = t_v^* - t_v = \lfloor k/v \rfloor - t_v, \qquad v = 1, \ldots, L.$$

The **total dimension defect** is:

$$\Delta = \sum_{v=1}^{L} d(v).$$

A generator is **maximally equidistributed (ME)** when $\Delta = 0$,
meaning $d(v) = 0$ for every resolution.

### Tempering types

Two tempering families are supported:

**Matsumoto-Kurita (MK) tempering** — applies shift-and-mask operations to
the output word.  Type I:

$$y \leftarrow y \oplus \bigl((y \ll \eta) \;\&\; b\bigr), \qquad
  y \leftarrow y \oplus \bigl((y \ll \mu) \;\&\; c\bigr)$$

Type II (MT19937-style) adds right-shift steps before and after.
The free parameters are the bitmasks $b$ and $c$.

**Lagged tempering** (MELG-style) — uses a lookahead word from the
state array:

$$y \leftarrow y \oplus (y \ll \sigma), \qquad
  y \leftarrow y \oplus (w_{i+L} \;\&\; b)$$

The free parameter is the bitmask $b$; the structural parameters
$\sigma$ and $L$ are fixed before optimization.


## 2. Safe Masks

### The problem

At resolution $v$, the dual-lattice method computes $d(v)$ using the
generating polynomials $h_0, h_1, \ldots, h_{v-1}$.  The polynomial
$h_j$ depends on the $j$-th output bit of the tempered generator.

If we perturb a bitmask bit that affects output bits $0, \ldots, v-2$,
the polynomials $h_0, \ldots, h_{v-2}$ change, and the cached
lattice basis for resolutions $1, \ldots, v-1$ becomes invalid.
To avoid a full rebuild, we must restrict perturbations at resolution
$v$ to bits that **only** affect output bits $v-1, v, \ldots, w-1$.

### Definition

The **safe mask** for parameter $p$ at resolution $v$ is a bitmask
$S_v(p)$ where bit $j$ is set iff flipping bit $j$ of parameter $p$
does not change any of the output bits $0, 1, \ldots, v-2$.

For a $w$-bit parameter, the base safe mask at resolution $v$ is:

$$S_v = (1 \ll (w - v + 1)) - 1$$

This allows only the lowest $w - v + 1$ bits to be perturbed.  These
correspond to the bit positions $v-1, v, \ldots, w-1$ of the output
(positions that do not affect the already-resolved resolutions).

### MK tempering adjustment

For MK type I tempering, the $b$-step feeds into the $c$-step via
the left shift by $\mu$.  A change to $b[p]$ propagates through the
$c$-step to output bit $p + \mu$.  If $p + \mu < v - 1$, the cached
basis is invalidated.  The safe mask for parameter $b$ therefore
clears additional bits:

```
for p in range(mu, min(v + mu - 1, w)):
    safe_mask_b[v] &= ~(1 << (w - 1 - p))
```

This ensures that no perturbation of $b$ propagates through $c$ to
an output bit below position $v - 1$.

### Properties

- At $v = 1$: the safe mask is $(1 \ll w) - 1$ (all bits safe), because
  there are no lower resolutions to protect.
- At $v = L$: only one bit is safe (the LSB), severely constraining
  the search.
- The number of safe bits decreases monotonically with $v$.


## 3. Incremental Lattice Basis (StackBase)

Computing $d(v)$ from scratch requires building the full dual-lattice
basis and running Lenstra's reduction — cost $O(v^2 \cdot k / 64)$ per
resolution.  For optimization, we need to evaluate $d(v)$ hundreds
of times at each resolution with different bitmask values.

The **StackBase** mechanism (implemented in `LatticeOptCache`)
maintains an incremental lattice basis with snapshots, enabling
O(1)-resolution recomputation:

### Operations

**`reset_step()`** — Initialize the incremental state.  Sets
`previous_v = 0` and clears the snapshot stack.

**`step(v)`** — Compute $d(v)$ for the current bitmask parameters.
Three cases depending on the relationship between `v` and `previous_v`:

1. **$v = 1$ (restart):**
   Rebuild from scratch.  Recompute $g_0^{-1} \bmod P$
   (`refresh_inv_g0`), recompute the generating polynomial $h_0$
   for the current bitmask, initialize the dual-lattice basis
   with the single row $(P(z))$, and run Lenstra at resolution 1.

2. **$v > \text{previous\_v}$ (forward):**
   Save a snapshot of the current basis at position $v-1$.
   Recompute only the generating polynomial $h_{v-1}$ for the
   current bitmask.  Call `dual_base_increase` to extend the
   basis with the new row $(h_{v-1}, 0, \ldots, 0, 1)$.  Run
   Lenstra at resolution $v$.

3. **$v \le \text{previous\_v}$ (retry):**
   Restore the snapshot saved at position $v-1$ (reverting the
   basis to its state before $h_{v-1}$ was added).  Recompute
   $h_{v-1}$ for the new bitmask, extend the basis, and run
   Lenstra at resolution $v$.

### Why this is efficient

Each call to `step(v)` recomputes **one** generating polynomial
(cost: $O(k)$ generator steps) and runs **one** Lenstra reduction
(cost: $O(v^2 \cdot k / 64)$).  The previously resolved resolutions
$1, \ldots, v-1$ are not recomputed — their correctness is
guaranteed by the safe mask constraint.

The snapshot stack uses copy-on-save semantics: `opt_stack[v-2]`
holds the basis state just before resolution $v$ was added, so
retrying at $v$ after a perturbation only needs to restore from the
snapshot and redo the single-resolution computation.


## 4. The Recursive Optimization Algorithm

### Single pass: `_run_once(delta, mse)`

The core algorithm is a recursive descent over resolutions $v = 1, \ldots, L$.
At each level $v$, it tries several random multi-bit perturbations within
the safe mask, evaluates $d(v)$ via `step(v)`, and recurses deeper if the
gap is within tolerance.

```
function optimize(v):
    Lim := 5 if v < L/2 else 2       // more tries at low resolutions

    repeat Lim times (or until budget exhausted):
        essais := essais + 1

        // Perturb all bitmask parameters within their safe masks
        for each parameter p with width w:
            mask := safe_mask[v][p]
            if mask == 0: skip
            perturbation := random_bits(w) & mask
            p := p XOR perturbation

        // Evaluate
        d(v) := cache.step(v)
        se   := d(1) + d(2) + ... + d(v)

        // Track best solution seen so far
        if se < best_se[v] and v >= max_v:
            essais := 0                // reset budget on improvement
            save best parameters
            best_se[v] := se
            max_v := v

        // Recurse if within tolerance
        if d(v) <= delta[v] and se <= mse:
            if v >= L: break           // reached deepest level
            else: optimize(v + 1)      // go deeper

        // Early exit if target met at level L
        if d(L) >= 0 and d(L) <= delta[L] and best_se[L] <= mse:
            break
```

**Key design choices:**

- **Multi-bit perturbation:** Unlike the backtracking algorithm (which
  flips one bit at a time), this algorithm XORs a random vector masked
  to the safe bits.  This explores the parameter space more aggressively
  and is better suited to large bitmask parameters.

- **Budget with reset:** The global counter `essais` limits total work.
  It resets to 0 whenever a new best is found, giving the algorithm more
  time to explore promising regions.

- **Asymmetric effort:** Low resolutions ($v < L/2$) get 5 tries per
  recursion level; high resolutions get 2.  Low resolutions have larger
  safe masks and more freedom, so extra tries are more productive.

- **Greedy depth tracking:** `max_v` tracks the deepest resolution
  reached with improvement.  A solution at depth $v$ is only accepted
  as "best" if $v \ge \text{max\_v}$, preventing a shallow solution
  from overwriting a deeper (more complete) one.


## 5. Operating Modes

The optimizer supports three modes, selected by the arguments to
`run(delta, mse, n_restarts)`:

### Mode 1: ME search (`delta=None, mse=None, n_restarts=1`)

Sets `delta = [0, 0, ..., 0]` and `mse = 0`.  The recursion only
advances past resolution $v$ when $d(v) = 0$, and the overall
target is $\Delta = 0$.  This is the strictest mode: it searches for
a **maximally equidistributed** tempering.

### Mode 2: Target search (`delta=[...], mse=N, n_restarts=1`)

The user provides per-resolution gap bounds `delta[v]` and a total
defect bound `mse`.  The recursion advances when $d(v) \le \delta[v]$
and $\text{se} \le \text{mse}$.  This is useful for finding
"almost ME" solutions when exact ME is too expensive or impossible.

### Mode 3: Minimize (`n_restarts > 1`)

Iteratively tightens the target.  Ignores user-supplied `delta` and
`mse`.

```
function run_minimize(n_restarts):
    best_result := None
    failures := 0

    while failures < n_restarts:
        randomize all bitmask parameters

        if best_result is None:
            delta := [INT_MAX, ...], mse := INT_MAX    // first pass: unconstrained
        else:
            delta := best_result.gaps                   // tighten to current best
            mse   := best_result.se

        result := run_once(delta, mse)

        if result.se < best_result.se:
            best_result := result
            failures := 0                               // reset on improvement
        else:
            failures := failures + 1                    // count consecutive failures

    return best_result
```

The first pass with `delta = INT_MAX` runs a single deep sweep that
finds any solution.  Subsequent passes use the best gaps found so
far as the target, forcing the algorithm to find strictly better
solutions.  The search stops after `n_restarts` consecutive failures
to improve.


## 6. Putting It All Together

### API

The `TemperingOptimizer` is a reusable configuration object.
Create it once with settings, then call `run(comb)` on any
`Combinaison`:

```python
opt = TemperingOptimizer(max_essais=400, delta=[0]*33, mse=0)
result = opt.run(comb)       # ME mode
result = opt.run(other_comb) # reuse same settings
```

### What `run(comb)` does

1. Enumerate all optimizable bitmask parameters from the transformation
   chain: `(component_index, transform_index, param_name) -> bit_width`.
2. Compute safe masks analytically for each parameter at each resolution.
3. Build the `LatticeOptCache` from the C++ generator and transformation
   objects.
4. `reset_step()` initializes the incremental lattice state.
5. `optimize(1)` starts the recursive descent at resolution 1.
6. At each resolution $v$:
   - Generate random perturbation within the safe mask.
   - XOR perturbation into each bitmask parameter.
   - Call `step(v)` to compute $d(v)$ incrementally.
   - If within tolerance, recurse to $v+1$.
7. After the recursion completes, restore the best parameters found.
8. Run a final `compute_all()` verification to confirm the gaps.


## 7. Complexity

### Per perturbation trial at resolution $v$

| Operation | Cost |
|-----------|------|
| Generate random perturbation | $O(w)$ |
| Update bitmask parameter + push to C++ | $O(w)$ |
| `step(v)`: recompute one generating polynomial | $O(k)$ generator steps |
| `step(v)`: Lenstra reduction at resolution $v$ | $O(v^2 \cdot k / 64)$ word ops |

The dominant cost is the Lenstra reduction.  At low resolutions
($v \ll L$), this is very fast; at high resolutions ($v \approx L$),
it is $O(L^2 \cdot k / 64)$.

### Total for one `_run_once` pass

With `max_essais` trials and $L$ resolutions, the worst case is:

$$O\!\left(\text{max\_essais} \cdot L^2 \cdot \frac{k}{64}\right)$$

In practice, the budget resets on improvement and trials are spread
across resolutions, so the actual work is much less.

### Comparison with the backtracking algorithm

| | Backtracking (old) | Safe-mask random (current) |
|---|---|---|
| Perturbation | Single-bit flip | Random multi-bit XOR |
| Candidate selection | Probe-based ($2 \lvert\mathcal{B}\rvert$ evals) | Analytical safe mask |
| Lattice update | Full rebuild per candidate | Incremental `step(v)` |
| Parallelism | Sequential candidates | Independent random trials |
| Large $k$ | Expensive (rebuild cost) | Efficient (incremental) |


## 8. Example: MT19937

For MT19937 with MK type II tempering ($w = 32$, $k = 19937$):

- Parameters: $b$ (32 bits) and $c$ (32 bits), with shifts
  $\eta = 7$, $\mu = 15$, $u = 11$, $l = 18$.
- Safe mask at $v = 1$: 32 + 32 = 64 bits (all free).
- Safe mask at $v = 16$: approximately 17 + 17 = 34 bits.
- Safe mask at $v = 32$: 1 + 1 = 2 bits.
- With `max_essais = 400`: typically reaches ME ($\Delta = 0$) in
  under a minute using the lattice method.

For MELG19937-64 with lagged tempering ($w = 64$, $k = 19937$):

- Parameter: $b$ (64 bits), with structural $\sigma$ and $L$.
- Safe mask at $v = 1$: 64 bits.
- Safe mask at $v = 64$: 1 bit.
- The minimize mode with `n_restarts = 10` typically finds ME by
  trying different random starting points for $b$.
