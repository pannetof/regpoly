# Tempering Optimization Algorithm

Optimization of tempering bitmask parameters to achieve maximally
equidistributed (ME) generators. Based on Harase (2009), generalized
to any transformation with bitmask parameters.

## Definitions

Let $G$ be a full-period $\mathbb{F}_2$-linear generator with state
size $k$ and word size $w$. Let $T$ be a tempering transformation with
bitmask parameters $b_1, b_2, \ldots, b_m \in \mathbb{F}_2^w$.

The dimension of equidistribution with $v$-bit accuracy is $k(v)$, with
upper bound $\lfloor k/v \rfloor$. The dimension gap at resolution $v$ is

$$d(v) = \lfloor k/v \rfloor - k(v), \qquad v = 1, \ldots, w.$$

The total dimension defect is

$$\Delta = \sum_{v=1}^{w} d(v).$$

The goal is to find $b_1, \ldots, b_m$ such that $\Delta = 0$.

## Notation

- $b_i[j]$ denotes bit $j$ of bitmask parameter $b_i$, for $j = 0, \ldots, w-1$.
- $b_i \oplus e_j$ denotes $b_i$ with bit $j$ flipped.
- $o_v(s)$ denotes the $v$ most significant bits of the output of $T$ applied to state $s$.
- $\mathcal{B} = \{(i, j) : i = 1, \ldots, m, \; j = 0, \ldots, w-1\}$ is the set of all flippable bits.

## Phase 1: Candidate Filtering

For each resolution $v$, compute the subset of bits that can affect
the $v$ most significant output bits.

Choose a random test state $s \in \mathbb{F}_2^k$, $s \neq 0$.
Compute the reference output $y = o_v(s)$.

For each $(i, j) \in \mathcal{B}$:

1. Set $b_i \leftarrow b_i \oplus e_j$.
2. Compute $y' = o_v(s)$.
3. If $y' \neq y$, then $(i, j)$ is a candidate for resolution $v$.
4. Revert: $b_i \leftarrow b_i \oplus e_j$.

This defines the candidate set

$$\mathcal{C}(v) = \{(i, j) \in \mathcal{B} : o_v(s; \, b_i \oplus e_j) \neq o_v(s; \, b_i)\}.$$

**Cost:** $2 \times |\mathcal{B}|$ applications of $T$ per resolution, negligible compared to the matrix computations in Phase 2.

## Phase 2: Backtracking Optimization

### Data structures

- Stack $S$: stores triples $(v, \, c, \, \hat{b})$ where $v$ is the resolution, $c \in \mathcal{C}(v)$ is the candidate that was flipped, and $\hat{b}$ is the saved parameter state before the flip.
- Counter $n_{\text{bt}}$: number of backtracks performed.
- Limit $N_{\text{bt}}$: maximum allowed backtracks.

### Algorithm

$$\text{Set } v \leftarrow 1, \quad c_{\text{start}} \leftarrow 0, \quad S \leftarrow \emptyset, \quad n_{\text{bt}} \leftarrow 0.$$

**While** $v \leq w$:

$\quad$ 1. Build the generator matrix and compute $d(v)$.

$\quad$ 2. **If** $d(v) = 0$: set $v \leftarrow v + 1$, $c_{\text{start}} \leftarrow 0$, **continue**.

$\quad$ 3. **For each** $(i, j) \in \mathcal{C}(v)$ with index $\geq c_{\text{start}}$:

$\qquad$ a. Save parameters: $\hat{b} \leftarrow (b_1, \ldots, b_m)$.

$\qquad$ b. Flip: $b_i \leftarrow b_i \oplus e_j$.

$\qquad$ c. Recompute $d(v')$ for all $v' = 1, \ldots, v$.

$\qquad$ d. **If** $d(v') = 0$ for all $v' \leq v$:

$\qquad\quad$ Push $(v, \text{index of } (i,j), \hat{b})$ onto $S$.

$\qquad\quad$ Set $v \leftarrow 1$, $c_{\text{start}} \leftarrow 0$.

$\qquad\quad$ **Go to** step 1.

$\qquad$ e. **Else:** revert $b_i \leftarrow \hat{b}_i$.

$\quad$ 4. No candidate fixed resolution $v$. **Backtrack:**

$\qquad$ **If** $S = \emptyset$ or $n_{\text{bt}} \geq N_{\text{bt}}$: **stop** (cannot reach $\Delta = 0$).

$\qquad$ Pop $(v_0, c_0, \hat{b})$ from $S$.

$\qquad$ Restore $(b_1, \ldots, b_m) \leftarrow \hat{b}$.

$\qquad$ Set $v \leftarrow v_0$, $c_{\text{start}} \leftarrow c_0 + 1$, $n_{\text{bt}} \leftarrow n_{\text{bt}} + 1$.

$\qquad$ **Go to** step 1.

### Termination

The algorithm terminates when either:

- $v > w$: all gaps are zero, $\Delta = 0$ (ME achieved).
- The stack is empty and no candidate fixes the current resolution.
- The backtrack limit $N_{\text{bt}}$ is reached.

## Complexity

At each resolution $v$, the algorithm tries at most $|\mathcal{C}(v)|$
candidates. Each candidate requires:

- Rebuilding the generator matrix via Gaussian elimination: $O(k^2 \cdot L)$.
- Computing $d(v')$ for $v' = 1, \ldots, v$: $O(k^2 \cdot v)$.

For $w = 64$ with MK tempering, $|\mathcal{C}(v)|$ is typically $20$-$40$.
With the candidate filtering, the total work is roughly

$$\sum_{v=1}^{w} |\mathcal{C}(v)| \cdot O(k^2 v)$$

which is much smaller than the naive $|\mathcal{B}| \cdot w \cdot O(k^2 w)$ without filtering.
