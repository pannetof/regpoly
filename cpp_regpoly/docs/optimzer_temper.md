# Influence Analysis for GF(2)-Linear Functions

## Problem Statement

Let $f : \text{GF}(2)^k \times \text{GF}(2)^{d_1} \times \cdots \times \text{GF}(2)^{d_n} \to \text{GF}(2)^L$ be a function that is **linear over GF(2)** in all its arguments jointly. Given a target output bit index $v$, determine which bits of each parameter vector $b_1, \ldots, b_n$ affect $y_v$, where $y = f(x, b_1, \ldots, b_n)$.

## Decomposition

Since $f$ is GF(2)-linear, it can be written as:

$$y = A x \oplus B_1 b_1 \oplus \cdots \oplus B_n b_n$$

where $A \in \text{GF}(2)^{L \times k}$ and $B_k \in \text{GF}(2)^{L \times d_k}$ are fixed binary matrices. The output bit $y_v$ is determined by row $v$ of each matrix:

$$y_v = A[v,:] \cdot x \;\oplus\; B_1[v,:] \cdot b_1 \;\oplus\; \cdots \;\oplus\; B_n[v,:] \cdot b_n$$

The influence of $b_k$ on $y_v$ is exactly the support of row $v$ of $B_k$:

$$b_{k,j} \text{ influences } y_v \iff B_k[v, j] = 1$$

## Key Identity

By linearity, zeroing all inputs except a single basis vector $e_j$ in position $b_k$ isolates column $j$ of $B_k$:

$$f(0,\, \ldots,\, 0,\, e_j,\, 0,\, \ldots,\, 0)_v = B_k[v,\, j]$$

This is the only property needed — no explicit matrix construction is required.

## Algorithm

```
Input:  f    — the GF(2)-linear function
        v    — target output bit index
        n    — number of parameter vectors
        d_k  — dimension of b_k, for k = 1..n

Output: influence[k] — set of bit indices j such that b_{k,j} affects y_v

for k = 1..n:
    influence[k] = {}
    for j = 0..d_k - 1:
        if f(0, 0, ..., e_j, ..., 0)[v] == 1:   // e_j in position b_k, all else 0
            influence[k].add(j)

return influence
```

**Cost:** $\displaystyle\sum_{k=1}^{n} d_k$ evaluations of $f$.

## Correctness

Linearity of $f$ gives $f(0, \ldots, 0) = 0$, and by the superposition principle each probe $f(0,\ldots,e_j,\ldots,0)$ returns exactly column $j$ of $B_k$. No assumption is needed about the internal structure of $f$ — only that the whole function is GF(2)-linear.

## Extension: Composed Structure

If $f$ is built as a composition, e.g.:

$$y = f_n(\cdots f_2(f_1(x,\, b_1),\, b_2) \cdots,\, b_n)$$

the same algorithm applies without modification, provided the **composed function** is GF(2)-linear overall. In that case the effective injection matrix for $b_k$ is:

$$\tilde{B}_k = A_{f_n} \cdots A_{f_{k+1}} \cdot B_{f_k}$$

where $A_{f_j}$ is the state-transition matrix of layer $j$ and $B_{f_k}$ is the direct injection matrix of $b_k$ at layer $k$. Parameters injected early are filtered through more layers; parameters injected late have more direct influence. The probing algorithm computes $\tilde{B}_k[v,:]$ implicitly without ever forming these products.

## Example: MELG Generator Output

The MELG64 output transformation (Harase & Kimoto 2018) is:

$$y = \underbrace{\bigl(w_i \oplus (w_i \ll \sigma_3)\bigr)}_{T_1 \cdot w_i(x)} \oplus \underbrace{(w_{i+L} \;\&\; b)}_{B(x)\, b}$$

where $w_i, w_{i+L}$ are words extracted from the generator state $x$, and $b \in \text{GF}(2)^{64}$ is the tunable mask parameter. This is **linear in $b$ for fixed $x$**, with injection matrix $B(x) = \text{diag}(w_{i+L}(x))$. Consequently:

$$b_j \text{ influences } y_v \iff v = j \text{ and } w_{i+L,\,v}(x) = 1$$

Each bit of $b$ can only affect the output bit at the same position, and only when the corresponding bit of the lookahead word $w_{i+L}$ is 1. This diagonal structure is what enables the bit-by-bit backtracking search used to achieve $\Delta = 0$.