# Influence Analysis for $\mathbb{F}_2$-Linear Functions

> **Notation.** This page follows the conventions in
> [Notation](notation.md): bit vectors are bold ($\mathbf{x}$,
> $\mathbf{y}$, $\mathbf{b}_k$, $\mathbf{e}_j$), matrices are plain
> capitals ($A$, $B_k$), scalars are plain italic ($k$, $L$, $j$, $v$).
> **Page-local override:** on this page $v$ denotes an **output bit
> index** $v \in \{0, \ldots, L-1\}$ (the $v$-th component of
> $\mathbf{y}$); this is distinct from the resolution $\ell$ used
> elsewhere. The letter $n$ here counts **parameter vectors**, not
> time steps. The state size $k$ from the notation page is reused
> unchanged, but the outer summation index $k = 1, \ldots, n$ is the
> parameter-vector index.

## Problem Statement

Let $f : \mathbb{F}_2^k \times \mathbb{F}_2^{d_1} \times \cdots \times \mathbb{F}_2^{d_n} \to \mathbb{F}_2^L$ be a function that is **linear over $\mathbb{F}_2$** in all its arguments jointly, where $k$ is the state size, $L$ the output bit width, and $d_k$ the dimension of the $k$-th parameter vector $\mathbf{b}_k \in \mathbb{F}_2^{d_k}$. Given a target output bit index $v$, determine which bits of each parameter vector $\mathbf{b}_1, \ldots, \mathbf{b}_n$ affect $y_v$, where $\mathbf{y} = f(\mathbf{x}, \mathbf{b}_1, \ldots, \mathbf{b}_n)$ and $\mathbf{x} \in \mathbb{F}_2^k$ is the state.

## Decomposition

Since $f$ is $\mathbb{F}_2$-linear, it can be written as:

$$\mathbf{y} = A\, \mathbf{x} \oplus B_1\, \mathbf{b}_1 \oplus \cdots \oplus B_n\, \mathbf{b}_n$$

where $A \in \mathbb{F}_2^{L \times k}$ and $B_k \in \mathbb{F}_2^{L \times d_k}$ are fixed binary matrices. The output bit $y_v$ is determined by row $v$ of each matrix:

$$y_v = A[v,:] \cdot \mathbf{x} \;\oplus\; B_1[v,:] \cdot \mathbf{b}_1 \;\oplus\; \cdots \;\oplus\; B_n[v,:] \cdot \mathbf{b}_n$$

The influence of $\mathbf{b}_k$ on $y_v$ is exactly the support of row $v$ of $B_k$:

$$b_{k,j} \text{ influences } y_v \iff B_k[v, j] = 1$$

where $b_{k,j}$ denotes the $j$-th (scalar) component of $\mathbf{b}_k$.

## Key Identity

By linearity, zeroing all inputs except a single basis vector $\mathbf{e}_j \in \mathbb{F}_2^{d_k}$ in position $\mathbf{b}_k$ isolates column $j$ of $B_k$:

$$f(\mathbf{0},\, \ldots,\, \mathbf{0},\, \mathbf{e}_j,\, \mathbf{0},\, \ldots,\, \mathbf{0})_v = B_k[v,\, j]$$

This is the only property needed — no explicit matrix construction is required.

## Algorithm

```
Input:  f    — the F_2-linear function
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

Linearity of $f$ gives $f(\mathbf{0}, \ldots, \mathbf{0}) = \mathbf{0}$, and by the superposition principle each probe $f(\mathbf{0},\ldots,\mathbf{e}_j,\ldots,\mathbf{0})$ returns exactly column $j$ of $B_k$. No assumption is needed about the internal structure of $f$ — only that the whole function is $\mathbb{F}_2$-linear.

## Extension: Composed Structure

If $f$ is built as a composition, e.g.:

$$\mathbf{y} = f_n(\cdots f_2(f_1(\mathbf{x},\, \mathbf{b}_1),\, \mathbf{b}_2) \cdots,\, \mathbf{b}_n)$$

the same algorithm applies without modification, provided the **composed function** is $\mathbb{F}_2$-linear overall. In that case the effective injection matrix for $\mathbf{b}_k$ is:

$$\tilde{B}_k = A_{f_n} \cdots A_{f_{k+1}} \cdot B_{f_k}$$

where $A_{f_j}$ is the state-transition matrix of layer $j$ and $B_{f_k}$ is the direct injection matrix of $\mathbf{b}_k$ at layer $k$. Parameters injected early are filtered through more layers; parameters injected late have more direct influence. The probing algorithm computes $\tilde{B}_k[v,:]$ implicitly without ever forming these products.

## Example: MELG Generator Output

The MELG64 output transformation (Harase & Kimoto 2018) is:

$$\mathbf{y} = \underbrace{\bigl(\mathbf{w}_i \oplus (\mathbf{w}_i \ll \sigma_3)\bigr)}_{T_1 \cdot \mathbf{w}_i(\mathbf{x})} \oplus \underbrace{(\mathbf{w}_{i+L} \;\&\; \mathbf{b})}_{B(\mathbf{x})\, \mathbf{b}}$$

where $\mathbf{w}_i, \mathbf{w}_{i+L} \in \mathbb{F}_2^{64}$ are 64-bit words extracted from the generator state $\mathbf{x}$, $\sigma_3$ is a fixed shift amount (a scalar integer), and $\mathbf{b} \in \mathbb{F}_2^{64}$ is the tunable mask parameter. This is **linear in $\mathbf{b}$ for fixed $\mathbf{x}$**, with injection matrix $B(\mathbf{x}) = \operatorname{diag}(\mathbf{w}_{i+L}(\mathbf{x}))$. Consequently:

$$b_j \text{ influences } y_v \iff v = j \text{ and } w_{i+L,\,v}(\mathbf{x}) = 1$$

where $w_{i+L,v}(\mathbf{x})$ is the $v$-th (scalar) bit of the word $\mathbf{w}_{i+L}(\mathbf{x})$. Each bit of $\mathbf{b}$ can only affect the output bit at the same position, and only when the corresponding bit of the lookahead word $\mathbf{w}_{i+L}$ is 1. This diagonal structure is what enables the bit-by-bit backtracking search used to achieve $\Delta = 0$.

## See also

- [Notation](notation.md) — canonical symbol table for the whole
  theory section.