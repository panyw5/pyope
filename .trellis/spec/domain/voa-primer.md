# VOA Primer for pyope Development

> Minimal VOA/CFT concepts needed to understand and develop pyope. Not a textbook — just enough to write correct code and tests.

---

## What is a VOA?

A **Vertex Operator Algebra** (VOA) is an algebraic structure where:
- **States** are local operators $A(z), B(z), \ldots$ (fields on a complex plane)
- **Products** are defined via the **Operator Product Expansion** (OPE)
- The algebra satisfies **associativity** (encoded as the Jacobi identity)

In pyope, a VOA is specified by declaring **generators** and their **OPE rules**.

---

## Key Concepts

### 1. Conformal Weight

Every operator $A$ has a **conformal weight** (or **conformal dimension**) $h_A \in \frac{1}{2}\mathbb{Z}_{\geq 0}$.

- $T$ (stress tensor): $h = 2$
- $J$ (current): $h = 1$
- $\partial^n A$: $h = h_A + n$
- $NO(A, B)$: $h = h_A + h_B$

### 2. Parity (Statistics)

Operators are either **bosonic** (parity 0) or **fermionic** (parity 1).

- Bosonic: commuting ($AB = BA$ up to singular terms)
- Fermionic: anticommuting ($AB = -BA$ up to singular terms)

This affects sign rules in OPE computation and normal ordering.

### 3. OPE (Operator Product Expansion)

The OPE of two operators $A(z)$ and $B(w)$ is:

$$A(z) B(w) \sim \sum_{n=1}^{N} \frac{\{AB\}_n(w)}{(z-w)^n} + \text{regular}$$

- $\{AB\}_n$ is the **n-th pole coefficient** (a local operator at $w$)
- $N$ is the **maximum pole order** (`max_pole`)
- The regular (non-singular) part is the **normal-ordered product** $NO(A,B)$

In pyope: `OPE(A, B).pole(n)` gives $\{AB\}_n$.

### 4. Normal-Ordered Product

$NO(A, B)$ (also written $(AB)$ or $:AB:$) is the **regular part** of the OPE:

$$NO(A, B)(w) = \{AB\}_0(w)$$

Properties:
- $\text{wt}(NO(A,B)) = \text{wt}(A) + \text{wt}(B)$
- NOT commutative in general: $NO(A,B) \neq NO(B,A)$
- NOT associative: $NO(A, NO(B,C)) \neq NO(NO(A,B), C)$ in general

### 5. Derivatives

$\partial A$ is the **conformal derivative** of $A$. It satisfies:
- $\text{wt}(\partial A) = \text{wt}(A) + 1$
- $\{A, \partial B\}_n = (n-1)\{A, B\}_{n-1} + \partial\{A,B\}_n$ (Thielemans eq 3.3.2)
- $\{\partial A, B\}_n = -(n-1)\{A, B\}_{n-1}$ (Thielemans eq 3.3.1)

### 6. Central Charge

The **central charge** $c$ appears in the Virasoro OPE:

$$T(z)T(w) \sim \frac{c/2}{(z-w)^4} + \frac{2T}{(z-w)^2} + \frac{\partial T}{(z-w)}$$

It characterizes the VOA. Common values: $c = 1$ (free boson), $c = 26$ (critical string), $c = 1/2$ (Ising model).

### 7. Primary Fields

An operator $\phi$ of weight $h$ is **primary** if:

$$T(z)\phi(w) \sim \frac{h\phi}{(z-w)^2} + \frac{\partial\phi}{(z-w)}$$

No higher-order poles. The stress tensor $T$ itself is **quasi-primary** (not primary — it has a 4th-order pole in $T$-$T$ OPE).

---

## Advanced Concepts (Research Layer)

### C2 Space and Associated Variety

The **C2 space** $C_2(V)$ is spanned by operators of the form $\{AB\}_n$ with $n \geq 2$. The **quotient** $R(V) = V / C_2(V)$ is a commutative algebra (Zhu's $C_2$ algebra).

The **associated variety** $X_V$ is the spectrum of $R(V)$. Its dimension is a key invariant:
- $\dim X_V = 0$: VOA is "C2-cofinite" (rational, well-behaved)
- $\dim X_V > 0$: VOA has more complex representation theory

In pyope: `C2Space` constructs the C2 quotient, `GenericC2Reducer` checks C2 membership without realization.

### Null States

A **null state** is a nonzero linear combination of operators that:
1. Is in $C_2(V)$ (i.e., zero in the quotient $R(V)$), AND
2. Lifts to an actual relation among descendants

Finding null states constrains the algebra structure. In pyope: `C2NullSearcher` implements this search.

### Free Field Realization

A **realization** expresses abstract generators in terms of free fields (bosons/fermions). This allows:
- Explicit computation of all OPEs
- Verification of abstract algebra relations
- Dimensional analysis of $R(V)$

In pyope: `RealizedGenerator` and `realize_and_simplify()` handle this.

---

## Standard Algebras in pyope

| Algebra | Generators | Central charge | Key reference |
|---------|-----------|---------------|---------------|
| Virasoro | $T$ | $c$ | Di Francesco et al., Ch. 6 |
| $\hat{\mathfrak{sl}}(2)_k$ | $J^+, J^0, J^-$ | $3k/(k+2)$ | Di Francesco et al., Ch. 15 |
| $W_3$ | $T, W$ | $c$ | Zamolodchikov (1985) |
| Free boson | $\partial\phi$ | $1$ | — |
| $bc$ ghosts | $b, c$ | $-2(6\lambda^2 - 6\lambda + 1)$ | Polchinski, Ch. 2 |
| $\beta\gamma$ ghosts | $\beta, \gamma$ | $2(6\lambda^2 - 6\lambda + 1)$ | Polchinski, Ch. 2 |
| $\mathcal{N}=3$ SCFT | $T, J, W, G, \tilde{G}, GW, \bar{W}, G\bar{W}$ | $c$ | Beem-Rastelli |
