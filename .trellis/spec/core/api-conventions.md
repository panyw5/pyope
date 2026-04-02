# API Conventions

> How to use and extend the pyope public API correctly.

---

## Core API Functions

### 1. Operator Declaration

```python
from pyope import BasisOperator, Bosonic, Fermionic

# Declare generators with conformal weight
T = BasisOperator("T", conformal_weight=2)
W = BasisOperator("W", conformal_weight=3)
J = BasisOperator("J", conformal_weight=1)

# Register parity (REQUIRED before any OPE computation)
Bosonic(T, W, J)           # bosonic operators (parity=0)
Fermionic(b, c)            # fermionic operators (parity=1)
```

### 2. OPE Definition via MakeOPE

```python
from pyope import OPE, MakeOPE, One, Zero, d

# Pole list: [highest_pole, ..., simple_pole]
# Virasoro: T(z)T(w) ~ c/2·𝟙/(z-w)⁴ + 2T/(z-w)² + ∂T/(z-w)
OPE[T, T] = MakeOPE([
    c/2 * One,   # (z-w)^{-4} pole
    Zero,         # (z-w)^{-3} pole (zero)
    2 * T,        # (z-w)^{-2} pole
    d(T),         # (z-w)^{-1} pole
])
```

**Pole ordering**: List goes from **highest** to **lowest** pole order. The list length determines `max_pole`.

**Scalar poles**: Use `One` (not `1`) for scalar (c-number) poles. Use `Zero` (not `0`) for explicit zero poles, though `0` also works in MakeOPE.

### 3. OPE Computation

```python
result = OPE(A, B)          # returns OPEData
result.max_pole              # highest pole order
result.pole(n)               # n-th pole coefficient (operator expression)
```

### 4. Normal-Ordered Product

```python
from pyope import NO

AB = NO(A, B)               # normal-ordered product (AB)(z)
ABC = NO(A, NO(B, C))       # nested: (A(BC))(z)
```

**FORBIDDEN**: Direct operator multiplication `A * B`. This raises `IllegalOperatorProductError`. Always use `NO(A, B)`.

### 5. Bracket (Mode Extraction)

```python
from pyope import bracket

# bracket(A, B, n) = n-th pole of OPE(A, B)
# WARNING: this is the n-th POLE, NOT the n-th MODE
pole_2 = bracket(T, T, 2)   # = 2*T
```

### 6. Derivatives

```python
from pyope import d, dn

dT = d(T)          # ∂T  (conformal weight = h+1)
d2T = dn(2, T)     # ∂²T (conformal weight = h+2)
d3T = dn(3, T)     # ∂³T
# Equivalent: d(T, 2) also gives ∂²T
```

### 7. Simplification

```python
from pyope import simplify

expr = simplify(some_operator_expression)
# Canonicalizes: expands Wick contractions, sorts normal-ordered products
# expand_derivatives=True by default (applies Leibniz rule)
```

---

## Conformal Weight Rules

| Expression | Weight |
|-----------|--------|
| `BasisOperator("T", conformal_weight=h)` | $h$ |
| `d(A)` where $\text{wt}(A) = h$ | $h + 1$ |
| `dn(n, A)` | $h + n$ |
| `NO(A, B)` | $\text{wt}(A) + \text{wt}(B)$ |

---

## OPE Definition Patterns

### Virasoro Algebra

```python
# T(z)T(w) ~ c/2 / (z-w)⁴ + 2T / (z-w)² + ∂T / (z-w)
OPE[T, T] = MakeOPE([c/2 * One, Zero, 2*T, d(T)])
```

### Primary Field of Weight h

```python
# T(z)φ(w) ~ h·φ / (z-w)² + ∂φ / (z-w)
OPE[T, phi] = MakeOPE([h * phi, d(phi)])
```

### Kac-Moody Current Algebra

```python
# Jᵃ(z)Jᵇ(w) ~ k·δᵃᵇ / (z-w)² + ifᵃᵇc·Jᶜ / (z-w)
OPE[Ja, Jb] = MakeOPE([k * delta_ab * One, f_abc * Jc])
```

### Fermionic (Ghost) Systems

```python
b = BasisOperator("b", conformal_weight=2)
c_ghost = BasisOperator("c", conformal_weight=-1)
Fermionic(b, c_ghost)

# b(z)c(w) ~ 𝟙/(z-w)
OPE[b, c_ghost] = MakeOPE([One])
```

---

## Common Pitfalls

### 1. Missing parity registration

```python
# BAD: forgot Bosonic() → computation will fail
T = BasisOperator("T", conformal_weight=2)
OPE[T, T] = MakeOPE([...])  # Error or wrong result

# GOOD: always register parity before defining OPE
T = BasisOperator("T", conformal_weight=2)
Bosonic(T)
OPE[T, T] = MakeOPE([...])
```

### 2. Direct operator multiplication

```python
# BAD: raises IllegalOperatorProductError
result = T * W

# GOOD: use NO()
result = NO(T, W)
```

### 3. Confusing bracket with mode

```python
# bracket(A, B, n) is the n-th POLE coefficient, not the n-th mode
# For Virasoro: bracket(T, T, 2) = 2*T (the 2nd-order pole coefficient)
```

### 4. Wrong pole ordering in MakeOPE

```python
# BAD: lowest-to-highest (wrong order)
MakeOPE([d(T), 2*T, Zero, c/2*One])

# GOOD: highest-to-lowest
MakeOPE([c/2*One, Zero, 2*T, d(T)])
```
