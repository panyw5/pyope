# OPEdefs.m → pyope Mapping

> Complete reference for translating between Thielemans' OPEdefs.m (Mathematica) and pyope (Python).

---

## Operator Declaration

| OPEdefs.m | pyope | Notes |
|-----------|-------|-------|
| `Boson[T, 2]` | `T = BasisOperator("T", conformal_weight=2); Bosonic(T)` | Two-step in Python |
| `Fermion[b, 2]` | `b = BasisOperator("b", conformal_weight=2); Fermionic(b)` | |
| `Boson[T, W, 2, 3]` | `T = BasisOperator("T", conformal_weight=2); W = BasisOperator("W", conformal_weight=3); Bosonic(T, W)` | OPEdefs allows batch |

---

## OPE Definition and Computation

| OPEdefs.m | pyope | Notes |
|-----------|-------|-------|
| `OPE[A, B] = MakeOPE[p1, p2, ...]` | `OPE[A, B] = MakeOPE([p1, p2, ...])` | Python needs `[]` list |
| `OPE[A, B]` (compute) | `OPE(A, B)` | `()` = compute, `[]` = define |
| `OPE[A, B][[n]]` | `OPE(A, B).pole(n)` | Both 1-indexed |
| `Poles[OPE[A, B]]` | `OPE(A, B).max_pole` | Highest pole order |

### Pole Ordering

Both systems use **highest-pole-first** ordering:

```
MakeOPE[c/2, 0, 2T, T']     (* Mathematica: 4th, 3rd, 2nd, 1st pole *)
MakeOPE([c/2*One, Zero, 2*T, d(T)])   # Python: same order
```

### Scalar Poles

```
(* Mathematica: plain numbers for c-number poles *)
MakeOPE[c/2, 0, 2T, T']

# Python: use One for c-number poles
MakeOPE([c/2 * One, Zero, 2*T, d(T)])
```

---

## Normal Ordering

| OPEdefs.m | pyope | Notes |
|-----------|-------|-------|
| `NO[A, B]` | `NO(A, B)` | |
| `NO[A, B, C]` | `NO(A, NO(B, C))` | pyope is binary only, nest for multi |
| `NO[A, B] // SimplifyNO` | `simplify(NO(A, B))` | |
| `SimplifyNO[expr]` | `simplify(expr)` | or `expr.simplifyNO()` |

---

## Derivatives

| OPEdefs.m | pyope | Notes |
|-----------|-------|-------|
| `A'` | `d(A)` | First derivative |
| `D[A]` | `d(A)` | Alternative Mathematica syntax |
| `Dn[2, A]` | `dn(2, A)` or `d(A, 2)` | n-th derivative |
| `A''` | `dn(2, A)` | |

---

## Bracket / Mode Extraction

| OPEdefs.m | pyope | Notes |
|-----------|-------|-------|
| `OPE[A, B][[n]]` | `bracket(A, B, n)` | n-th **pole**, not n-th mode |

**Important**: In both systems, `bracket(A, B, n)` or `OPE[A,B][[n]]` returns the **n-th pole coefficient**, which is the coefficient of $(z-w)^{-n}$. This is related to but **not identical** to the n-th mode $A_n$ acting on $B$.

---

## Simplification

| OPEdefs.m | pyope | Notes |
|-----------|-------|-------|
| `SimplifyNO[expr]` | `simplify(expr)` | Main simplifier |
| `ExpandNO[expr]` | `simplify(expr, expand_derivatives=True)` | With Leibniz rule |
| (no equivalent) | `simplify(expr, preserve_nested_structure=True)` | Keep nested NO |

---

## Advanced: C2 and Operator Spaces

| OPEdefs.m concept | pyope equivalent | Notes |
|-------------------|-----------------|-------|
| Manual basis enumeration | `LocalOperatorBasis(generators, weight)` | Automatic enumeration |
| Manual C2 quotient | `C2Space(generators, weight)` | Automatic |
| Manual realization check | `realize_and_simplify(expr)` | Via free field backend |
| (no equivalent) | `GenericC2Reducer` | Abstract C2 check without realization |
| (no equivalent) | `C2NullSearcher` | Automated null state search |

---

## Translation Examples

### Example 1: Virasoro Algebra

**OPEdefs.m:**
```mathematica
Boson[T, 2];
OPE[T, T] = MakeOPE[c/2, 0, 2 T, T'];
result = OPE[T, T'] // SimplifyNO;
Print[result[[1]]]; (* highest pole *)
```

**pyope:**
```python
T = BasisOperator("T", conformal_weight=2)
Bosonic(T)
c = sp.Symbol("c")
OPE[T, T] = MakeOPE([c/2 * One, Zero, 2*T, d(T)])
result = OPE(T, d(T))
print(result.pole(result.max_pole))  # highest pole
```

### Example 2: W₃ with Composite Operator

**OPEdefs.m:**
```mathematica
Lambda = NO[T, T] - 3/10 Dn[2, T];
OPE[W, W] = MakeOPE[c, 0, 2T, T', 2 beta Lambda + 3/10 Dn[2,T],
                     beta D[Lambda] + 1/15 Dn[3,T]];
```

**pyope:**
```python
Lambda = NO(T, T) - sp.Rational(3, 10) * dn(2, T)
OPE[W, W] = MakeOPE([
    c * One,
    Zero,
    2 * T,
    d(T),
    2 * beta * Lambda + sp.Rational(3, 10) * dn(2, T),
    beta * d(Lambda) + sp.Rational(1, 15) * dn(3, T),
])
```

### Example 3: Ghost System

**OPEdefs.m:**
```mathematica
Fermion[b, 2];
Fermion[c, -1];
OPE[b, c] = MakeOPE[1];
```

**pyope:**
```python
b = BasisOperator("b", conformal_weight=2)
c_ghost = BasisOperator("c", conformal_weight=-1)
Fermionic(b, c_ghost)
OPE[b, c_ghost] = MakeOPE([One])
```
