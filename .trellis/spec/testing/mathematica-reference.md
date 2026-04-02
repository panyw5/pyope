# Mathematica Reference Testing

> How to generate, store, and use OPEdefs.m reference data for pyope tests.

---

## The Golden Rule

> **Every assert on an OPE/bracket/NO computation result MUST have a traceable reference.**
>
> Acceptable sources (in priority order):
> 1. OPEdefs.m direct computation output
> 2. Published paper formula with explicit citation (author, year, equation number)
> 3. Mathematical identity derivable from VOA axioms (with proof sketch)

**NEVER**: Guess, hand-calculate, or "it should be this" without verification.

---

## OPEdefs.m Reference Workflow

### Step 1: Write the Mathematica Script

Create a temporary `.wls` (Wolfram Language Script) file:

```mathematica
(* tmp_verify_<test_name>.wls *)
(* Reference computation for test_<test_name>.py *)

SetDirectory[NotebookDirectory[]];
Get["OPEdefs/OPEdefs.m"];

(* Define the algebra *)
Boson[T, 2];
OPE[T, T] = MakeOPE[c/2, 0, 2 T, T'];

(* Compute and print *)
result = OPE[T, \[PartialD]T];
Print["max_pole = ", Length[result]];
Do[
  Print["pole(", i, ") = ", result[[i]]],
  {i, Length[result]}
];
```

### Step 2: Run in Mathematica

```bash
wolframscript -file tmp_verify_<test_name>.wls
```

### Step 3: Record Results in Test

```python
def test_t_derivative_ope(self):
    """
    OPE[T, ∂T] verification.

    Reference: OPEdefs.m direct computation
    Mathematica command:
        Boson[T, 2];
        OPE[T, T] = MakeOPE[c/2, 0, 2 T, T'];
        OPE[T, T'] // Poles
    Result:
        pole(5) = 2*c
        pole(4) = 0
        pole(3) = 4*T
        pole(2) = 3*T'
        pole(1) = T''
    """
    # ... setup ...
    result = OPE(T, d(T))
    assert simplify(result.pole(5) - 2*c*One) == 0
    assert result.pole(4) == 0
    assert simplify(result.pole(3) - 4*T) == 0
    assert simplify(result.pole(2) - 3*d(T)) == 0
    assert simplify(result.pole(1) - dn(2, T)) == 0
```

---

## OPEdefs.m ↔ pyope Function Mapping

| OPEdefs.m (Mathematica) | pyope (Python) | Notes |
|-------------------------|---------------|-------|
| `Boson[T, 2]` | `T = BasisOperator("T", conformal_weight=2); Bosonic(T)` | |
| `Fermion[b, 2]` | `b = BasisOperator("b", conformal_weight=2); Fermionic(b)` | |
| `OPE[A, B] = MakeOPE[...]` | `OPE[A, B] = MakeOPE([...])` | Python uses list `[]` |
| `OPE[A, B]` (compute) | `OPE(A, B)` | `__call__` vs `__getitem__` |
| `NO[A, B]` | `NO(A, B)` | |
| `OPE[A, B][[n]]` | `OPE(A, B).pole(n)` | 1-indexed in both |
| `A'` / `D[A]` | `d(A)` | Derivative |
| `Dn[n, A]` | `dn(n, A)` | n-th derivative |
| `SimplifyNO[expr]` | `simplify(expr)` | Normal ordering simplification |

### Key Differences

1. **MakeOPE pole order**: Both use highest-to-lowest. But OPEdefs.m uses `MakeOPE[a, b, c]` while pyope uses `MakeOPE([a, b, c])` (list).

2. **Scalar poles**: OPEdefs.m uses plain numbers (`c/2`). pyope uses `c/2 * One` for c-number poles.

3. **Zero poles**: OPEdefs.m uses `0`. pyope accepts both `0` and `Zero`.

---

## Acceptable Reference Sources

### Tier 1: OPEdefs.m Direct Computation (Preferred)

The gold standard. Run the equivalent computation in OPEdefs.m and compare output.

### Tier 2: Published Literature

Acceptable with proper citation:

```python
def test_virasoro_ope(self):
    """
    Reference: Di Francesco, Mathieu, Sénéchal,
    "Conformal Field Theory" (1997), eq. (6.23)
    """
```

### Tier 3: VOA Axiom Derivation

For structural properties (not numerical values):

```python
def test_conformal_weight_additivity(self):
    """
    Follows from VOA axiom: wt(NO(A,B)) = wt(A) + wt(B)
    No Mathematica reference needed for structural properties.
    """
```

---

## What Does NOT Need Mathematica Reference

- Operator declaration tests (weight, parity, name)
- Type system tests (isinstance checks, error raising)
- Structural properties (weight additivity, derivative weight)
- API behavior tests (MakeOPE from list, pole access)
- Regression tests for Python-specific bugs

---

## Storing Reference Data

For complex algebras with many reference values, create a dedicated file:

```
tests/
├── test_w3_algebra_ref.py         # Uses reference data
└── parse_mathematica_results.py   # Helper to parse Mathematica output
```

Large reference datasets can be stored as Python dicts in the test file:

```python
# Reference data from OPEdefs.m computation (2024-03-19)
# Script: tmp_verify_w3_ww_ope.wls
W3_WW_REFERENCE = {
    "max_pole": 6,
    "poles": {
        6: "c",
        5: "0",
        4: "2*T",
        3: "T'",
        2: "2*beta*Lambda + 3/10*T''",
        1: "beta*Lambda' + 1/15*T'''",
    }
}
```
