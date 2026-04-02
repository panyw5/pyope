# Architecture

> Module organization, type hierarchy, and dependency flow for pyope.

---

## Layer Architecture

pyope has three layers. Dependencies flow **downward only**.

```
┌─────────────────────────────────────────────────────────┐
│  Research Layer (operator_spaces, c2, null_search, ...)  │
│  → Builds on Core to study algebraic structures          │
├─────────────────────────────────────────────────────────┤
│  Core Computation Layer (api, simplify, jacobi)          │
│  → OPE computation, normal ordering, verification        │
├─────────────────────────────────────────────────────────┤
│  Foundation Layer (operators, registry, ope_data, ...)   │
│  → Type system, global state, data structures            │
└─────────────────────────────────────────────────────────┘
```

### Foundation Layer

| Module | Responsibility | Key exports |
|--------|---------------|-------------|
| `operators.py` | Operator type hierarchy (SymPy Symbol subclasses) | `Operator`, `BasisOperator`, `DerivativeOperator`, `NormalOrderedOperator`, `d()`, `dn()` |
| `registry.py` | Global OPE/parity registration | `OPERegistry`, `ope_registry`, `Bosonic()`, `Fermionic()` |
| `ope_data.py` | OPE result container (pole list) | `OPEData`, `OPEData.from_list()`, `.pole(n)`, `.max_pole` |
| `constants.py` | Identity/zero/delta constants | `One`, `Zero`, `Delta`, `ConstantOperator` |
| `local_operator.py` | Local operator protocol + algebra | `LocalOperator`, `OperatorSum`, `OperatorProduct` |
| `exceptions.py` | Custom error types | `IllegalOperatorProductError` |
| `cache.py` | Computation caching | `OPECache`, `get_ope_cache()` |

### Core Computation Layer

| Module | Responsibility | Key exports |
|--------|---------------|-------------|
| `api.py` | All OPE computation | `OPE()`, `NO()`, `bracket()`, `MakeOPE()`, `normal_product()` |
| `simplify.py` | Expression canonicalization | `simplify()`, Wick contraction, nested NO expansion |
| `jacobi.py` | Algebra consistency checks | `check_jacobi_identity()`, `verify_jacobi_identity()` |
| `quasiprimary.py` | Quasiprimary product formula | `qp()`, `quasiprimary_product()` |
| `compact_ope.py` | Compact sl(2)-covariant notation | `compact_family_poles()` |

### Research Layer

| Module | Responsibility | Key exports |
|--------|---------------|-------------|
| `operator_spaces.py` | Operator space enumeration and linear algebra | `LocalOperatorBasis`, `C2Space`, `RealizedGenerator`, `realize_and_simplify`, `list_zero_relations` |
| `c2.py` | Abstract C2 quotient reduction | `GenericC2Reducer`, `AbstractC2Reducer`, `C2ReductionWitness` |
| `null_search.py` | Null state search orchestration | `C2NullSearcher`, `NullSearchResult` |
| `descendants.py` | Descendant state generation | `DescendantSpace` |
| `singularity.py` | Singular vector analysis | `SingularVectorAnalyzer` |
| `realizations.py` | Realization backend plugins | `RealizationBackend`, `IdentityRealizationBackend` |
| `free_field_c2.py` | Free-field C2 reduction | `FreeFieldC2Reducer` |

---

## Type Hierarchy

```
sympy.Symbol
└── Operator                    (base: has .parity, .conformal_weight)
    ├── BasisOperator           (user-defined generator, e.g. T, W, J)
    ├── DerivativeOperator      (∂ⁿA, created by d(A) / dn(n, A))
    └── NormalOrderedOperator   (NO(A, B), created by NO())
```

**Critical rule**: All operators are SymPy Symbols. This means:
- They participate in SymPy expression arithmetic (`Add`, `Mul`, `Pow`)
- Scalar coefficients are SymPy `Number`/`Symbol` instances
- `expr.has(Operator)` checks if an expression contains operators

---

## Global State

The library uses a **global registry** pattern:

```python
from pyope.registry import ope_registry  # Singleton

# User defines algebra:
Bosonic(T, W)                    # registers parity
OPE[T, T] = MakeOPE([...])      # registers OPE data

# Registry is consulted during computation:
OPE(T, T)  # looks up ope_registry internally
```

**Test isolation**: Every test must clear the registry. The `conftest.py` `disable_cache_for_tests` fixture handles this automatically via `ope_registry.clear()`.

---

## Dependency Rules

1. **Foundation → nothing** (no upward imports)
2. **Core → Foundation only**
3. **Research → Core + Foundation**
4. **No circular imports** — use lazy imports (`from .xxx import` inside functions) if needed
5. **`__init__.py` re-exports everything** — users import from `pyope` directly

---

## Adding a New Module

1. Determine which layer it belongs to
2. Create `src/pyope/<module>.py`
3. If it has a public API, create `src/pyope/<module>.pyi` type stub
4. Add exports to `src/pyope/__init__.py`
5. Add tests to `tests/test_<module>.py`
6. Follow the patterns of the closest existing module in the same layer
