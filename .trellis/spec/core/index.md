# Core Library Guidelines

> Architecture, API conventions, and coding patterns for pyope.

---

## Overview

pyope is a Python library for computing Operator Product Expansions (OPE) in Vertex Operator Algebras (VOA). It is a port of Thielemans' Mathematica package **OPEdefs.m** to Python/SymPy.

**Key constraint**: All symbolic algebra is built on **SymPy**. Operators are `sympy.Symbol` subclasses. Expressions are `sympy.Expr` trees.

---

## Guidelines Index

| Guide | Description | Status |
|-------|-------------|--------|
| [Architecture](./architecture.md) | Module organization, type hierarchy, dependency flow | Filled |
| [API Conventions](./api-conventions.md) | OPE/NO/bracket/MakeOPE usage patterns | Filled |
| [SymPy Patterns](./sympy-patterns.md) | SymPy-specific coding patterns and pitfalls | Filled |
| [Quality Guidelines](./quality-guidelines.md) | Code standards, forbidden patterns, review checklist | Filled |

---

## Pre-Development Checklist

Before writing any pyope code, read these in order:

1. **Always**: [Architecture](./architecture.md) — understand module boundaries
2. **Always**: [API Conventions](./api-conventions.md) — know the public API surface
3. **If touching SymPy expressions**: [SymPy Patterns](./sympy-patterns.md)
4. **Before PR/review**: [Quality Guidelines](./quality-guidelines.md)

---

## Quick Reference: Module Map

```
src/pyope/
├── operators.py         # Type system: Operator → BasisOperator, DerivativeOperator, NormalOrderedOperator
├── registry.py          # Global state: OPERegistry, Bosonic(), Fermionic()
├── ope_data.py          # Data: OPEData (pole storage)
├── api.py               # Core API: OPE(), NO(), bracket(), MakeOPE()
├── simplify.py          # Simplification: Wick contraction, normal ordering canonicalization
├── constants.py         # Special values: One, Zero, Delta
├── local_operator.py    # LocalOperator protocol, OperatorSum, OperatorProduct
├── cache.py             # Performance: OPECache, math function caches
├── jacobi.py            # Verification: Jacobi identity checks
├── operator_spaces.py   # Research: LocalOperatorBasis, C2Space, realization, zero relations
├── c2.py                # Research: AbstractC2Reducer, GenericC2Reducer (realization-free)
├── null_search.py       # Research: C2NullSearcher (quotient precheck + descendant lift)
├── descendants.py       # Research: DescendantSpace (fixed-weight descendant generation)
├── singularity.py       # Research: SingularVectorAnalyzer
├── quasiprimary.py      # Utility: quasiprimary products, qp()
├── compact_ope.py       # Utility: sl(2)-covariant compact notation
├── realizations.py      # Backend: RealizationBackend plugin system
├── free_field_c2.py     # Backend: Free-field C2 reducers
├── exceptions.py        # Errors: IllegalOperatorProductError
└── utils.py             # Helpers: sympify_coefficient, math utilities
```
