# Testing Guidelines

> How to write and organize tests for pyope, with emphasis on Mathematica reference validation.

---

## Overview

pyope testing has a **unique constraint**: computational results must be validated against the reference implementation (**Thielemans' OPEdefs.m**) or published CFT/VOA literature. Hand-written expected values for OPE computations are **forbidden**.

---

## Guidelines Index

| Guide | Description | Status |
|-------|-------------|--------|
| [Test Conventions](./test-conventions.md) | Test structure, fixtures, naming, markers | Filled |
| [Mathematica Reference](./mathematica-reference.md) | How to generate and use OPEdefs.m reference data | Filled |
| [Algebra Fixtures](./algebra-fixtures.md) | Standard VOA fixtures and reuse patterns | Filled |

---

## Pre-Development Checklist

Before writing any test:

1. **Always**: [Test Conventions](./test-conventions.md) — understand test structure
2. **If testing OPE/NO/bracket computations**: [Mathematica Reference](./mathematica-reference.md) — mandatory
3. **If defining a VOA algebra**: [Algebra Fixtures](./algebra-fixtures.md) — reuse existing fixtures

---

## Quick Reference: Test Categories

| Category | What it tests | Reference source | Example |
|----------|--------------|-----------------|---------|
| **Definition** | Operator declaration, parity, weight | VOA axioms | `test_stress_tensor_declaration` |
| **Computation** | OPE, NO, bracket results | **OPEdefs.m** (mandatory) | `test_t_derivative_ope` |
| **Identity** | Jacobi, associativity, Ward | Mathematical theorem | `test_jacobi_identity` |
| **Property** | Primary field condition, weight additivity | CFT definition | `test_primary_field_condition` |
| **Regression** | Bug fix, edge case | Previous bug report | `test_q0_case_fix` |
| **Research** | C2 space, null states, realization | OPEdefs.m + papers | `test_c2_null_searcher` |

---

## The Golden Rule

> **Computation tests MUST cite their reference.**
>
> Every `assert` on an OPE/bracket/NO result must have a traceable source:
> either a Mathematica OPEdefs.m computation or a published formula with citation.
