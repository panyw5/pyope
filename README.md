# PyOPE

PyOPE is a Python library for symbolic Operator Product Expansion (OPE) calculations in Vertex Operator Algebras (VOA) and 2d Conformal Field Theory (CFT).

It is built on top of Python and SymPy, with the goal of providing a programmable, testable, and extensible environment for OPE computations while staying close to the Mathematica reference implementation `OPEdefs.m`.

## Status

- Version: `0.1.0.post1`
- Development status: **Alpha**
- Python: `>=3.8`
- Main dependencies: `sympy`, `numpy`

The project currently supports:

- Registration and evaluation of OPEs between basic generators
- OPE rules for derivatives, linear combinations, and normal-ordered operators
- Construction and simplification of `NO(...)` and `normal_product(...)`
- Extraction of pole coefficients via `bracket(A, B, n)`
- Jacobi identity checks
- A set of experimental tools for C2 spaces, realizations, and null-state searches

The C2, null-state, and realization-related interfaces are still evolving and may change in later releases.

## Installation

Install from PyPI:

```bash
pip install pyope-voa
```

The distribution name is `pyope-voa`, but the import name remains:

```python
import pyope
```

Install from source:

```bash
pip install -e .
```

Install development dependencies:

```bash
pip install -e ".[dev]"
```

Install notebook-related dependencies:

```bash
pip install -e ".[jupyter]"
```

## Quick Start

Here is a minimal Virasoro example defining the OPE of $T(z)T(w)$:

```python
import sympy as sp

from pyope import BasicOperator, Bosonic, MakeOPE, OPE
from pyope import One, Zero, NO, bracket, d

T = BasicOperator("T", conformal_weight=2)
Bosonic(T)

c = sp.Symbol("c")

OPE[T, T] = MakeOPE(
    [
        sp.Rational(1, 2) * c * One,
        Zero,
        2 * T,
        d(T),
    ]
)

tt = OPE(T, T)

print("max pole =", tt.max_pole)
print("{TT}_4 =", bracket(T, T, 4))
print("{TT}_2 =", bracket(T, T, 2))
print("(TT) =", NO(T, T))
```

`MakeOPE([...])` follows the same convention as the Mathematica package: entries are listed from the highest pole down to the $(z-w)^{-1}$ term.

In the example above, the list corresponds to:

- the $(z-w)^{-4}$ coefficient
- the $(z-w)^{-3}$ coefficient
- the $(z-w)^{-2}$ coefficient
- the $(z-w)^{-1}$ coefficient

## Core Concepts

### 1. Basis operators

In typical usage, you first define generators and then declare their statistics:

```python
from pyope import BasicOperator, Bosonic

J = BasicOperator("J", conformal_weight=1)
Bosonic(J)
```

### 2. Registering OPE data

Define an OPE with:

```python
OPE[A, B] = MakeOPE([...])
```

Evaluate it with:

```python
result = OPE(A, B)
```

The result is an `OPEData` object, and individual poles can be accessed with `.pole(n)`.

### 3. Normal-ordered products

PyOPE intentionally rejects direct multiplication such as `A * B` for local operators, to avoid confusing VOA operator syntax with ordinary multiplication.

Use one of the following instead:

- `NO(A, B)`
- `normal_product(A, B, C, ...)`

Example:

```python
from pyope import BasicOperator, Bosonic, normal_product, simplify

A = BasicOperator("A")
B = BasicOperator("B")
Bosonic(A, B)

expr = normal_product(B, A, B)
print(simplify(expr))
```

### 4. Derivatives and brackets

- `d(A)` denotes the first derivative
- `dn(n, A)` denotes the $n$th derivative
- `bracket(A, B, n)` extracts the $n$th bracket / pole coefficient

```python
from pyope import bracket, d, dn

print(bracket(T, T, 4))
print(bracket(d(T), T, 3))
print(dn(2, T))
```

### 5. Jacobi identity checks

```python
from pyope import verify_jacobi_identity

print(verify_jacobi_identity(T, T, T))
```

## Available API

Frequently used public interfaces include:

- `OPE`, `MakeOPE`, `NO`, `NO_product`, `normal_product`, `bracket`
- `BasicOperator`, `Operator`, `d`, `dn`
- `One`, `Zero`, `Delta`
- `simplify`
- `check_jacobi_identity`, `verify_jacobi_identity`

The package also exports a number of more research-oriented and experimental tools, including:

- `C2Space`, `C2NullSearcher`, `GenericC2Reducer`
- `DescendantSpace`
- `SingularVectorAnalyzer`
- `RealizationBackend` and related realization helpers

For the full export list, see `src/pyope/__init__.py`.

## Examples And References

The repository already contains a number of examples and reference materials:

- GitHub repository: <https://github.com/panyw5/pyope>
- Demos and notebooks: <https://github.com/panyw5/pyope/tree/main/demo>
- Test framework notes: <https://github.com/panyw5/pyope/blob/main/tests/TEST_FRAMEWORK.md>
- Mathematica reference materials: <https://github.com/panyw5/pyope/tree/main/OPEdefs>

Current examples cover:

- Virasoro
- Kac-Moody
- Jacobi identities
- Several W-algebra and null-state experiments

## Running Tests

Run the full test suite:

```bash
python -m pytest
```

Run only Mathematica-reference tests:

```bash
python -m pytest -m mathematica_ref
```

Skip slow tests:

```bash
python -m pytest -m "not slow"
```

## Packaging Notes

The current PyPI release is centered on the core package under `src/pyope`.

- the `wheel` contains the core library and package metadata
- notebooks, `.wls` files, and temporary research scripts are not included in the current wheel

That means users installing with `pip install pyope-voa` get the core library rather than the full research repository.

## Citation And Background

- K. Thielemans, "An Algorithmic Approach to Operator Product Expansions, W-algebras and W-strings", arXiv:hep-th/9506159

## License

MIT
