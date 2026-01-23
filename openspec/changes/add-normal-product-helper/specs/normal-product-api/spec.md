# Spec: normal_product() API

## ADDED Requirements

### Requirement: Support variable-argument normal product calculation

The system MUST provide a public API function `normal_product()` that computes the nested normal ordered product of multiple operators. The function MUST accept zero or more operator arguments and return their right-associative normal ordered product.

#### Scenario: Compute normal product of zero operators

**Given** the pyope library is imported
**When** the user calls `normal_product()` with no arguments
**Then** the function returns `One` (the identity operator)

**Example**:
```python
from pyope import normal_product, One

result = normal_product()
assert result == One
```

#### Scenario: Compute normal product of single operator

**Given** the pyope library is imported
**And** a basis operator `T` exists
**When** the user calls `normal_product(T)` with one argument
**Then** the function returns the operator itself unchanged

**Example**:
```python
from pyope import BasisOperator, normal_product

T = BasisOperator("T", 2)
result = normal_product(T)
assert result == T
```

#### Scenario: Compute normal product of two operators

**Given** the pyope library is imported
**And** two operators `A` and `B` exist
**When** the user calls `normal_product(A, B)` with two arguments
**Then** the function returns `NO(A, B)` (equivalent to the two-argument NO function)

**Example**:
```python
from pyope import BasisOperator, normal_product, NO

A = BasisOperator("A", 1)
B = BasisOperator("B", 1)
result = normal_product(A, B)
expected = NO(A, B)
assert result == expected
```

#### Scenario: Compute normal product of three operators

**Given** the pyope library is imported
**And** three operators `A`, `B`, `C` exist
**When** the user calls `normal_product(A, B, C)` with three arguments
**Then** the function returns `NO(A, NO(B, C))` (right-associative nesting)

**Example**:
```python
from pyope import BasisOperator, normal_product, NO

T = BasisOperator("T", 2)
J = BasisOperator("J", 1)
W = BasisOperator("W", 3)
result = normal_product(T, J, W)
expected = NO(T, NO(J, W))
assert result == expected
```

#### Scenario: Compute normal product of four or more operators

**Given** the pyope library is imported
**And** multiple operators `A, B, C, D, ...` exist
**When** the user calls `normal_product(A, B, C, D, ...)` with four or more arguments
**Then** the function returns the right-associative nesting `NO(A, NO(B, NO(C, D)))`

**Example**:
```python
from pyope import BasisOperator, normal_product, NO

T = BasisOperator("T", 2)
J = BasisOperator("J", 1)
W = BasisOperator("W", 3)
L = BasisOperator("L", 1)
result = normal_product(T, J, W, L)
expected = NO(T, NO(J, NO(W, L)))
assert result == expected
```

#### Scenario: Handle normal product with Zero operator

**Given** the pyope library is imported
**And** an operator `A` exists
**When** the user calls `normal_product(A, Zero)` or `normal_product(Zero, A)`
**Then** the function returns `Zero` (zero multiplication rule)

**Example**:
```python
from pyope import BasisOperator, normal_product, Zero

T = BasisOperator("T", 2)
result1 = normal_product(T, Zero)
result2 = normal_product(Zero, T)
assert result1 == Zero
assert result2 == Zero
```

#### Scenario: Handle normal product with One operator

**Given** the pyope library is imported
**And** two operators `A` and `B` exist
**When** the user calls `normal_product(A, One, B)`
**Then** the function returns `NO(A, B)` (One is automatically filtered out)

**Example**:
```python
from pyope import BasisOperator, normal_product, NO, One

T = BasisOperator("T", 2)
J = BasisOperator("J", 1)
result = normal_product(T, One, J)
expected = NO(T, J)
assert result == expected
```

#### Scenario: Handle normal product with derivative operators

**Given** the pyope library is imported
**And** derivative operators exist
**When** the user calls `normal_product(d(A), B, C)`
**Then** the function correctly nests derivative operators in the normal product

**Example**:
```python
from pyope import BasisOperator, normal_product, NO, d

T = BasisOperator("T", 2)
J = BasisOperator("J", 1)
result = normal_product(d(T), J, T)
expected = NO(d(T), NO(J, T))
assert result == expected
```

#### Scenario: Handle normal product with scalar coefficients

**Given** the pyope library is imported
**And** operators with scalar coefficients exist
**When** the user calls `normal_product(c*A, k*B)` where `c` and `k` are scalars
**Then** the function returns the product of scalars times the normal product of operators

**Example**:
```python
import sympy as sp
from pyope import BasisOperator, normal_product

c = sp.Symbol("c")
k = sp.Symbol("k")
T = BasisOperator("T", 2)
J = BasisOperator("J", 1)

result = normal_product(c * T, k * J)
# Result should be c*k*NO(T, J) (scalars automatically factored out)
```

#### Scenario: Function must be exported from public API

**Given** the pyope library is installed
**When** the user imports `normal_product` from pyope
**Then** the function is available without importing internal modules

**Example**:
```python
# This should work without importing api.py
from pyope import normal_product
```

#### Scenario: Function must have proper type annotations

**Given** the pyope library is imported
**When** the user inspects the `normal_product` function signature
**Then** the function has type annotations for all parameters and return value

**Example**:
```python
from typing import Any
def normal_product(*operators: Any) -> Any:
    ...
```

#### Scenario: Function must have comprehensive docstring

**Given** the pyope library is imported
**When** the user calls `help(normal_product)`
**Then** the documentation includes:
- A clear description of what the function does
- Parameter documentation
- Return value documentation
- At least 3 usage examples
- Explanation of the conversion rule

**Example**:
```python
help(normal_product)
# Output should show comprehensive documentation
```

#### Scenario: Maintain mathematical correctness with Mathematica

**Given** a set of test operators
**When** computing `normal_product(A, B, C, D)` in Python
**And** computing the equivalent expression in Wolfram Language (OPEdefs.m)
**Then** both results represent the same mathematical object

**Validation**:
- Create parallel test files in Python and Wolfram Script
- Compare outputs for various operator combinations
- Ensure algebraic equivalence

## Related Requirements

### Dependencies
- **NO() function** (existing): `normal_product` uses `NO()` internally
- **One, Zero constants** (existing): For handling special cases
- **BasisOperator class** (existing): For creating test operators

### Cross-References
- This requirement extends the public API defined in `src/pyope/api.py`
- Related to the normal ordering functionality already present
- Supports future work on W-algebras and composite operators

## Non-Requirements

The following are explicitly **out of scope** for this requirement:

1. **Left-associative nesting**: Only right-associative nesting is supported
2. **Custom nesting strategies**: No option to change nesting order
3. **Performance optimization**: Performance should be equivalent to manual NO nesting
4. **Visualization**: No special visualization for normal_product expressions
5. **Operator simplification**: Relies on existing simplify() mechanisms

## Success Metrics

- **Functional correctness**: All test scenarios pass
- **API usability**: Users can use the function without reading internal code
- **Documentation quality**: Docstring and standalone documentation are complete
- **Mathematical correctness**: Results match Mathematica computations
- **Code quality**: Passes mypy, ruff, and black checks
- **Test coverage**: 100% coverage of the new function
