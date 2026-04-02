# Test Conventions

> Test structure, fixtures, naming, and markers for pyope.

---

## Directory Structure

```
tests/
├── conftest.py                    # Global fixtures: cache disable, registry clear, sl2/w3 algebras
├── test_virasoro_voa.py           # Virasoro algebra (complete test suite pattern)
├── test_w3_algebra.py             # W₃ algebra
├── test_w3_algebra_ref.py         # W₃ with Mathematica reference data
├── test_sl2_k.py                  # sl(2) Kac-Moody
├── test_bc_betagamma.py           # bc/βγ ghost systems
├── test_simplify.py               # Simplification engine
├── test_jacobi.py                 # Jacobi identity verification
├── test_operator_spaces.py        # Operator space enumeration, C2, realization
├── test_c2_space.py               # C2 quotient space
├── test_c2_null_searcher.py       # Null state search
├── test_descendants.py            # Descendant generation
├── test_singularity.py            # Singular vector analysis
└── ...
```

---

## Test Class Organization

Group tests by **aspect**, not by function:

```python
class TestW3AlgebraDefinition:
    """Operator declarations and OPE rule setup."""

class TestW3AlgebraComputations:
    """Computed OPE results (must reference OPEdefs.m)."""

class TestW3AlgebraProperties:
    """Algebraic identities and consistency checks."""
```

---

## Naming Convention

```python
def test_<subject>_<aspect>(self):
    """
    <Human-readable description>

    Reference: <OPEdefs.m file, line, or published paper>
    """
```

Examples:
- `test_virasoro_ope_with_central_charge` — definition test
- `test_t_derivative_ope` — computation test (needs Mathematica ref)
- `test_primary_field_condition` — property test
- `test_conformal_weight_additivity` — property test

---

## Fixtures

### Global (conftest.py, autouse)

```python
@pytest.fixture(autouse=True, scope="function")
def disable_cache_for_tests():
    """Disables cache and clears registry before/after each test."""
    cache = get_ope_cache()
    cache.disable()
    cache.clear()
    ope_registry.clear()
    yield
    cache.enable()
    cache.clear()
    ope_registry.clear()
```

**This runs for every test automatically.** You do NOT need to clear the registry manually.

### Algebra Fixtures (conftest.py)

```python
@pytest.fixture
def sl2_algebra():
    """sl(2) Kac-Moody at level k. Returns dict with Jplus, Jzero, Jminus, k, clear."""

@pytest.fixture
def w3_algebra():
    """W₃ algebra. Returns dict with T, W, c, beta, Lambda, clear."""
```

### Per-File Fixtures

For algebras specific to one test file, define a local `autouse` fixture:

```python
@pytest.fixture(autouse=True)
def clear_registry():
    from pyope.registry import ope_registry
    ope_registry.clear()
    yield
    ope_registry.clear()
```

---

## Test Markers

```python
@pytest.mark.slow              # Long-running tests (use -m "not slow" to skip)
@pytest.mark.mathematica_ref   # Tests with OPEdefs.m reference data
@pytest.mark.requires_derivative  # Tests that need d()/dn() support
```

---

## Standard Test Template

```python
"""
<Algebra Name> Tests

Reference: OPEdefs/<file>.nb line <N>-<M>
"""

import pytest
import sympy as sp
from sympy import Symbol, Rational, simplify

from pyope import (
    BasisOperator, Bosonic, Fermionic,
    OPE, NO, bracket, MakeOPE,
    d, dn, One, Zero,
)


class TestAlgebraDefinition:
    def test_operator_declarations(self):
        """Verify operator types, weights, and parities."""
        ...

class TestAlgebraComputations:
    def test_some_ope(self):
        """
        OPE(A, B) verification.

        Reference: OPEdefs.m computation
        Mathematica: OPE[A, B] // Simplify
        Expected: [c-number pole, ..., operator pole]
        """
        ...
        # Assert against Mathematica reference
        assert simplify(result.pole(n) - expected) == 0

class TestAlgebraProperties:
    def test_some_identity(self):
        """Verify algebraic identity (e.g., Jacobi, associativity)."""
        ...
```

---

## Running Tests

```bash
# All tests
pytest tests/ -v

# Fast tests only
pytest tests/ -v -m "not slow"

# Single file
pytest tests/test_virasoro_voa.py -v -s

# With Mathematica reference tests highlighted
pytest tests/ -v -m "mathematica_ref"
```
