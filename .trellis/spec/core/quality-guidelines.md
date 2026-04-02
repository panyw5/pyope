# Quality Guidelines

> Code standards, forbidden patterns, and review checklist for pyope.

---

## Code Style

- **Formatter**: `black` (line length 88)
- **Linter**: `ruff`
- **Type checker**: `mypy`
- **Test runner**: `pytest`

```bash
# Run quality checks
black --check src/ tests/
ruff check src/ tests/
mypy src/pyope/
pytest tests/ -x
```

---

## Forbidden Patterns

### F1: Direct operator multiplication

```python
# FORBIDDEN: raises IllegalOperatorProductError at runtime
result = T * W
result = A * B * C

# CORRECT: use NO() for normal-ordered products
result = NO(T, W)
result = NO(A, NO(B, C))
```

**Why**: In VOA, operator "multiplication" is not well-defined without specifying ordering. `NO(A, B)` makes the normal ordering explicit.

### F2: Hand-written expected values in computation tests

```python
# FORBIDDEN: inventing expected OPE results
def test_some_ope():
    result = OPE(W, W)
    assert result.pole(4) == "something I think is right"  # NO!

# CORRECT: derive from OPEdefs.m or known CFT identities
# See testing/mathematica-reference.md for the proper workflow
```

**Why**: VOA computations are subtle. Expected values must come from Thielemans' OPEdefs.m computations or published literature.

### F3: Assuming cache state

```python
# FORBIDDEN: directly manipulating cache internals
from pyope.cache import _cache_dict
_cache_dict.clear()

# CORRECT: use the public API
cache = get_ope_cache()
cache.clear()
```

### F4: Importing across layers upward

```python
# FORBIDDEN: Foundation importing from Core
# In operators.py:
from .api import OPE  # NO! operators.py is Foundation, api.py is Core

# If truly needed, use lazy import inside a function body
```

### F5: Using Python float for mathematical constants

```python
# FORBIDDEN
coeff = 0.3  # floating point

# CORRECT
coeff = sp.Rational(3, 10)  # exact
```

---

## Required Patterns

### R1: Register parity before OPE definition

```python
T = BasisOperator("T", conformal_weight=2)
Bosonic(T)                              # ← REQUIRED
OPE[T, T] = MakeOPE([...])
```

### R2: Test isolation via registry clear

Every test file should either:
- Use the `autouse` fixture from `conftest.py` (preferred)
- Or manually clear: `ope_registry.clear()` in setup/teardown

### R3: Document Mathematica reference in test docstrings

```python
def test_w3_ww_ope(self):
    """
    W(z)W(w) OPE verification.

    Reference: OPEdefs/ope-examples.nb line 716-971
    Mathematica command: OPE[W, W]
    """
```

### R4: Type stubs for public modules

Any module that exports public API should have a `.pyi` stub file.

---

## Review Checklist

Before merging any change:

- [ ] `black --check` passes
- [ ] `ruff check` passes
- [ ] `mypy src/pyope/` passes (or known exceptions documented)
- [ ] `pytest tests/ -x` passes
- [ ] No new forbidden patterns introduced
- [ ] If new public API: `.pyi` stub updated
- [ ] If new computation: Mathematica reference documented
- [ ] If touching operator types: `IllegalOperatorProductError` guard preserved
