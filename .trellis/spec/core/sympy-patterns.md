# SymPy Patterns

> SymPy-specific coding patterns, idioms, and pitfalls for pyope development.

---

## Core Principle

All pyope operators are **SymPy Symbol subclasses**. Expressions are SymPy `Expr` trees. This gives us free symbolic algebra but imposes SymPy's rules.

---

## Coefficient Handling

### Use `sympy.Rational`, not Python fractions

```python
# BAD: Python float (loses precision)
coeff = 3/10

# GOOD: SymPy exact rational
coeff = sp.Rational(3, 10)
```

### Use `sp.Symbol` for parameters

```python
# Central charge, level, etc. are sympy Symbols
c = sp.Symbol("c")
k = sp.Symbol("k")
beta = sp.Symbol("beta")
```

### Sympify external values

```python
from pyope.utils import sympify_coefficient

# When receiving coefficients from external sources:
coeff = sympify_coefficient(user_input)
```

---

## Expression Comparison

### Equality with simplification

```python
# BAD: direct == may fail for equivalent expressions
assert result == expected  # Can fail if forms differ

# GOOD: use sympy.simplify for structural differences
from sympy import simplify, expand
assert simplify(result - expected) == 0

# GOOD: expand first for polynomial expressions
assert expand(result - expected) == 0
```

### Checking for zero

```python
# For operator expressions:
if expr == 0:      # works for literal zero
if expr is Zero:   # works for pyope Zero constant
if expr == Zero:   # also works

# For sympy expressions with parameters:
from sympy import simplify
assert simplify(expr) == 0
```

---

## Performance Pitfalls

### Matrix rank with sympy vs numpy

```python
# BAD: sympy Matrix.rank() can HANG on large rational matrices
M = sp.Matrix(large_data)
rank = M.rank()  # ← may never return

# GOOD: convert to numpy for numerical rank
import numpy as np
M_np = np.array(M.tolist(), dtype=float)
rank = np.linalg.matrix_rank(M_np)
```

### Avoid deep expression nesting

```python
# Deeply nested NO expressions (e.g., NO(G, NO(GW, NO(...))))
# become extremely expensive to realize (100s+ per monomial at weight 8+)
# Strategy: use C2(F) syntactic checks before full realization when possible
```

### Cache-aware coding

```python
from pyope.cache import get_ope_cache

cache = get_ope_cache()
# Cache is enabled by default for production use
# Cache is disabled in tests (via conftest.py fixture)
# Never assume cache state — always go through get_ope_cache()
```

---

## Creating New Operator Types

If you need a new operator subclass:

```python
class MyOperator(Operator):
    """Custom operator type."""

    def __new__(cls, name, my_param, **assumptions):
        obj = Symbol.__new__(cls, name, **assumptions)
        obj._my_param = my_param
        return obj

    @property
    def parity(self):
        return ope_registry.get_parity(self)

    @property
    def conformal_weight(self):
        return self._my_param  # or compute from stored data
```

**Key**: Override `__new__`, not `__init__` (SymPy Symbol convention). Store custom data as `_prefixed` attributes.

---

## Expression Tree Walking

```python
# Check if expression contains operators
expr.has(Operator)

# Extract scalar coefficient from scalar * operator
from pyope.local_operator import extract_scalar_operator
scalar, op = extract_scalar_operator(term)

# Check if an expression is a local operator
from pyope.local_operator import is_local_operator
if is_local_operator(expr):
    ...
```

---

## Common SymPy Gotchas in pyope

1. **`Mul` evaluation**: SymPy auto-simplifies `Mul(2, T)` to `2*T` but `Mul(T, W)` stays as product. Since `T * W` is forbidden (raises error), scalar-operator products are always `scalar * single_operator`.

2. **`Add` with Zero**: `T + 0` gives `T` in SymPy, but `T + Zero` may not simplify automatically since `Zero` is a special `ConstantOperator`. The `__add__` override on `Operator` handles this.

3. **`==` vs `is`**: SymPy symbols with the same name are `==` but not `is`. Always use `==` for comparison.

4. **`sp.S.Zero` vs `Zero`**: `sp.S.Zero` is SymPy's numeric zero. `Zero` is pyope's operator zero (`ConstantOperator`). They are `==` but different types.
