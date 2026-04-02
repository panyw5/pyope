# Algebra Fixtures

> Standard VOA algebra definitions and fixture reuse patterns.

---

## Available Fixtures (conftest.py)

### sl2_algebra

**sl(2) Kac-Moody** at level $k$. Three generators $J^+, J^0, J^-$ with:

$$J^0(z) J^0(w) \sim \frac{k/2}{(z-w)^2}, \quad J^+(z) J^-(w) \sim \frac{k}{(z-w)^2} + \frac{2J^0}{(z-w)}$$

```python
def test_something(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    J0 = sl2_algebra['Jzero']
    Jm = sl2_algebra['Jminus']
    k  = sl2_algebra['k']
```

### w3_algebra

**$W_3$ algebra** with stress tensor $T$ (weight 2) and primary $W$ (weight 3).

Includes auxiliary operator $\Lambda = NO(T,T) - \frac{3}{10}\partial^2 T$.

```python
def test_something(w3_algebra):
    T      = w3_algebra['T']
    W      = w3_algebra['W']
    c      = w3_algebra['c']
    beta   = w3_algebra['beta']
    Lambda = w3_algebra['Lambda']
```

---

## Defining New Algebra Fixtures

### When to Create a Fixture

- The algebra is used in **3+ tests** across the file
- The algebra setup is **more than 10 lines**
- The algebra is a **standard VOA** (Virasoro, Kac-Moody, W-algebra, ghost system, ...)

### Pattern

```python
@pytest.fixture
def my_algebra():
    """
    <Algebra name> fixture.

    Generators: <list>
    OPE rules: <brief description>
    Reference: <source>
    """
    # 1. Declare operators
    A = BasisOperator("A", conformal_weight=h_a)
    B = BasisOperator("B", conformal_weight=h_b)

    # 2. Register parity
    Bosonic(A, B)  # or Fermionic(...)

    # 3. Define parameters
    c = sp.Symbol("c")

    # 4. Define OPE rules
    OPE[A, B] = MakeOPE([...])

    # 5. Return dict
    return {
        "A": A,
        "B": B,
        "c": c,
    }
```

### Fixture Scope

- **Default** (`scope="function"`): Each test gets a fresh algebra. Safe but slow for large algebras.
- **Do NOT use** `scope="session"` or `scope="module"` — the autouse `disable_cache_for_tests` fixture clears the registry per function.

---

## Standard Algebras Reference

For implementing new fixtures, here are the standard VOA definitions:

### Virasoro

- Generator: $T$ (weight 2, bosonic)
- OPE: $T(z)T(w) \sim \frac{c/2}{(z-w)^4} + \frac{2T}{(z-w)^2} + \frac{\partial T}{(z-w)}$

### Free Boson

- Generator: $\partial\phi$ (weight 1, bosonic)
- OPE: $\partial\phi(z)\partial\phi(w) \sim \frac{1}{(z-w)^2}$

### $bc$ Ghost System

- Generators: $b$ (weight $\lambda$, fermionic), $c$ (weight $1-\lambda$, fermionic)
- OPE: $b(z)c(w) \sim \frac{1}{z-w}$

### $\beta\gamma$ Ghost System

- Generators: $\beta$ (weight $\lambda$, bosonic), $\gamma$ (weight $1-\lambda$, bosonic)
- OPE: $\beta(z)\gamma(w) \sim \frac{-1}{z-w}$
- Note: the minus sign for bosonic ghosts

### sl(2) Kac-Moody

- Generators: $J^+, J^0, J^-$ (all weight 1, bosonic)
- OPEs: see `conftest.py` `sl2_algebra` fixture

### $W_3$ Algebra

- Generators: $T$ (weight 2), $W$ (weight 3)
- OPEs: see `conftest.py` `w3_algebra` fixture

---

## Anti-Patterns

### Don't duplicate algebra setup

```python
# BAD: repeating Virasoro setup in every test
class TestSomething:
    def test_a(self):
        T = BasisOperator("T", conformal_weight=2)
        Bosonic(T)
        c = Symbol("c")
        OPE[T, T] = MakeOPE([c/2*One, Zero, 2*T, d(T)])
        # ... test ...

    def test_b(self):
        T = BasisOperator("T", conformal_weight=2)  # same setup again!
        # ...

# GOOD: use a fixture or class-level setup
@pytest.fixture
def virasoro():
    T = BasisOperator("T", conformal_weight=2)
    Bosonic(T)
    c = sp.Symbol("c")
    OPE[T, T] = MakeOPE([c/2*One, Zero, 2*T, d(T)])
    return {"T": T, "c": c}
```

### Don't test OPE definition correctness — test computation

```python
# WEAK: just checking what you put in comes back out
OPE[T, T] = MakeOPE([c/2*One, Zero, 2*T, d(T)])
assert OPE(T, T).pole(4) == c/2*One  # This just tests MakeOPE roundtrip

# STRONG: test a DERIVED computation
# OPE(T, d(T)) is computed from OPE(T,T) via the derivative formula
# Its result is non-trivial and validates the computation engine
result = OPE(T, d(T))
assert simplify(result.pole(5) - 2*c*One) == 0  # Derived, not stored
```
