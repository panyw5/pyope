"""
Local operator hierarchy for PyOPE.

This module defines a unified LocalOperator class hierarchy that represents
all local operators in conformal field theory and vertex operator algebras.

Design:
    LocalOperator (abstract base)
    ├── BasisOperator (fundamental operators/generators)
    ├── DerivativeOperator (derivatives of operators)
    ├── NormalOrderedOperator (normal ordered products)
    ├── ConstantOperator (constant operators like One)
    └── CompositeOperator (linear combinations)
        ├── OperatorSum (A + B)
        └── OperatorProduct (c*A, where c is scalar)
"""

from __future__ import annotations
from typing import Union, Optional, Tuple, Any
import sympy as sp
from abc import ABC, abstractmethod


class LocalOperator(sp.Expr, ABC):
    """
    Abstract base class for all local operators.

    A local operator is a field operator in conformal field theory that
    acts at a specific point. This includes:
    - Fundamental operators (generators)
    - Derivatives of operators
    - Normal ordered products
    - Linear combinations of operators

    All local operators support:
    - Algebraic operations (+, -, *, /)
    - OPE computation
    - Parity (bosonic/fermionic)
    """

    def __new__(cls, *args, **kwargs):
        """Create a new local operator."""
        obj = sp.Expr.__new__(cls)
        return obj

    @property
    @abstractmethod
    def parity(self) -> Union[int, sp.Symbol]:
        """
        Get operator parity.

        Returns:
            0 for bosonic, 1 for fermionic, or symbolic expression
        """
        pass

    @property
    def is_bosonic(self) -> bool:
        """Check if operator is bosonic."""
        parity = self.parity
        if isinstance(parity, int):
            return parity % 2 == 0
        return False  # Symbolic parity

    @property
    def is_fermionic(self) -> bool:
        """Check if operator is fermionic."""
        parity = self.parity
        if isinstance(parity, int):
            return parity % 2 == 1
        return False  # Symbolic parity

    def __add__(self, other) -> LocalOperator:
        """Add two operators: A + B."""
        if other == 0:
            return self
        if isinstance(other, LocalOperator):
            return OperatorSum(self, other)
        # Scalar + Operator
        return OperatorSum(self, other)

    def __radd__(self, other) -> LocalOperator:
        """Right addition: scalar + A."""
        if other == 0:
            return self
        return OperatorSum(other, self)

    def __sub__(self, other) -> LocalOperator:
        """Subtract operators: A - B."""
        return self + (-1) * other

    def __rsub__(self, other) -> LocalOperator:
        """Right subtraction: scalar - A."""
        return other + (-1) * self

    def __mul__(self, other) -> LocalOperator:
        """Multiply operator by scalar: c * A."""
        if other == 0:
            return sp.S.Zero
        if other == 1:
            return self
        if isinstance(other, (int, float, sp.Number, sp.Symbol)):
            return OperatorProduct(other, self)
        # For operator * operator, we don't define it here
        # (use NO for normal ordered product)
        return NotImplemented

    def __rmul__(self, other) -> LocalOperator:
        """Right multiplication: scalar * A."""
        return self.__mul__(other)

    def __neg__(self) -> LocalOperator:
        """Negate operator: -A."""
        return (-1) * self

    def _latex(self, printer) -> str:
        """LaTeX representation (to be overridden by subclasses)."""
        return str(self)


class CompositeOperator(LocalOperator):
    """
    Base class for composite operators (sums and products).
    """

    @property
    def parity(self) -> Union[int, sp.Symbol]:
        """Parity of composite operator (to be computed by subclasses)."""
        return 0


class OperatorSum(CompositeOperator):
    """
    Sum of operators: A + B + C + ...

    Represents linear combinations of local operators.
    """

    def __new__(cls, *terms):
        """
        Create operator sum.

        Args:
            *terms: Operators or scalars to sum
        """
        # Flatten nested sums and collect terms
        flattened = []
        for term in terms:
            if isinstance(term, OperatorSum):
                flattened.extend(term.terms)
            elif term != 0:
                flattened.append(term)

        # If empty, return zero
        if not flattened:
            return sp.S.Zero

        # If single term, return it directly
        if len(flattened) == 1:
            return flattened[0]

        # Create sum
        obj = LocalOperator.__new__(cls)
        obj._terms = tuple(flattened)
        return obj

    @property
    def terms(self) -> Tuple:
        """Get the terms of the sum."""
        return getattr(self, '_terms', ())

    @property
    def args(self) -> Tuple:
        """Get args for SymPy compatibility."""
        return self.terms

    @property
    def parity(self) -> Union[int, sp.Symbol]:
        """
        Parity of sum.

        For a sum to have definite parity, all terms must have the same parity.
        Otherwise, it's undefined (return 0 as default).
        """
        if not self.terms:
            return 0

        # Get parity of first term
        first_parity = None
        for term in self.terms:
            if isinstance(term, LocalOperator):
                first_parity = term.parity
                break

        if first_parity is None:
            return 0

        # Check if all terms have same parity
        for term in self.terms:
            if isinstance(term, LocalOperator):
                if term.parity != first_parity:
                    return 0  # Mixed parity

        return first_parity

    def __str__(self) -> str:
        """String representation."""
        return " + ".join(str(term) for term in self.terms)

    def __repr__(self) -> str:
        """Repr representation."""
        return f"OperatorSum({', '.join(repr(term) for term in self.terms)})"

    def _latex(self, printer) -> str:
        """LaTeX representation."""
        return " + ".join(printer._print(term) for term in self.terms)

    def _sympystr(self, printer) -> str:
        """SymPy string representation."""
        return str(self)


class OperatorProduct(CompositeOperator):
    """
    Scalar multiple of operator: c * A

    Represents a scalar coefficient times an operator.
    """

    def __new__(cls, coeff, operator):
        """
        Create operator product.

        Args:
            coeff: Scalar coefficient
            operator: Operator to multiply
        """
        # Handle special cases
        if coeff == 0:
            return sp.S.Zero
        if coeff == 1:
            return operator

        # If operator is already a product, combine coefficients
        if isinstance(operator, OperatorProduct):
            new_coeff = coeff * operator.coeff
            return cls.__new__(cls, new_coeff, operator.operator)

        # If operator is a sum, distribute
        if isinstance(operator, OperatorSum):
            return OperatorSum(*(coeff * term for term in operator.terms))

        # Create product
        obj = LocalOperator.__new__(cls)
        obj._coeff = coeff
        obj._operator = operator
        return obj

    @property
    def coeff(self):
        """Get the coefficient."""
        return getattr(self, '_coeff', sp.S.One)

    @property
    def operator(self):
        """Get the operator."""
        return getattr(self, '_operator', None)

    @property
    def args(self) -> Tuple:
        """Get args for SymPy compatibility."""
        # Convert coefficient to SymPy object if needed
        coeff = self.coeff
        if not isinstance(coeff, sp.Basic):
            coeff = sp.sympify(coeff)
        return (coeff, self.operator)

    @property
    def parity(self) -> Union[int, sp.Symbol]:
        """Parity of product (same as operator parity)."""
        if isinstance(self.operator, LocalOperator):
            return self.operator.parity
        return 0

    def __str__(self) -> str:
        """String representation."""
        return f"{self.coeff}*{self.operator}"

    def __repr__(self) -> str:
        """Repr representation."""
        return f"OperatorProduct({repr(self.coeff)}, {repr(self.operator)})"

    def _latex(self, printer) -> str:
        """LaTeX representation."""
        coeff_str = printer._print(self.coeff)
        op_str = printer._print(self.operator)
        return f"{coeff_str} {op_str}"

    def _sympystr(self, printer) -> str:
        """SymPy string representation."""
        return str(self)


# Helper function to check if something is a local operator
def is_local_operator(obj) -> bool:
    """
    Check if an object is a local operator.

    Args:
        obj: Object to check

    Returns:
        True if obj is a LocalOperator
    """
    return isinstance(obj, LocalOperator)


# Helper function to extract scalar and operator from expression
def extract_scalar_operator(expr) -> Tuple[sp.Expr, Optional[LocalOperator]]:
    """
    Extract scalar coefficient and operator from expression.

    Args:
        expr: Expression to analyze

    Returns:
        Tuple of (scalar, operator) or (expr, None) if not a product
    """
    if isinstance(expr, OperatorProduct):
        return expr.coeff, expr.operator
    elif isinstance(expr, LocalOperator):
        return sp.S.One, expr
    else:
        return expr, None
