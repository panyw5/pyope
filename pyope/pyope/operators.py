"""
Operator classes for PyOPE.

This module defines the operator hierarchy used in OPE calculations.
All operators inherit from LocalOperator for unified type system.
"""

from __future__ import annotations
from typing import Optional, Union, Tuple, Any
import sympy as sp
from .local_operator import LocalOperator, OperatorSum, OperatorProduct


class Operator(sp.Symbol, LocalOperator):
    """
    Base class for all operators in VOA.

    Inherits from sympy.Symbol to enable seamless integration with SymPy's
    expression system. Operators are non-commutative symbols with additional
    properties like parity and position.
    """

    # Class-level storage for operator properties
    _operator_properties: dict = {}

    def __new__(cls, name: str, parity: Union[int, sp.Symbol] = 0, **kwargs):
        """
        Create a new operator.

        Args:
            name: Name of the operator
            parity: Parity of the operator (0=bosonic, 1=fermionic, or symbolic)
            **kwargs: Additional arguments passed to Symbol
        """
        # Create non-commutative symbol
        obj = sp.Symbol.__new__(cls, name, commutative=False, **kwargs)

        # Store properties in class-level dict (since Symbol is immutable)
        obj_id = id(obj)
        cls._operator_properties[obj_id] = {
            'parity': parity,
            'position': None,
        }

        return obj

    @property
    def parity(self) -> Union[int, sp.Symbol]:
        """Get operator parity."""
        return self._operator_properties.get(id(self), {}).get('parity', 0)

    @property
    def position(self) -> Optional[int]:
        """Get operator position (for ordering)."""
        return self._operator_properties.get(id(self), {}).get('position')

    @position.setter
    def position(self, value: int) -> None:
        """Set operator position."""
        obj_id = id(self)
        if obj_id not in self._operator_properties:
            self._operator_properties[obj_id] = {'parity': 0, 'position': None}
        self._operator_properties[obj_id]['position'] = value

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

    # Override arithmetic operations to return LocalOperator types
    def __add__(self, other):
        """Add operators: A + B."""
        if other == 0:
            return self
        return OperatorSum(self, other)

    def __radd__(self, other):
        """Right addition: scalar + A."""
        if other == 0:
            return self
        return OperatorSum(other, self)

    def __sub__(self, other):
        """Subtract operators: A - B."""
        return self + (-1) * other

    def __rsub__(self, other):
        """Right subtraction: scalar - A."""
        return other + (-1) * self

    def __mul__(self, other):
        """Multiply operator by scalar: c * A."""
        if other == 0:
            return sp.S.Zero
        if other == 1:
            return self
        # Only allow scalar multiplication
        if isinstance(other, (int, float, sp.Number, sp.Symbol, sp.Expr)) and not isinstance(other, Operator):
            return OperatorProduct(other, self)
        # For operator * operator, don't define (use NO for normal ordered product)
        return NotImplemented

    def __rmul__(self, other):
        """Right multiplication: scalar * A."""
        return self.__mul__(other)

    def __neg__(self):
        """Negate operator: -A."""
        return (-1) * self


class BasisOperator(Operator):
    """
    Basic operator in the VOA.

    This represents fundamental operators like T (stress-energy tensor),
    J (current), ψ (fermion), etc.
    """

    def __new__(
        cls,
        name: str,
        bosonic: bool = True,
        indexed: bool = False,
        indices: Optional[Tuple[Union[int, sp.Symbol], ...]] = None,
        base_name: Optional[str] = None,
        **kwargs
    ):
        """
        Create a basis operator.

        Args:
            name: Name of the operator (may include indices for display)
            bosonic: True if bosonic, False if fermionic
            indexed: True if this operator can have indices (like J[i])
            indices: Tuple of indices if this is an indexed operator instance
            base_name: Base name without indices (for indexed operators)
            **kwargs: Additional arguments
        """
        parity = 0 if bosonic else 1

        # Create the operator
        obj = Operator.__new__(cls, name, parity, **kwargs)

        # Store additional properties
        obj_id = id(obj)
        cls._operator_properties[obj_id].update({
            'indexed': indexed,
            'indices': indices or (),
            'base_name': base_name or name,
        })

        return obj

    @property
    def indexed(self) -> bool:
        """Check if operator supports indexing."""
        return self._operator_properties.get(id(self), {}).get('indexed', False)

    @property
    def indices(self) -> Tuple[Union[int, sp.Symbol], ...]:
        """Get operator indices."""
        return self._operator_properties.get(id(self), {}).get('indices', ())

    @property
    def base_name(self) -> str:
        """Get base name without indices."""
        return self._operator_properties.get(id(self), {}).get('base_name', self.name)

    def __getitem__(self, index: Union[int, sp.Symbol, Tuple]) -> BasisOperator:
        """
        Create indexed operator: J[i] or J[i, j].

        Args:
            index: Single index or tuple of indices

        Returns:
            New BasisOperator with indices

        Raises:
            TypeError: If operator does not support indexing
        """
        if not self.indexed:
            raise TypeError(f"Operator {self.base_name} does not support indexing")

        # Convert single index to tuple
        if not isinstance(index, tuple):
            index = (index,)

        # Create new operator with indices
        # Use a modified name to distinguish indexed versions
        indices_str = ",".join(str(i) for i in index)
        indexed_name = f"{self.base_name}[{indices_str}]"

        return BasisOperator(
            name=indexed_name,
            bosonic=self.is_bosonic,
            indexed=True,
            indices=index,
            base_name=self.base_name,  # Preserve base name
        )


class DerivativeOperator(Operator):
    """
    Derivative of an operator: ∂A, ∂²A, etc.

    Represents the derivative ∂^n A where n is the order.
    """

    def __new__(cls, base: Operator, order: int = 1, **kwargs):
        """
        Create a derivative operator.

        Args:
            base: The base operator to differentiate
            order: Order of differentiation (default: 1)
            **kwargs: Additional arguments
        """
        if order < 1:
            raise ValueError("Derivative order must be >= 1")

        # Create name for derivative
        if order == 1:
            name = f"∂{base.name}"
        else:
            name = f"∂^{order}{base.name}"

        # Create the operator with same parity as base
        obj = Operator.__new__(cls, name, base.parity, **kwargs)

        # Store derivative-specific properties
        obj_id = id(obj)
        cls._operator_properties[obj_id].update({
            'base': base,
            'order': order,
        })

        return obj

    @property
    def base(self) -> Operator:
        """Get the base operator."""
        return self._operator_properties.get(id(self), {}).get('base')

    @property
    def order(self) -> int:
        """Get the derivative order."""
        return self._operator_properties.get(id(self), {}).get('order', 1)

    def _latex(self, printer) -> str:
        """
        Custom LaTeX representation for derivative operators.

        ∂A is displayed as \partial A
        ∂²A is displayed as \partial^2 A
        """
        base_latex = printer._print(self.base)
        if self.order == 1:
            return f"\\partial {base_latex}"
        else:
            return f"\\partial^{{{self.order}}} {base_latex}"


class NormalOrderedOperator(Operator):
    """
    Normal ordered product of two operators: NO[A, B] = (AB)(z).

    Represents the normal ordered product using point-splitting convention.
    """

    def __new__(cls, left: Operator, right: Operator, **kwargs):
        """
        Create a normal ordered operator.

        Args:
            left: Left operator
            right: Right operator
            **kwargs: Additional arguments
        """
        name = f"NO[{left.name},{right.name}]"

        # Parity of NO[A,B] is sum of parities
        if isinstance(left.parity, int) and isinstance(right.parity, int):
            parity = (left.parity + right.parity) % 2
        else:
            parity = left.parity + right.parity  # Symbolic

        # Create the operator
        obj = Operator.__new__(cls, name, parity, **kwargs)

        # Store NO-specific properties
        obj_id = id(obj)
        cls._operator_properties[obj_id].update({
            'left': left,
            'right': right,
        })

        return obj

    @property
    def left(self) -> Operator:
        """Get the left operator."""
        return self._operator_properties.get(id(self), {}).get('left')

    @property
    def right(self) -> Operator:
        """Get the right operator."""
        return self._operator_properties.get(id(self), {}).get('right')

    @property
    def factors(self) -> Tuple[Operator, Operator]:
        """Get both factors as a tuple."""
        return (self.left, self.right)


# Helper functions for creating operators

def d(op, order: int = 1):
    """
    Create derivative operator or differentiate an expression.

    - If op is an Operator: d(A) = ∂A, d(A, 2) = ∂²A
    - If op is a DerivativeOperator: d(∂^n A) = ∂^(n+order) A (合并阶数)
    - If op is an expression: d(expr) applies operator_derivative to expr

    Args:
        op: Operator or expression to differentiate
        order: Order of differentiation (default: 1)

    Returns:
        DerivativeOperator if op is Operator, otherwise differentiated expression
    """
    # Import here to avoid circular dependency
    from .api import operator_derivative

    if isinstance(op, DerivativeOperator):
        # 合并导数阶数: d(∂^n A, m) = ∂^(n+m) A
        return DerivativeOperator(op.base, op.order + order)
    elif isinstance(op, Operator):
        return DerivativeOperator(op, order)
    else:
        # For expressions, use operator_derivative
        return operator_derivative(op, order)


def dn(order: int, op: Operator) -> DerivativeOperator:
    """
    Create derivative operator with order first: dn(2, A) = ∂²A.

    Args:
        order: Order of differentiation
        op: Operator to differentiate

    Returns:
        DerivativeOperator
    """
    return DerivativeOperator(op, order)
