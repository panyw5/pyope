"""
Operator classes for PyOPE.

This module defines the operator hierarchy used in OPE calculations.
All operators inherit from sympy.Symbol for seamless integration with SymPy.
"""

from __future__ import annotations
from typing import Optional, Union, Tuple, Any
import sympy as sp


class Operator(sp.Symbol):
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

def d(op: Operator, order: int = 1) -> DerivativeOperator:
    """
    Create derivative operator: d(A) = ∂A, d(A, 2) = ∂²A.

    Args:
        op: Operator to differentiate
        order: Order of differentiation (default: 1)

    Returns:
        DerivativeOperator
    """
    return DerivativeOperator(op, order)


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
