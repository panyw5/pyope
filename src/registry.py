"""
Registry system for operators and OPEs.

This module provides:
- OperatorRegistry: Manages all declared operators
- OPERegistry: Stores and retrieves OPE definitions
"""

from __future__ import annotations
from typing import Dict, Optional, Tuple, Callable
from .operators import Operator, BasisOperator
from .ope_data import OPEData


class OperatorRegistry:
    """
    Registry for all declared operators.

    Maintains a global list of operators with their positions for ordering.
    """

    def __init__(self) -> None:
        """Initialize the operator registry."""
        self._operators: Dict[str, Operator] = {}
        self._position_counter: int = 0

    def register(self, operator: Operator) -> None:
        """
        Register an operator.

        Args:
            operator: Operator to register
        """
        key = self._make_key(operator)

        if key in self._operators:
            # Operator already registered
            return

        # Assign position for ordering
        operator.position = self._position_counter
        self._position_counter += 1

        self._operators[key] = operator

    def _make_key(self, operator: Operator) -> str:
        """
        Create a unique key for an operator.

        Args:
            operator: Operator to create key for

        Returns:
            Unique string key
        """
        if isinstance(operator, BasisOperator) and operator.indices:
            # For indexed operators, use name without indices
            return operator.name
        return operator.name

    def get(self, name: str) -> Optional[Operator]:
        """
        Get an operator by name.

        Args:
            name: Operator name

        Returns:
            Operator if found, None otherwise
        """
        return self._operators.get(name)

    def clear(self) -> None:
        """Clear all registered operators."""
        self._operators.clear()
        self._position_counter = 0

    def __len__(self) -> int:
        """Get number of registered operators."""
        return len(self._operators)


class OPERegistry:
    """
    Registry for OPE definitions.

    Stores user-defined OPEs and provides lookup functionality.
    """

    def __init__(self) -> None:
        """Initialize the OPE registry."""
        self._opes: Dict[Tuple[int, int], OPEData] = {}
        self._ope_functions: Dict[Tuple[int, int], Callable] = {}

    def define(
        self,
        A: Operator,
        B: Operator,
        ope_data: OPEData,
    ) -> None:
        """
        Define an OPE between two operators.

        Args:
            A: First operator
            B: Second operator
            ope_data: OPE result
        """
        key = self._make_key(A, B)
        self._opes[key] = ope_data

    def define_function(
        self,
        A: Operator,
        B: Operator,
        func: Callable,
    ) -> None:
        """
        Define an OPE using a function (for indexed operators).

        Args:
            A: First operator (may have symbolic indices)
            B: Second operator (may have symbolic indices)
            func: Function that computes the OPE given indices
        """
        key = self._make_key(A, B)
        self._ope_functions[key] = func

    def lookup(self, A: Operator, B: Operator) -> Optional[OPEData]:
        """
        Look up an OPE definition.

        Args:
            A: First operator
            B: Second operator

        Returns:
            OPEData if defined, None otherwise
        """
        key = self._make_key(A, B)

        # First check direct definitions
        if key in self._opes:
            return self._opes[key]

        # Then check function definitions (for indexed operators)
        if key in self._ope_functions:
            func = self._ope_functions[key]
            # Call function with operator indices
            result = func(A, B)
            return result

        return None

    def _make_key(self, A: Operator, B: Operator) -> Tuple[int, int]:
        """
        Create a unique key for an operator pair.

        Args:
            A: First operator
            B: Second operator

        Returns:
            Tuple of hashes
        """
        return (hash(A), hash(B))

    def is_defined(self, A: Operator, B: Operator) -> bool:
        """
        Check if an OPE is defined.

        Args:
            A: First operator
            B: Second operator

        Returns:
            True if OPE is defined
        """
        return self.lookup(A, B) is not None

    def clear(self) -> None:
        """Clear all OPE definitions."""
        self._opes.clear()
        self._ope_functions.clear()

    def __len__(self) -> int:
        """Get number of defined OPEs."""
        return len(self._opes) + len(self._ope_functions)


# Global registries
operator_registry = OperatorRegistry()
ope_registry = OPERegistry()


# Helper functions for declaring operators

def Bosonic(*operators: str) -> None:
    """
    Declare one or more bosonic operators.

    Args:
        *operators: Operator names to declare as bosonic

    Example:
        >>> Bosonic('T', 'J')
    """
    for name in operators:
        op = BasisOperator(name, bosonic=True)
        operator_registry.register(op)


def Fermionic(*operators: str) -> None:
    """
    Declare one or more fermionic operators.

    Args:
        *operators: Operator names to declare as fermionic

    Example:
        >>> Fermionic('ψ', 'χ')
    """
    for name in operators:
        op = BasisOperator(name, bosonic=False)
        operator_registry.register(op)
