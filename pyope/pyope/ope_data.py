"""
OPEData class for storing Operator Product Expansion results.

This module defines the OPEData class which stores the poles of an OPE
in a sparse format for efficient computation.
"""

from __future__ import annotations
from typing import Dict, Optional, Union, Callable, Any
import sympy as sp
from .operators import Operator


class OPEData:
    """
    Stores the result of an Operator Product Expansion.

    An OPE A(z)B(w) is represented as a Laurent series:
        A(z)B(w) = sum_{n=1}^{max_pole} [AB]_n / (z-w)^n + regular terms

    This class stores the poles [AB]_n in a sparse dictionary format.
    """

    def __init__(self, poles: Optional[Dict[int, sp.Expr]] = None):
        """
        Initialize OPEData.

        Args:
            poles: Dictionary mapping pole order to operator expression
                   {n: [AB]_n} where n >= 1
        """
        self._poles: Dict[int, sp.Expr] = poles or {}
        self._max_pole: Optional[int] = None
        self._update_max_pole()

    def _update_max_pole(self) -> None:
        """Update the cached maximum pole order."""
        if self._poles:
            self._max_pole = max(self._poles.keys())
        else:
            self._max_pole = None

    @classmethod
    def from_list(cls, pole_list: list) -> OPEData:
        """
        Create OPEData from a list of poles (Mathematica style).

        The list should be ordered from highest pole to lowest:
        [pole_n, pole_{n-1}, ..., pole_1]

        Args:
            pole_list: List of pole expressions

        Returns:
            OPEData object
        """
        n = len(pole_list)
        poles = {}
        for i, pole in enumerate(pole_list):
            pole_order = n - i
            if pole != 0:  # Only store non-zero poles
                poles[pole_order] = pole
        return cls(poles)

    def pole(self, n: int) -> sp.Expr:
        """
        Get the n-th order pole.

        Args:
            n: Pole order (n >= 1 for poles, n <= 0 for regular terms)

        Returns:
            Operator expression at pole n, or 0 if not present
        """
        if n <= 0:
            # Regular terms not yet implemented
            return sp.S.Zero
        return self._poles.get(n, sp.S.Zero)

    def set_pole(self, n: int, value: sp.Expr) -> None:
        """
        Set the n-th order pole.

        Args:
            n: Pole order (must be >= 1)
            value: Operator expression
        """
        if n < 1:
            raise ValueError("Pole order must be >= 1")

        if value == 0:
            # Remove zero poles
            self._poles.pop(n, None)
        else:
            self._poles[n] = value

        self._update_max_pole()

    @property
    def max_pole(self) -> int:
        """Get the maximum pole order."""
        return self._max_pole if self._max_pole is not None else 0

    @property
    def poles(self) -> Dict[int, sp.Expr]:
        """Get the poles dictionary (read-only view)."""
        return self._poles.copy()

    def is_zero(self) -> bool:
        """Check if this OPE is zero (no poles)."""
        return len(self._poles) == 0

    # Arithmetic operations

    def __add__(self, other: Union[OPEData, int, float]) -> OPEData:
        """
        Add two OPEData objects or add a scalar.

        Args:
            other: Another OPEData or scalar

        Returns:
            New OPEData with combined poles
        """
        if isinstance(other, (int, float)) and other == 0:
            return self

        if not isinstance(other, OPEData):
            raise TypeError(f"Cannot add OPEData with {type(other)}")

        # Combine poles
        result_poles = self._poles.copy()
        for n, pole in other._poles.items():
            if n in result_poles:
                # Use expand instead of simplify for non-commutative expressions
                combined = result_poles[n] + pole
                try:
                    combined = sp.expand(combined)
                except:
                    pass  # Keep unsimplified if expand fails

                if combined == 0:
                    del result_poles[n]
                else:
                    result_poles[n] = combined
            else:
                result_poles[n] = pole

        return OPEData(result_poles)

    def __radd__(self, other: Union[int, float]) -> OPEData:
        """Right addition (for sum() support)."""
        return self.__add__(other)

    def __mul__(self, scalar: Union[int, float, sp.Expr]) -> OPEData:
        """
        Multiply OPEData by a scalar.

        Args:
            scalar: Scalar multiplier

        Returns:
            New OPEData with scaled poles
        """
        if scalar == 0:
            return OPEData()

        if scalar == 1:
            return self

        result_poles = {}
        for n, pole in self._poles.items():
            scaled = scalar * pole
            try:
                scaled = sp.expand(scaled)
            except:
                pass  # Keep unsimplified if expand fails
            result_poles[n] = scaled
        return OPEData(result_poles)

    def __rmul__(self, scalar: Union[int, float, sp.Expr]) -> OPEData:
        """Right multiplication."""
        return self.__mul__(scalar)

    def __neg__(self) -> OPEData:
        """Negation."""
        return self.__mul__(-1)

    def __sub__(self, other: OPEData) -> OPEData:
        """Subtraction."""
        return self.__add__(other.__neg__())

    # Comparison

    def __eq__(self, other: Any) -> bool:
        """Equality comparison."""
        if not isinstance(other, OPEData):
            return False

        # Check if all poles are equal
        all_orders = set(self._poles.keys()) | set(other._poles.keys())
        for n in all_orders:
            diff = self.pole(n) - other.pole(n)
            try:
                diff = sp.expand(diff)
            except:
                pass  # Keep unsimplified if expand fails
            if diff != 0:
                return False
        return True

    # Utility methods

    def simplify(self, func: Optional[Callable] = None) -> OPEData:
        """
        Simplify all poles using a given function.

        Args:
            func: Simplification function (default: sympy.simplify)

        Returns:
            New OPEData with simplified poles
        """
        if func is None:
            func = sp.simplify

        result_poles = {}
        for n, pole in self._poles.items():
            simplified = func(pole)
            if simplified != 0:
                result_poles[n] = simplified

        return OPEData(result_poles)

    def to_series(
        self, z: sp.Symbol, w: sp.Symbol, max_order: int = 0
    ) -> sp.series:
        """
        Convert to a SymPy series expansion.

        Args:
            z: First variable
            w: Second variable
            max_order: Maximum order to include (0 = only poles)

        Returns:
            SymPy series object
        """
        if self.is_zero():
            return sp.O((z - w) ** max_order)

        # Build series from poles
        terms = []
        for n in sorted(self._poles.keys(), reverse=True):
            pole = self._poles[n]
            # Substitute w for the argument in operators
            # This is a simplified version - full implementation needs operator evaluation
            terms.append(pole / (z - w) ** n)

        result = sum(terms) if terms else sp.S.Zero
        return result + sp.O((z - w) ** max_order)

    def __repr__(self) -> str:
        """String representation."""
        if self.is_zero():
            return "OPEData({})"

        pole_strs = []
        for n in sorted(self._poles.keys(), reverse=True):
            pole_strs.append(f"{n}: {self._poles[n]}")

        return f"OPEData({{{', '.join(pole_strs)}}})"

    def __str__(self) -> str:
        """User-friendly string representation."""
        if self.is_zero():
            return "<< 0 >>"

        pole_strs = []
        for n in sorted(self._poles.keys(), reverse=True):
            pole_strs.append(f"{n}|| {self._poles[n]}")

        return f"<< {' '.join(pole_strs)} >>"
