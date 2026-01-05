"""
Core API functions for PyOPE.

This module provides the main user-facing functions:
- OPE: Compute operator product expansions
- NO: Normal ordered products
- bracket: Commutator/anticommutator
"""

from __future__ import annotations
from typing import Union, Optional
import sympy as sp
from .operators import Operator, NormalOrderedOperator, DerivativeOperator
from .ope_data import OPEData
from .registry import ope_registry


class OPEFunction:
    """
    Main OPE computation function.

    This class implements the OPE computation with automatic linearity expansion.
    """

    def __call__(
        self, A: Union[Operator, sp.Expr], B: Union[Operator, sp.Expr]
    ) -> OPEData:
        """
        Compute OPE of A(z) with B(w).

        Args:
            A: First operator or expression
            B: Second operator or expression

        Returns:
            OPEData containing the poles of the OPE
        """
        # Handle linearity: OPE(A, B + C) = OPE(A, B) + OPE(A, C)
        if isinstance(B, sp.Add):
            return sum(self(A, term) for term in B.args)

        if isinstance(A, sp.Add):
            return sum(self(term, B) for term in A.args)

        # Handle scalar multiplication: OPE(c*A, B) = c * OPE(A, B)
        if isinstance(B, sp.Mul):
            scalar, op = self._extract_scalar(B)
            if scalar is not None:
                return scalar * self(A, op)

        if isinstance(A, sp.Mul):
            scalar, op = self._extract_scalar(A)
            if scalar is not None:
                return scalar * self(op, B)

        # Handle zero
        if A == 0 or B == 0:
            return OPEData()

        # Now A and B should be operators
        if not isinstance(A, Operator):
            raise TypeError(f"Expected Operator, got {type(A)}")
        if not isinstance(B, Operator):
            raise TypeError(f"Expected Operator, got {type(B)}")

        # Dispatch to appropriate computation method
        return self._compute(A, B)

    def _extract_scalar(
        self, expr: sp.Expr
    ) -> tuple[Optional[sp.Expr], Optional[Operator]]:
        """
        Extract scalar coefficient from expression.

        Args:
            expr: Expression to analyze

        Returns:
            Tuple of (scalar, operator) or (None, None)
        """
        if not isinstance(expr, sp.Mul):
            return None, None

        # Separate scalar and operator parts
        scalar_parts = []
        operator_part = None

        for arg in expr.args:
            if isinstance(arg, Operator):
                if operator_part is not None:
                    # Multiple operators, can't extract scalar
                    return None, None
                operator_part = arg
            else:
                scalar_parts.append(arg)

        if operator_part is None:
            return None, None

        scalar = sp.Mul(*scalar_parts) if scalar_parts else sp.S.One
        return scalar, operator_part

    def _compute(self, A: Operator, B: Operator) -> OPEData:
        """
        Compute OPE between two operators.

        Args:
            A: First operator
            B: Second operator

        Returns:
            OPEData result
        """
        # Check for defined OPE
        result = ope_registry.lookup(A, B)
        if result is not None:
            return result

        # Default: OPE is zero (regular)
        return OPEData()

    def __setitem__(
        self, key: tuple[Operator, Operator], value: Union[OPEData, list]
    ) -> None:
        """
        Define an OPE: OPE[A, B] = ...

        Args:
            key: Tuple of (A, B)
            value: OPEData or list of poles
        """
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("Key must be a tuple of two operators")

        A, B = key

        if isinstance(value, list):
            value = OPEData.from_list(value)

        if not isinstance(value, OPEData):
            raise TypeError(f"Value must be OPEData or list, got {type(value)}")

        ope_registry.define(A, B, value)

    @staticmethod
    def make(poles: Union[list, OPEData]) -> OPEData:
        """
        Create OPEData from poles list (Mathematica-style).

        Args:
            poles: List of poles [pole_n, ..., pole_1] or OPEData

        Returns:
            OPEData object
        """
        if isinstance(poles, OPEData):
            return poles
        return OPEData.from_list(poles)


class NOFunction:
    """
    Normal ordered product function.

    Computes NO[A, B] = (AB)(z) using point-splitting convention.
    """

    def __call__(
        self, A: Union[Operator, sp.Expr], B: Union[Operator, sp.Expr]
    ) -> Union[NormalOrderedOperator, sp.Expr]:
        """
        Compute normal ordered product of A and B.

        Args:
            A: First operator or expression
            B: Second operator or expression

        Returns:
            NormalOrderedOperator or simplified expression
        """
        # Handle linearity
        if isinstance(B, sp.Add):
            return sum(self(A, term) for term in B.args)

        if isinstance(A, sp.Add):
            return sum(self(term, B) for term in A.args)

        # Handle scalar multiplication
        if isinstance(B, sp.Mul):
            scalar, op = self._extract_scalar(B)
            if scalar is not None:
                return scalar * self(A, op)

        if isinstance(A, sp.Mul):
            scalar, op = self._extract_scalar(A)
            if scalar is not None:
                return scalar * self(op, B)

        # Handle zero
        if A == 0 or B == 0:
            return sp.S.Zero

        # Now A and B should be operators
        if not isinstance(A, Operator):
            raise TypeError(f"Expected Operator, got {type(A)}")
        if not isinstance(B, Operator):
            raise TypeError(f"Expected Operator, got {type(B)}")

        # Create normal ordered operator
        return NormalOrderedOperator(A, B)

    def _extract_scalar(
        self, expr: sp.Expr
    ) -> tuple[Optional[sp.Expr], Optional[Operator]]:
        """Extract scalar coefficient from expression."""
        if not isinstance(expr, sp.Mul):
            return None, None

        scalar_parts = []
        operator_part = None

        for arg in expr.args:
            if isinstance(arg, Operator):
                if operator_part is not None:
                    return None, None
                operator_part = arg
            else:
                scalar_parts.append(arg)

        if operator_part is None:
            return None, None

        scalar = sp.Mul(*scalar_parts) if scalar_parts else sp.S.One
        return scalar, operator_part


def bracket(
    A: Operator, B: Operator, anticommutator: bool = False
) -> sp.Expr:
    """
    Compute commutator or anticommutator.

    [A, B] = AB - BA (commutator)
    {A, B} = AB + BA (anticommutator)

    Args:
        A: First operator
        B: Second operator
        anticommutator: If True, compute {A,B}, otherwise [A,B]

    Returns:
        Expression for the (anti)commutator
    """
    # This is a placeholder - full implementation requires OPE computation
    # For now, just return the symbolic expression
    AB = NO(A, B)
    BA = NO(B, A)

    if anticommutator:
        return AB + BA
    else:
        return AB - BA


# Create global instances
OPE = OPEFunction()
NO = NOFunction()


__all__ = ["OPE", "NO", "bracket"]
