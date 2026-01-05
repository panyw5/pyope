"""
Special constants and operators.

This module defines special operators like One (unit operator) and Delta (Kronecker delta).
"""

from __future__ import annotations
from typing import Any
import sympy as sp
from .operators import Operator


class ConstantOperator(Operator):
    """
    Constant operator (like the unit operator One).

    Constant operators have zero derivative and regular OPE with all operators.
    """

    def __new__(cls, name: str, **kwargs):
        """
        Create a constant operator.

        Args:
            name: Name of the constant operator
            **kwargs: Additional arguments
        """
        # Create operator with parity=0 (always bosonic)
        return Operator.__new__(cls, name, parity=0, **kwargs)


# The unit operator (identity)
One = ConstantOperator("One")


class DeltaFunction(sp.Function):
    """
    Kronecker delta function: δ_{ij}.

    Returns 1 if i == j, 0 otherwise.
    For symbolic indices, remains unevaluated.
    """

    @classmethod
    def eval(cls, i: Any, j: Any) -> sp.Expr:
        """
        Evaluate the delta function.

        Args:
            i: First index
            j: Second index

        Returns:
            1 if i == j, 0 if i != j and both are numbers, None otherwise
        """
        # If indices are equal
        if i == j:
            return sp.S.One

        # If both are numbers and different
        if isinstance(i, (int, sp.Integer)) and isinstance(j, (int, sp.Integer)):
            return sp.S.Zero

        # Otherwise, keep symbolic
        return None

    def _eval_simplify(self, **kwargs: Any) -> sp.Expr:
        """Simplify the delta function."""
        i, j = self.args
        if i == j:
            return sp.S.One
        return self


# Create a convenient Delta function
def Delta(i: Any, j: Any) -> sp.Expr:
    """
    Kronecker delta function.

    Args:
        i: First index
        j: Second index

    Returns:
        1 if i == j, 0 if different numbers, Delta(i,j) if symbolic
    """
    return DeltaFunction(i, j)
