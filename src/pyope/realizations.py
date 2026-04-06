"""Optional realization backends for accelerated quotient computations."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import sympy as sp

from .api import normal_product
from .constants import Zero
from .operator_spaces import realize_and_simplify
from .operators import DerivativeOperator, NormalOrderedOperator, Operator


class RealizationBackend(ABC):
    """Abstract interface for realization-aware local-operator backends."""

    @abstractmethod
    def realize(self, expr: Any) -> Any:
        raise NotImplementedError

    @abstractmethod
    def sparse_terms(self, expr: Any) -> dict[Any, Any]:
        raise NotImplementedError

    @abstractmethod
    def quotient_normal_form(self, expr: Any) -> Any:
        raise NotImplementedError


class IdentityRealizationBackend(RealizationBackend):
    """Minimal backend that delegates to the canonicalizer without acceleration."""

    def __init__(self, canonicalizer: Any):
        self.canonicalizer = canonicalizer

    def realize(self, expr: Any) -> Any:
        return self.canonicalizer.canonicalize(expr)

    def sparse_terms(self, expr: Any) -> dict[Any, Any]:
        return self.canonicalizer.sparse_terms(expr)

    def quotient_normal_form(self, expr: Any) -> Any:
        return self.canonicalizer.canonicalize(expr)


class DerivativeKillingRealizationBackend(IdentityRealizationBackend):
    """First free-field-style quotient rule: all derivative factors vanish."""

    def realize(self, expr: Any) -> Any:
        return realize_and_simplify(expr)

    def quotient_normal_form(self, expr: Any) -> Any:
        realized = self.realize(expr)
        reduced = self._kill_derivatives(realized)
        return self.canonicalizer.canonicalize(reduced)

    def _kill_derivatives(self, expr: Any) -> Any:
        if expr == Zero:
            return Zero
        if isinstance(expr, DerivativeOperator):
            return Zero
        if isinstance(expr, NormalOrderedOperator):
            left = self._kill_derivatives(expr.left)
            right = self._kill_derivatives(expr.right)
            if left == Zero or right == Zero:
                return Zero
            return normal_product(left, right)
        if isinstance(expr, Operator):
            return expr
        if isinstance(expr, sp.Add):
            return sp.Add(*[self._kill_derivatives(arg) for arg in expr.args])
        if isinstance(expr, sp.Mul):
            coeff, rest = expr.as_coeff_Mul()
            reduced_rest = self._kill_derivatives(rest)
            if reduced_rest == Zero:
                return Zero
            return coeff * reduced_rest
        return expr
