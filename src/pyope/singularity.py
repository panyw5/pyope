"""Singular-vector and positive-mode constraint utilities."""

from __future__ import annotations

from typing import Any, Iterable, Optional

import sympy as sp

from .api import OPE
from .constants import Zero
from .local_operator import collect_operator_terms
from .operator_spaces import (
    LocalOperatorBasis,
    _get_conformal_weight,
    _normalize_weight,
    _weights_equal,
)


class SingularVectorAnalyzer:
    """Analyze singular / primary constraints from OPE pole data."""

    def __init__(
        self,
        basis_builder: LocalOperatorBasis,
        generators: Optional[Iterable[Any]] = None,
        stress_tensor: Any = None,
    ):
        self.basis_builder = basis_builder
        self.generators = (
            tuple(generators)
            if generators is not None
            else tuple(basis_builder.generators)
        )
        self.stress_tensor = stress_tensor

    def positive_mode_constraints(
        self, expr: Any, generators: Optional[Iterable[Any]] = None
    ) -> dict[Any, dict[int, Any]]:
        """Return forbidden positive-mode pole coefficients for each constraint generator."""
        constraint_generators = (
            tuple(generators) if generators is not None else self.generators
        )
        canonical_expr = self.basis_builder.canonicalize(expr)
        constraints: dict[Any, dict[int, Any]] = {}

        for generator in constraint_generators:
            ope = OPE(generator, canonical_expr)
            allowed_max_pole = self._allowed_max_pole(generator)
            violating: dict[int, Any] = {}
            for pole in range(allowed_max_pole + 1, ope.max_pole + 1):
                coeff = self.basis_builder.canonicalize(ope.pole(pole))
                if coeff != Zero:
                    violating[pole] = coeff
            constraints[generator] = violating

        return constraints

    def is_singular(
        self, expr: Any, generators: Optional[Iterable[Any]] = None
    ) -> bool:
        """Check whether all forbidden positive-mode poles vanish."""
        constraints = self.positive_mode_constraints(expr, generators=generators)
        return all(not violations for violations in constraints.values())

    def find_singular_vectors(
        self,
        weight: Any,
        ansatz: Optional[Iterable[Any]] = None,
        generators: Optional[Iterable[Any]] = None,
    ) -> list[dict[str, Any]]:
        """Solve linear ansaetze for singular vectors at fixed weight."""
        normalized_weight = _normalize_weight(weight)
        ansatz_basis = (
            list(ansatz) if ansatz is not None else self.basis_builder.list(weight)
        )
        ansatz_basis = [
            self.basis_builder.canonicalize(expr)
            for expr in ansatz_basis
            if _weights_equal(_get_conformal_weight(expr), normalized_weight)
        ]
        if not ansatz_basis:
            return []

        constraint_generators = (
            tuple(generators) if generators is not None else self.generators
        )
        equations = []
        symbols = sp.symbols(f"c0:{len(ansatz_basis)}")
        trial_expr = sp.Add(
            *[symbol * expr for symbol, expr in zip(symbols, ansatz_basis)]
        )

        for generator in constraint_generators:
            ope = OPE(generator, trial_expr)
            allowed_max_pole = self._allowed_max_pole(generator)
            for pole in range(allowed_max_pole + 1, ope.max_pole + 1):
                coeff = self.basis_builder.canonicalize(ope.pole(pole))
                if coeff == Zero:
                    continue

                coeff_terms = collect_operator_terms(coeff)
                for operator_coeff in coeff_terms.values():
                    equations.append(sp.sympify(operator_coeff))

        if not equations:
            results = []
            for index, basis_expr in enumerate(ansatz_basis):
                coeffs = {
                    symbol: (sp.Integer(1) if i == index else sp.Integer(0))
                    for i, symbol in enumerate(symbols)
                }
                results.append({"vector": basis_expr, "coefficients": coeffs})
            return results

        matrix, _ = sp.linear_eq_to_matrix(equations, symbols)
        nullspace = matrix.nullspace()
        results = []
        for basis_vector in nullspace:
            expr = self.basis_builder.canonicalize(
                sp.Add(
                    *[
                        coeff * basis_expr
                        for coeff, basis_expr in zip(basis_vector, ansatz_basis)
                        if coeff != 0
                    ]
                )
            )
            if expr != Zero:
                results.append(
                    {
                        "vector": expr,
                        "coefficients": {
                            symbol: coeff
                            for symbol, coeff in zip(symbols, basis_vector)
                        },
                    }
                )

        return results

    def _allowed_max_pole(self, generator: Any) -> int:
        if self.stress_tensor is not None and generator == self.stress_tensor:
            return 2
        return 0
