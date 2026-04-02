"""High-level null-search APIs built on sparse $C_2$ reduction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Optional

import sympy as sp

from .c2 import AbstractC2Reducer, C2ReductionWitness
from .constants import Zero
from .operator_spaces import SparseLinearContext, _get_conformal_weight


@dataclass
class NullSearchResult:
    status: str
    target: Any
    quotient_remainder: Any
    obstruction: Any
    witness: C2ReductionWitness
    descendant_basis: list[Any]
    descendant_coeffs: list[sp.Expr]
    null_operator: Any
    c2_remainder: Any
    target_weight: Any
    singular_violations: Any = None
    is_singular: bool | None = None
    n: int | None = None

    def legacy_payload(self, basis_builder: Any) -> dict[str, Any]:
        """Adapt the structured result to the legacy dict payload shape."""
        target_weight = self.target_weight
        if target_weight is None:
            target_vector = sp.zeros(0, 1)
            descendant_matrix = sp.zeros(0, len(self.descendant_basis))
            c2_matrix = sp.zeros(0, len(self.witness.generators))
        else:
            target_vector = basis_builder.coordinates(self.target, target_weight)
            descendant_matrix = _coordinate_matrix_for_basis(
                basis_builder,
                self.descendant_basis,
                target_weight,
            )
            c2_matrix = _coordinate_matrix_for_basis(
                basis_builder,
                self.witness.generators,
                target_weight,
            )

        solution_entries = [*self.descendant_coeffs, *self.witness.coefficients]
        solution_vector = sp.Matrix([[sp.sympify(entry)] for entry in solution_entries])

        payload = {
            "weight": target_weight,
            "target": basis_builder.canonicalize(self.target),
            "target_vector": target_vector,
            "descendant_basis": list(self.descendant_basis),
            "descendant_matrix": descendant_matrix,
            "descendant_coeffs": list(self.descendant_coeffs),
            "c2_basis": list(self.witness.generators),
            "c2_matrix": c2_matrix,
            "c2_coeffs": list(self.witness.coefficients),
            "solution_vector": solution_vector,
            "null_operator": basis_builder.canonicalize(self.null_operator),
            "c2_remainder": basis_builder.canonicalize(self.c2_remainder),
            "status": self.status,
            "is_singular": self.is_singular,
            "singular_violations": self.singular_violations,
            "quotient_remainder": basis_builder.canonicalize(self.quotient_remainder),
            "obstruction": basis_builder.canonicalize(self.obstruction),
        }
        if self.n is not None:
            payload["n"] = self.n
        return payload


class C2NullSearcher:
    """Orchestrate quotient prechecks and descendant lifts for null-state search."""

    def __init__(
        self,
        canonicalizer: Any,
        linear_context: Any = None,
        descendants: Any = None,
        singular_constraints: Any = None,
        c2_reducer: AbstractC2Reducer | None = None,
    ):
        self.canonicalizer = canonicalizer
        self.linear_context = linear_context or SparseLinearContext(canonicalizer)
        self.descendants = descendants
        self.singular_constraints = singular_constraints
        if c2_reducer is None:
            raise ValueError("c2_reducer must be provided")
        self.c2_reducer = c2_reducer

    def quotient_precheck(self, target_expr: Any) -> NullSearchResult:
        target = self.canonicalizer.canonicalize(target_expr)
        witness = self.c2_reducer.solve_c2_witness(target)
        obstruction = witness.remainder
        status = "needs_lift" if obstruction == 0 else "obstructed"
        return NullSearchResult(
            status=status,
            target=target,
            quotient_remainder=obstruction,
            obstruction=obstruction,
            witness=witness,
            descendant_basis=[],
            descendant_coeffs=[],
            null_operator=Zero,
            c2_remainder=witness.c2_part,
            target_weight=witness.target_weight,
        )

    def search_from_sources(
        self, target_weight: Any, sources: Iterable[Any], target_expr: Any
    ) -> Optional[NullSearchResult]:
        if self.descendants is None:
            raise ValueError("descendants must be provided")

        target = self.canonicalizer.canonicalize(target_expr)
        precheck = self.quotient_precheck(target)
        normalized_target = self._quotient_normal_form(target)

        descendant_basis = list(self.descendants.span(sources, target_weight))
        quotient_descendants = [
            self._quotient_normal_form(expr) for expr in descendant_basis
        ]

        solution = self._solve_descendant_lift(quotient_descendants, normalized_target)
        if solution is None:
            if precheck.status == "obstructed":
                return precheck
            return None

        descendant_coeffs = solution
        null_operator = self.canonicalizer.canonicalize(
            self._linear_combination(descendant_basis, descendant_coeffs)
        )
        witness = self.c2_reducer.solve_c2_witness(null_operator - target)
        c2_remainder = self.canonicalizer.canonicalize(target - null_operator)
        singular_violations = None
        is_singular = None
        status = "solved"

        if self.singular_constraints is not None:
            if hasattr(self.singular_constraints, "positive_mode_constraints"):
                singular_violations = (
                    self.singular_constraints.positive_mode_constraints(null_operator)
                )
                is_singular = all(
                    not violations for violations in singular_violations.values()
                )
            elif hasattr(self.singular_constraints, "is_singular"):
                is_singular = bool(self.singular_constraints.is_singular(null_operator))
                singular_violations = {} if is_singular else {"unknown": True}

            if is_singular is False:
                status = "failed_singularity"

        return NullSearchResult(
            status=status,
            target=target,
            quotient_remainder=precheck.quotient_remainder,
            obstruction=precheck.obstruction,
            witness=witness,
            descendant_basis=descendant_basis,
            descendant_coeffs=descendant_coeffs,
            null_operator=null_operator,
            c2_remainder=c2_remainder,
            target_weight=target_weight,
            singular_violations=singular_violations,
            is_singular=is_singular,
        )

    def search_stress_tensor_nilpotency(
        self, n: int, sources: Any = None
    ) -> NullSearchResult | None:
        target = self.canonicalizer.nested_stress_tensor(n)
        if sources is None:
            result = self.quotient_precheck(target)
            result.n = n
            return result
        target_weight = _get_conformal_weight(target)
        if target_weight is None:
            raise ValueError("Could not infer target weight for stress tensor product")
        result = self.search_from_sources(target_weight, sources, target)
        if result is not None:
            result.n = n
            return result
        fallback = self.quotient_precheck(target)
        fallback.n = n
        return fallback

    def _quotient_normal_form(self, expr: Any) -> Any:
        return self.canonicalizer.canonicalize(
            self.c2_reducer.quotient_normal_form(expr)
        )

    def _solve_descendant_lift(
        self, quotient_descendants: Iterable[Any], target: Any
    ) -> list[sp.Expr] | None:
        target_terms = self.linear_context.sparse_terms(target)
        descendant_terms = [
            self.linear_context.sparse_terms(expr) for expr in quotient_descendants
        ]

        monomials = sorted(
            {monomial for terms in descendant_terms for monomial in terms}
            | set(target_terms),
            key=sp.srepr,
        )
        if not monomials:
            return [] if target == Zero else None

        columns = []
        for terms in descendant_terms:
            columns.append(
                sp.Matrix(
                    [sp.sympify(terms.get(monomial, 0)) for monomial in monomials]
                )
            )
        rhs = sp.Matrix(
            [sp.sympify(target_terms.get(monomial, 0)) for monomial in monomials]
        )

        if not columns:
            return None

        matrix = sp.Matrix.hstack(*columns)
        if matrix.rank() != sp.Matrix.hstack(matrix, rhs).rank():
            return None

        solution = matrix.gauss_jordan_solve(rhs)
        if isinstance(solution, tuple):
            vector = solution[0]
            parameters = solution[1]
            if parameters.shape[0] > 0:
                vector = vector.xreplace({param: 0 for param in parameters})
        else:
            vector = solution

        return [sp.sympify(vector[i, 0]) for i in range(vector.rows)]

    def _linear_combination(self, basis: Iterable[Any], coeffs: Iterable[Any]) -> Any:
        terms = []
        for coeff, expr in zip(coeffs, basis):
            coeff = sp.sympify(coeff)
            if coeff != 0:
                terms.append(coeff * expr)
        if not terms:
            return Zero
        return sp.Add(*terms)


def _coordinate_matrix_for_basis(
    basis_builder: Any, expressions: Iterable[Any], weight: Any
) -> sp.Matrix:
    columns = [basis_builder.coordinates(expr, weight=weight) for expr in expressions]
    if not columns:
        basis_dim = len(basis_builder.list(weight))
        return sp.zeros(basis_dim, 0)
    return sp.Matrix.hstack(*columns)
