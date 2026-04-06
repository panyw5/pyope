"""Operator-space utilities for fixed-weight local-operator computations."""

from __future__ import annotations

from collections import deque
import inspect
from functools import lru_cache
from typing import Any, Iterable, MutableMapping, Optional

import sympy as sp

from .api import OPE, normal_product
from .constants import One, Zero
from .local_operator import (
    collect_operator_terms,
    extract_scalar_operator,
    get_operator_parity,
)
from .operators import BasisOperator, Operator, d
from .registry import ope_registry
from .simplify import simplify


class RealizedGenerator(BasisOperator):
    """Named strong generator with an optional composite realization."""

    def __new__(
        cls,
        name: str,
        realization: Any,
        conformal_weight: Optional[Any] = None,
        fermionic: Optional[bool] = None,
        **assumptions,
    ):
        if conformal_weight is None:
            conformal_weight = _get_conformal_weight(realization)
        if conformal_weight is None:
            raise ValueError(
                "RealizedGenerator requires a well-defined conformal weight"
            )

        if fermionic is None:
            fermionic = bool(get_operator_parity(realization))

        obj = BasisOperator.__new__(
            cls,
            name,
            fermionic=fermionic,
            conformal_weight=conformal_weight,
            **assumptions,
        )
        obj._realization = realization
        return obj

    @property
    def realization(self) -> Any:
        return getattr(self, "_realization")


def _resolve_calling_namespace() -> MutableMapping[str, Any]:
    frame = inspect.currentframe()
    if frame is None or frame.f_back is None or frame.f_back.f_back is None:
        raise RuntimeError("Could not access calling frame for make_realized")
    caller = frame.f_back.f_back
    return (
        caller.f_locals if caller.f_locals is not caller.f_globals else caller.f_globals
    )


def _lookup_binding_name(namespace: MutableMapping[str, Any], expr: Any) -> str:
    matches = [name for name, value in namespace.items() if value is expr]
    if not matches:
        raise ValueError(
            "Could not infer variable name for expression in make_realized"
        )
    return matches[-1]


def make_realized(expressions: Any, **assumptions) -> list[RealizedGenerator]:
    """Convert bound expressions to same-named `RealizedGenerator` objects.

    Typical usage is to first define composite expressions naturally and then
    promote them by passing the variables themselves:

        J0 = NO(b, c) + 2 * NO(beta, gamma)
        Jminus = -NO(beta, NO(gamma, gamma)) + sp.Rational(3, 2) * d(gamma)
        J0, Jminus = make_realized([J0, Jminus])

    Args:
        expressions: List/tuple/set of expressions already bound to variables
            in the caller namespace.
        **assumptions: Extra keyword arguments forwarded to `RealizedGenerator`.

    Returns:
        Promoted `RealizedGenerator` objects in the same order as the input.
    """
    if not isinstance(expressions, (list, tuple, set)):
        raise TypeError("make_realized expects a list, tuple, or set of expressions")

    exprs = list(expressions)

    if not exprs:
        raise ValueError("make_realized requires at least one expression")

    namespace = _resolve_calling_namespace()
    promoted: list[RealizedGenerator] = []

    for expr in exprs:
        name = _lookup_binding_name(namespace, expr)
        realization = expr.realization if isinstance(expr, RealizedGenerator) else expr
        generator = RealizedGenerator(name, realization=realization, **assumptions)
        promoted.append(generator)

    return promoted


def _normalize_weight(value: Any) -> sp.Expr:
    """Normalize a conformal weight into a stable SymPy scalar."""
    return sp.nsimplify(value)


def _normalize_max_occurence(value: Any) -> int | None:
    """Normalize an occurrence cutoff for nonpositive bosonic atoms."""
    if value is None:
        return None

    normalized = sp.nsimplify(value)
    if normalized.is_integer is False:
        raise ValueError("max_occurence must be a nonnegative integer")

    try:
        max_occurence = int(normalized)
    except TypeError as exc:
        raise ValueError("max_occurence must be a nonnegative integer") from exc

    if max_occurence < 0:
        raise ValueError("max_occurence must be a nonnegative integer")

    return max_occurence


def _independent_from_vectors(
    expressions: Iterable[Any], vector_getter: Any
) -> list[Any]:
    ordered = sorted(set(expressions), key=sp.srepr)
    independent = []

    if ordered:
        sample_vector = vector_getter(ordered[0])
        if isinstance(sample_vector, dict):
            eliminator = _SparseIndependentEliminator()
            if sample_vector:
                if eliminator.insert(sample_vector):
                    independent.append(ordered[0])

            for expr in ordered[1:]:
                vector = vector_getter(expr)
                if eliminator.insert(vector):
                    independent.append(expr)
            return independent

    columns = []

    for expr in ordered:
        vector = vector_getter(expr)
        is_zero_vector = (
            not vector
            if isinstance(vector, dict)
            else vector == sp.zeros(*vector.shape)
        )
        if not columns:
            if not is_zero_vector:
                independent.append(expr)
                columns.append(vector)
            continue

        if isinstance(vector, dict):
            current, _ = _matrix_from_sparse_vectors(columns)
            candidate, _ = _matrix_from_sparse_vectors([*columns, vector])
        else:
            current = sp.Matrix.hstack(*columns)
            candidate = sp.Matrix.hstack(*columns, vector)
        if candidate.rank() > current.rank():
            independent.append(expr)
            columns.append(vector)

    return independent


class _SparseIndependentEliminator:
    """Incremental sparse elimination for independence checks."""

    def __init__(self) -> None:
        self._pivot_vectors: dict[Any, dict[Any, sp.Expr]] = {}
        self._pivot_order: list[Any] = []

    def _pivot_key(self, monomial: Any) -> str:
        return sp.srepr(monomial)

    def _leading_term(self, vector: dict[Any, sp.Expr]) -> tuple[Any, sp.Expr] | None:
        if not vector:
            return None
        monomial = max(vector, key=self._pivot_key)
        return monomial, vector[monomial]

    def reduce(self, vector: dict[Any, sp.Expr]) -> dict[Any, sp.Expr]:
        reduced, _ = self.reduce_with_coefficients(vector)
        return reduced

    def reduce_with_coefficients(
        self, vector: dict[Any, sp.Expr]
    ) -> tuple[dict[Any, sp.Expr], dict[Any, sp.Expr]]:
        reduced = {
            monomial: sp.sympify(coeff)
            for monomial, coeff in vector.items()
            if coeff != 0
        }
        coefficients: dict[Any, sp.Expr] = {}

        while reduced:
            leading = self._leading_term(reduced)
            if leading is None:
                break

            pivot_monomial, pivot_coeff = leading
            pivot_vector = self._pivot_vectors.get(pivot_monomial)
            if pivot_vector is None:
                break

            factor = sp.sympify(pivot_coeff / pivot_vector[pivot_monomial])
            coefficients[pivot_monomial] = sp.sympify(
                coefficients.get(pivot_monomial, 0) + factor
            )
            for monomial, coeff in pivot_vector.items():
                updated = sp.sympify(reduced.get(monomial, 0) - factor * coeff)
                if updated == 0:
                    reduced.pop(monomial, None)
                else:
                    reduced[monomial] = updated

        coefficients = {
            monomial: coeff for monomial, coeff in coefficients.items() if coeff != 0
        }
        return reduced, coefficients

    def insert(self, vector: dict[Any, sp.Expr]) -> bool:
        reduced = self.reduce(vector)
        return self.insert_reduced(reduced)

    def insert_reduced(self, reduced: dict[Any, sp.Expr]) -> bool:
        leading = self._leading_term(reduced)
        if leading is None:
            return False

        pivot_monomial, pivot_coeff = leading
        normalized = {
            monomial: sp.sympify(coeff / pivot_coeff)
            for monomial, coeff in reduced.items()
        }
        self._pivot_vectors[pivot_monomial] = normalized
        self._pivot_order.append(pivot_monomial)
        return True

    @property
    def pivot_order(self) -> list[Any]:
        return list(self._pivot_order)


def _canonical_terms(expr: Any, canonicalizer: Any) -> dict[Any, sp.Expr]:
    """Return a sparse canonical expansion of an operator expression."""
    canonical = canonicalizer(expr)
    if canonical == Zero:
        return {}

    terms = {}
    for operator, coeff in collect_operator_terms(canonical).items():
        coeff = sp.sympify(coeff)
        if coeff == 0:
            continue

        monomial = One if operator == 1 else operator
        terms[monomial] = sp.sympify(terms.get(monomial, 0) + coeff)

    return {monomial: coeff for monomial, coeff in terms.items() if coeff != 0}


def _matrix_from_sparse_vectors(
    vectors: Iterable[dict[Any, sp.Expr]],
) -> tuple[sp.Matrix, list[Any]]:
    """Build a dense matrix from sparse vectors over canonical monomials."""
    vector_list = list(vectors)
    if not vector_list:
        return sp.zeros(0, 0), []

    monomials = sorted(
        {monomial for vector in vector_list for monomial in vector}, key=sp.srepr
    )
    if not monomials:
        return sp.zeros(0, len(vector_list)), []

    row_index = {monomial: idx for idx, monomial in enumerate(monomials)}
    matrix = sp.zeros(len(monomials), len(vector_list))

    for col, vector in enumerate(vector_list):
        for monomial, coeff in vector.items():
            matrix[row_index[monomial], col] = coeff

    return matrix, monomials


def _compressed_matrix_from_sparse_vectors(
    vectors: Iterable[dict[Any, sp.Expr]],
) -> tuple[sp.Matrix, list[Any]]:
    """Build a compressed coordinate matrix using incremental sparse pivots."""
    vector_list = list(vectors)
    if not vector_list:
        return sp.zeros(0, 0), []

    eliminator = _SparseIndependentEliminator()
    coordinate_columns: list[dict[Any, sp.Expr]] = []

    for vector in vector_list:
        reduced, coefficients = eliminator.reduce_with_coefficients(vector)
        leading = eliminator._leading_term(reduced)

        column = dict(coefficients)
        if leading is not None:
            pivot_monomial, pivot_coeff = leading
            eliminator.insert_reduced(reduced)
            column[pivot_monomial] = sp.sympify(
                column.get(pivot_monomial, 0) + pivot_coeff
            )

        coordinate_columns.append(
            {key: value for key, value in column.items() if value != 0}
        )

    pivots = eliminator.pivot_order
    if not pivots:
        return sp.zeros(0, len(vector_list)), []

    row_index = {pivot: idx for idx, pivot in enumerate(pivots)}
    matrix = sp.zeros(len(pivots), len(vector_list))
    for col, vector in enumerate(coordinate_columns):
        for pivot, coeff in vector.items():
            matrix[row_index[pivot], col] = coeff

    return matrix, pivots


def _zero_relations_from_vectors(
    expressions: Iterable[Any], vector_getter: Any
) -> list[dict[str, Any]]:
    ordered = list(expressions)
    if not ordered:
        return []

    columns = [vector_getter(expr) for expr in ordered]
    if not columns:
        return []

    first = columns[0]
    if isinstance(first, dict):
        matrix, _ = _compressed_matrix_from_sparse_vectors(columns)
    else:
        matrix = sp.Matrix.hstack(*columns)
    relations = []

    for basis_vector in matrix.nullspace():
        coefficients = [
            sp.sympify(basis_vector[index, 0]) for index in range(len(ordered))
        ]
        terms = [
            (expr, coeff) for expr, coeff in zip(ordered, coefficients) if coeff != 0
        ]
        relation = (
            sp.Add(*[coeff * expr for expr, coeff in terms], evaluate=False)
            if terms
            else Zero
        )
        relations.append(
            {
                "operators": ordered,
                "coefficients": coefficients,
                "coefficient_vector": basis_vector,
                "terms": terms,
                "relation": relation,
            }
        )

    return relations


def _validate_expression_weight(expr: Any, weight: Any = None) -> None:
    """Validate an expression against an optional requested conformal weight."""
    if weight is None:
        return

    expr_weight = _get_conformal_weight(expr)
    if expr_weight is not None and not _weights_equal(
        expr_weight, _normalize_weight(weight)
    ):
        raise ValueError("Expression weight does not match requested weight")


def _weights_equal(left: Any, right: Any) -> bool:
    return sp.simplify(_normalize_weight(left) - _normalize_weight(right)) == 0


def _get_conformal_weight(expr: Any) -> sp.Expr | None:
    """Return the conformal weight of an operator expression, if defined."""
    if expr is Zero:
        return None
    if expr is One:
        return sp.Integer(0)

    direct_weight = getattr(expr, "conformal_weight", None)
    if direct_weight is not None:
        return _normalize_weight(direct_weight)

    if isinstance(expr, Operator):
        weight = getattr(expr, "conformal_weight", None)
        return None if weight is None else _normalize_weight(weight)

    if isinstance(expr, sp.Add):
        weights = []
        for term in expr.args:
            term_weight = _get_conformal_weight(term)
            if term_weight is None:
                return None
            weights.append(term_weight)

        if not weights:
            return None

        first = weights[0]
        if any(not _weights_equal(first, weight) for weight in weights[1:]):
            raise ValueError(
                "Expression contains terms with inconsistent conformal weights"
            )
        return first

    if isinstance(expr, sp.Mul):
        _, operator = extract_scalar_operator(expr)
        return _get_conformal_weight(operator)

    if isinstance(expr, sp.Expr) and not expr.has(Operator):
        return sp.Integer(0)

    weight = getattr(expr, "conformal_weight", None)
    return None if weight is None else _normalize_weight(weight)


def _operator_sort_key(operator: Operator) -> tuple:
    """Provide a stable ordering consistent with compare_operators()."""
    base = getattr(operator, "base", operator)
    order = int(getattr(operator, "order", 0))

    position = ope_registry.get_position(base)
    if position is not None:
        base_key = (0, position)
    else:
        base_name = getattr(base, "name", str(base))
        base_key = (1, base_name)

    name = getattr(operator, "name", str(operator))
    return (base_key, -order, name)


def _underlying_basis_generators(expr: Any) -> set[Operator]:
    """Collect underlying basis generators appearing in an operator expression."""
    if not isinstance(expr, Operator):
        return set()

    base = getattr(expr, "base", None)
    if base is not None:
        return _underlying_basis_generators(base)

    left = getattr(expr, "left", None)
    right = getattr(expr, "right", None)
    if left is not None and right is not None:
        return _underlying_basis_generators(left) | _underlying_basis_generators(right)

    return {expr}


def _is_negative_single_letter_fermion_sector(expr: Any) -> bool:
    """Return whether an expression is built from one negative-weight fermionic letter."""
    generators = _underlying_basis_generators(expr)
    if len(generators) != 1:
        return False

    generator = next(iter(generators))
    generator_weight = _get_conformal_weight(generator)
    if generator_weight is None or generator_weight >= 0:
        return False

    return bool(get_operator_parity(generator))


def _combine_like_terms_preserving_metadata(expr: Any) -> Any:
    """Combine like terms without rebuilding operators from stripped keys."""
    terms = collect_operator_terms(expr)
    rebuilt_terms = []

    for operator, coeff in terms.items():
        coeff = sp.sympify(coeff)
        if coeff == 0:
            continue

        if operator == 1 or operator is One:
            rebuilt_terms.append(coeff)
        elif coeff == 1:
            rebuilt_terms.append(operator)
        else:
            rebuilt_terms.append(coeff * operator)

    if not rebuilt_terms:
        return Zero
    if len(rebuilt_terms) == 1:
        return rebuilt_terms[0]
    return sp.Add(*rebuilt_terms)


def _realize_expr(expr: Any) -> Any:
    """Recursively expand realized generators into their underlying expressions."""
    if isinstance(expr, RealizedGenerator):
        return _realize_expr(expr.realization)

    if expr is Zero or expr is One:
        return expr

    if isinstance(expr, sp.Add):
        return sp.Add(*[_realize_expr(arg) for arg in expr.args])

    if isinstance(expr, sp.Mul):
        coeff, operator = extract_scalar_operator(expr)
        realized_operator = _realize_expr(operator)
        if realized_operator == 1:
            return sp.sympify(coeff)
        if coeff == 1:
            return realized_operator
        return sp.sympify(coeff) * realized_operator

    if isinstance(expr, Operator):
        base = getattr(expr, "base", None)
        order = getattr(expr, "order", None)
        left = getattr(expr, "left", None)
        right = getattr(expr, "right", None)

        if base is not None and order is not None:
            return d(_realize_expr(base), order)
        if left is not None and right is not None:
            return normal_product(_realize_expr(left), _realize_expr(right))
        return expr

    return expr


def realize(expr: Any) -> Any:
    """Expand realized generators in an operator expression."""
    realize_method = getattr(expr, "realize", None)
    if callable(realize_method):
        return realize_method()
    return _realize_expr(expr)


def realize_and_simplify(expr: Any) -> Any:
    """Expand realized generators and canonicalize the resulting expression."""
    return _combine_like_terms_preserving_metadata(simplify(expr.realize()))


def realized_coordinates(
    expr: Any,
    free_field_basis: "LocalOperatorBasis",
    weight: Any = None,
    max_occurence: Any = None,
) -> sp.Matrix:
    """Project an expression to coordinates after realization expansion."""
    return free_field_basis.coordinates(
        realize_and_simplify(expr),
        weight=weight,
        max_occurence=max_occurence,
    )


def list_independent_ops(
    expressions: Iterable[Any],
    basis_builder: "LocalOperatorBasis",
    weight: Any = None,
    max_occurence: Any = None,
) -> list[Any]:
    """Return a linearly independent subset via local canonical expansions."""

    def vector_getter(expr: Any) -> dict[Any, sp.Expr]:
        _validate_expression_weight(expr, weight)
        return _canonical_terms(expr, basis_builder.canonicalize)

    return _independent_from_vectors(expressions, vector_getter)


def list_zero_relations(
    expressions: Iterable[Any],
    basis_builder: "LocalOperatorBasis",
    weight: Any = None,
    max_occurence: Any = None,
) -> list[dict[str, Any]]:
    """Return a basis of linear relations via local canonical expansions."""

    def vector_getter(expr: Any) -> dict[Any, sp.Expr]:
        _validate_expression_weight(expr, weight)
        return _canonical_terms(expr, basis_builder.canonicalize)

    return _zero_relations_from_vectors(expressions, vector_getter)


def independent_under_realization(
    expressions: Iterable[Any],
    free_field_basis: "LocalOperatorBasis",
    weight: Any = None,
    max_occurence: Any = None,
) -> list[Any]:
    """Return expressions that stay linearly independent after realization."""

    def vector_getter(expr: Any) -> sp.Matrix:
        target_weight = weight
        if target_weight is None:
            target_weight = _get_conformal_weight(expr)
            if target_weight is None:
                raise ValueError(
                    "Could not infer conformal weight under realization; provide weight"
                )

        return realized_coordinates(
            expr,
            free_field_basis,
            weight=target_weight,
            max_occurence=max_occurence,
        )

    return _independent_from_vectors(expressions, vector_getter)


class LocalOperatorBasis:
    """Build canonical fixed-weight bases of local operators."""

    def __init__(
        self,
        generators: Iterable[Operator],
        max_weight: Any = None,
        max_occurence: Any = None,
    ):
        self.generators = tuple(generators)
        self.max_weight = None if max_weight is None else _normalize_weight(max_weight)
        self.max_occurence = _normalize_max_occurence(max_occurence)

        if not self.generators:
            raise ValueError("LocalOperatorBasis requires at least one generator")

        has_nonpositive_generator = False
        for generator in self.generators:
            if not isinstance(generator, Operator):
                raise TypeError(
                    "LocalOperatorBasis generators must be Operator instances"
                )
            generator_weight = _get_conformal_weight(generator)
            if generator_weight is None:
                raise ValueError(
                    f"Generator {generator} must define a conformal weight"
                )
            if generator_weight <= 0:
                has_nonpositive_generator = True

        if has_nonpositive_generator and self.max_occurence is None:
            raise ValueError(
                "LocalOperatorBasis requires max_occurence when generators include "
                "nonpositive conformal weights"
            )

    def _resolve_max_occurence(self, value: Any = None) -> int | None:
        normalized = _normalize_max_occurence(value)
        if normalized is None:
            return self.max_occurence
        return normalized

    def canonicalize(self, expr: Any) -> Any:
        """Canonicalize an operator expression for basis comparisons."""
        simplified = self._canonicalize_cached(expr, ope_registry.version)
        combined = _combine_like_terms_preserving_metadata(simplified)
        if combined == 0:
            return Zero
        return combined

    @staticmethod
    @lru_cache(maxsize=None)
    def _canonicalize_cached(expr: Any, registry_version: int) -> Any:
        del registry_version
        return simplify(expr)

    def enumerate_candidates(self, weight: Any) -> list[Any]:
        """Enumerate raw canonical candidate expressions of fixed weight."""
        normalized_weight = _normalize_weight(weight)
        if self.max_weight is not None and normalized_weight > self.max_weight:
            raise ValueError("Requested weight exceeds configured max_weight")

        basis_terms = self._basis_terms(
            normalized_weight, self._resolve_max_occurence()
        )
        return list(basis_terms)

    def basis(self, weight: Any, max_occurence: Any = None) -> list[Any]:
        """Return a canonical basis of operator monomials at fixed weight."""
        normalized_weight = _normalize_weight(weight)
        normalized_max_occurence = self._resolve_max_occurence(max_occurence)
        return list(self._basis_terms(normalized_weight, normalized_max_occurence))

    def coordinates(
        self, expr: Any, weight: Any = None, max_occurence: Any = None
    ) -> sp.Matrix:
        """Return basis coordinates of an operator expression as a column vector."""
        canonical = self.canonicalize(expr)
        normalized_max_occurence = self._resolve_max_occurence(max_occurence)

        if canonical == Zero:
            if weight is None:
                inferred_weight = _get_conformal_weight(expr)
                if inferred_weight is None:
                    raise ValueError("Weight is required for the zero expression")
                weight = inferred_weight
            basis = self.basis(weight, max_occurence=normalized_max_occurence)
            return sp.zeros(len(basis), 1)

        expr_weight = _get_conformal_weight(canonical)
        if expr_weight is None and weight is None:
            raise ValueError(
                "Could not infer conformal weight; please provide weight explicitly"
            )

        target_weight = expr_weight if weight is None else _normalize_weight(weight)
        if expr_weight is not None and not _weights_equal(expr_weight, target_weight):
            raise ValueError("Expression weight does not match requested basis weight")

        basis = self.basis(target_weight, max_occurence=normalized_max_occurence)
        index = {op: i for i, op in enumerate(basis)}
        vector = sp.zeros(len(basis), 1)

        for operator, coeff in collect_operator_terms(canonical).items():
            if operator == 1:
                operator = One
            if operator not in index:
                raise ValueError(
                    f"Operator {operator} is not in the canonical basis at weight {target_weight}"
                )
            vector[index[operator], 0] += sp.sympify(coeff)

        return vector

    def realized_coordinates(
        self,
        expr: Any,
        free_field_basis: "LocalOperatorBasis",
        weight: Any = None,
        max_occurence: Any = None,
    ) -> sp.Matrix:
        """Project an expression after expanding realized generators."""
        return realized_coordinates(
            expr,
            free_field_basis,
            weight=weight,
            max_occurence=max_occurence,
        )

    def independent_under_realization(
        self,
        expressions: Iterable[Any],
        free_field_basis: "LocalOperatorBasis",
        weight: Any = None,
        max_occurence: Any = None,
    ) -> list[Any]:
        """Filter expressions by linear independence after realization."""
        return independent_under_realization(
            expressions,
            free_field_basis,
            weight=weight,
            max_occurence=max_occurence,
        )

    def list_independent_ops(
        self,
        expressions: Iterable[Any],
        weight: Any = None,
        max_occurence: Any = None,
    ) -> list[Any]:
        """Filter expressions by linear independence in this basis."""
        return list_independent_ops(
            expressions,
            self,
            weight=weight,
            max_occurence=max_occurence,
        )

    def list_zero_relations(
        self,
        expressions: Iterable[Any],
        weight: Any = None,
        max_occurence: Any = None,
    ) -> list[dict[str, Any]]:
        """Return a basis of linear relations among expressions in this basis."""
        return list_zero_relations(
            expressions,
            self,
            weight=weight,
            max_occurence=max_occurence,
        )

    @lru_cache(maxsize=None)
    def _basis_terms(
        self, weight: sp.Expr, max_occurence: int | None = None
    ) -> tuple[Any, ...]:
        atomic_operators = self._atomic_operators(weight, max_occurence)
        canonical_terms: set[Any] = set()
        nonpositive_bosonic_counts = tuple(0 for _ in atomic_operators)

        for factors in self._factor_tuples(
            weight,
            0,
            atomic_operators,
            max_occurence=max_occurence,
            nonpositive_bosonic_counts=nonpositive_bosonic_counts,
        ):
            if len(factors) == 0:
                canonical_terms.add(One)
                continue
            candidate = factors[0] if len(factors) == 1 else normal_product(*factors)
            canonical = self.canonicalize(candidate)
            for operator, coeff in collect_operator_terms(canonical).items():
                if coeff == 0:
                    continue
                operator_weight = _get_conformal_weight(operator)
                if operator_weight is None:
                    continue
                if _weights_equal(operator_weight, weight):
                    if weight == 0 and _is_negative_single_letter_fermion_sector(
                        operator
                    ):
                        continue
                    canonical_terms.add(operator)

        return tuple(sorted(canonical_terms, key=sp.srepr))

    def _atomic_operators(
        self, target_weight: sp.Expr, max_occurence: int | None = None
    ) -> list[tuple[Operator, sp.Expr, bool]]:
        atoms = []
        min_fermionic_compensation = self._min_fermionic_compensation()
        max_atomic_weight = sp.simplify(target_weight - min_fermionic_compensation)

        for generator in self.generators:
            generator_weight = _get_conformal_weight(generator)
            if generator_weight is None:
                raise ValueError(
                    f"Generator {generator} must define a conformal weight"
                )
            max_order_expr = sp.floor(max_atomic_weight - generator_weight)
            if max_order_expr.is_integer is False:
                continue

            max_order = int(max_order_expr)
            if max_order < 0:
                continue

            for order in range(max_order + 1):
                operator = generator if order == 0 else d(generator, order)
                operator_weight = _get_conformal_weight(operator)
                if operator_weight is None:
                    continue
                fermionic = bool(get_operator_parity(operator))
                if operator_weight <= 0 and not fermionic and max_occurence is None:
                    raise ValueError(
                        "LocalOperatorBasis cannot unrestrictedly enumerate fixed-weight spaces "
                        "with nonpositive bosonic atomic operators; provide max_occurence "
                        "to truncate nonpositive bosonic atoms; got "
                        f"{operator} with weight {operator_weight}"
                    )
                if operator_weight <= max_atomic_weight:
                    atoms.append((operator, operator_weight, fermionic))

        atoms.sort(
            key=lambda item: (_operator_sort_key(item[0]), item[1], sp.srepr(item[0]))
        )
        return atoms

    def _min_fermionic_compensation(self) -> sp.Expr:
        compensation = sp.Integer(0)

        for generator in self.generators:
            generator_weight = _get_conformal_weight(generator)
            if generator_weight is None:
                continue

            fermionic = bool(get_operator_parity(generator))
            if fermionic:
                max_nonpositive_order = int(sp.floor(-generator_weight))
                if max_nonpositive_order >= 0:
                    for order in range(max_nonpositive_order + 1):
                        operator = generator if order == 0 else d(generator, order)
                        operator_weight = _get_conformal_weight(operator)
                        if operator_weight is not None and operator_weight <= 0:
                            compensation += operator_weight

        return sp.simplify(compensation)

    def _factor_tuples(
        self,
        remaining_weight: sp.Expr,
        start_index: int,
        atomic_operators: list[tuple[Operator, sp.Expr, bool]],
        max_occurence: int | None = None,
        nonpositive_bosonic_counts: tuple[int, ...] = (),
    ) -> list[tuple[Operator, ...]]:
        atomic_ops = tuple(atomic_operators)

        @lru_cache(maxsize=None)
        def helper(
            remaining: sp.Expr,
            start: int,
            counts: tuple[int, ...],
        ) -> tuple[tuple[Operator, ...], ...]:
            tuples: list[tuple[Operator, ...]] = [tuple()] if remaining == 0 else []

            for index in range(start, len(atomic_ops)):
                operator, operator_weight, fermionic = atomic_ops[index]
                if operator_weight <= 0 and not fermionic:
                    if max_occurence is None:
                        continue
                    if counts[index] >= max_occurence:
                        continue

                next_remaining = sp.simplify(remaining - operator_weight)
                next_start = index + 1 if fermionic else index
                next_counts = counts
                if operator_weight <= 0 and not fermionic:
                    mutable_counts = list(counts)
                    mutable_counts[index] += 1
                    next_counts = tuple(mutable_counts)

                min_future_compensation = self._min_suffix_compensation(
                    next_start,
                    atomic_operators,
                    max_occurence=max_occurence,
                    nonpositive_bosonic_counts=next_counts,
                )
                if next_remaining < min_future_compensation:
                    continue

                for suffix in helper(next_remaining, next_start, next_counts):
                    tuples.append((operator,) + suffix)

            return tuple(tuples)

        return list(helper(remaining_weight, start_index, nonpositive_bosonic_counts))

    def _min_suffix_compensation(
        self,
        start_index: int,
        atomic_operators: list[tuple[Operator, sp.Expr, bool]],
        max_occurence: int | None = None,
        nonpositive_bosonic_counts: tuple[int, ...] = (),
    ) -> sp.Expr:
        compensation = sp.Integer(0)

        for index in range(start_index, len(atomic_operators)):
            _, operator_weight, fermionic = atomic_operators[index]
            if operator_weight >= 0:
                continue

            if fermionic:
                compensation += operator_weight
                continue

            if max_occurence is None:
                continue

            used = (
                nonpositive_bosonic_counts[index] if nonpositive_bosonic_counts else 0
            )
            remaining = max(0, max_occurence - used)
            compensation += remaining * operator_weight

        return sp.simplify(compensation)


class DescendantSpace:
    """Generate descendant subspaces from source operators at fixed weight."""

    def __init__(self, basis_builder: LocalOperatorBasis):
        self.basis_builder = basis_builder

    def generate(self, source: Any, target_weight: Any) -> list[Any]:
        """Generate canonical descendant candidates from one source."""
        normalized_target = _normalize_weight(target_weight)
        canonical_source = self.basis_builder.canonicalize(source)
        source_weight = _get_conformal_weight(canonical_source)

        if source_weight is None:
            raise ValueError(
                "Source operator must have a well-defined conformal weight"
            )
        if normalized_target < source_weight:
            return []

        queue = deque([canonical_source])
        seen = {canonical_source}
        generated: set[Any] = set()

        while queue:
            current = queue.popleft()
            current_weight = _get_conformal_weight(current)
            if current_weight is None or current_weight > normalized_target:
                continue

            if _weights_equal(current_weight, normalized_target):
                generated.add(current)
                continue

            derivative_descendant = self.basis_builder.canonicalize(d(current))
            self._enqueue_if_relevant(
                derivative_descendant,
                normalized_target,
                seen,
                queue,
            )

            for generator in self.basis_builder.generators:
                product_descendant = self.basis_builder.canonicalize(
                    normal_product(generator, current)
                )
                self._enqueue_if_relevant(
                    product_descendant,
                    normalized_target,
                    seen,
                    queue,
                )

        return sorted(generated, key=sp.srepr)

    def basis(self, source: Any, target_weight: Any) -> list[Any]:
        """Return a linearly independent spanning set for one source."""
        return self._independent_span(
            self.generate(source, target_weight), target_weight
        )

    def span(self, sources: Iterable[Any], target_weight: Any) -> list[Any]:
        """Return a linearly independent spanning set from multiple sources."""
        generated = []
        for source in sources:
            generated.extend(self.generate(source, target_weight))
        return self._independent_span(generated, target_weight)

    def _enqueue_if_relevant(
        self,
        expr: Any,
        target_weight: sp.Expr,
        seen: set[Any],
        queue: deque[Any],
    ) -> None:
        if expr == Zero:
            return

        expr_weight = _get_conformal_weight(expr)
        if expr_weight is None or expr_weight > target_weight:
            return
        if expr not in seen:
            seen.add(expr)
            queue.append(expr)

    def _independent_span(
        self, expressions: Iterable[Any], target_weight: Any
    ) -> list[Any]:
        normalized_target = _normalize_weight(target_weight)
        independent = []
        columns = []

        for expr in sorted(set(expressions), key=sp.srepr):
            vector = self.basis_builder.coordinates(expr, weight=normalized_target)
            if not columns:
                if vector != sp.zeros(*vector.shape):
                    independent.append(expr)
                    columns.append(vector)
                continue

            current = sp.Matrix.hstack(*columns)
            candidate = sp.Matrix.hstack(*columns, vector)
            if candidate.rank() > current.rank():
                independent.append(expr)
                columns.append(vector)

        return independent


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
            list(ansatz) if ansatz is not None else self.basis_builder.basis(weight)
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


class C2Space:
    """Construct the fixed-weight C2 subspace."""

    def __init__(self, basis_builder: LocalOperatorBasis):
        self.basis_builder = basis_builder

    def generators(self, weight: Any) -> list[Any]:
        """Generate canonical C2 candidates of the form NO(d(a), phi)."""
        normalized_weight = _normalize_weight(weight)
        generated = set()

        for generator in self.basis_builder.generators:
            generator_weight = _get_conformal_weight(generator)
            if generator_weight is None:
                continue

            phi_weight = sp.simplify(normalized_weight - generator_weight - 1)
            if phi_weight < 0:
                continue

            for phi in self.basis_builder.basis(phi_weight):
                candidate = self.basis_builder.canonicalize(
                    normal_product(d(generator), phi)
                )
                if candidate == Zero:
                    continue

                candidate_weight = _get_conformal_weight(candidate)
                if candidate_weight is None or not _weights_equal(
                    candidate_weight, normalized_weight
                ):
                    continue

                generated.add(candidate)

        return sorted(generated, key=sp.srepr)

    def basis(self, weight: Any) -> list[Any]:
        """Return a linearly independent basis for the fixed-weight C2 subspace."""
        normalized_weight = _normalize_weight(weight)
        independent = []
        columns = []

        for expr in self.generators(normalized_weight):
            vector = self.basis_builder.coordinates(expr, weight=normalized_weight)
            if not columns:
                if vector != sp.zeros(*vector.shape):
                    independent.append(expr)
                    columns.append(vector)
                continue

            current = sp.Matrix.hstack(*columns)
            candidate = sp.Matrix.hstack(*columns, vector)
            if candidate.rank() > current.rank():
                independent.append(expr)
                columns.append(vector)

        return independent

    def contains(self, expr: Any, weight: Any = None) -> bool:
        """Test whether an expression lies in the fixed-weight C2 subspace."""
        canonical = self.basis_builder.canonicalize(expr)
        if canonical == Zero:
            return True

        expr_weight = _get_conformal_weight(canonical)
        if expr_weight is None and weight is None:
            raise ValueError("Could not infer conformal weight; please provide weight")

        target_weight = expr_weight if weight is None else _normalize_weight(weight)
        if expr_weight is not None and not _weights_equal(expr_weight, target_weight):
            raise ValueError("Expression weight does not match requested C2 weight")

        basis_vectors = []
        for generator in self.basis(target_weight):
            basis_vectors.append(
                self.basis_builder.coordinates(generator, weight=target_weight)
            )

        expr_vector = self.basis_builder.coordinates(canonical, weight=target_weight)
        if not basis_vectors:
            return expr_vector == sp.zeros(*expr_vector.shape)

        matrix = sp.Matrix.hstack(*basis_vectors)
        augmented = sp.Matrix.hstack(matrix, expr_vector)
        return augmented.rank() == matrix.rank()


class C2NullSearcher:
    """Search for descendant operators equivalent to a target modulo C2."""

    def __init__(
        self,
        basis_builder: LocalOperatorBasis,
        descendants: Optional[DescendantSpace] = None,
        c2_space: Optional[C2Space] = None,
        stress_tensor: Any = None,
    ):
        self.basis_builder = basis_builder
        self.descendants = descendants or DescendantSpace(basis_builder)
        self.c2_space = c2_space or C2Space(basis_builder)
        self.stress_tensor = stress_tensor

    def search_from_sources(
        self,
        target_weight: Any,
        sources: Iterable[Any],
        target_expr: Any,
    ) -> Optional[dict[str, Any]]:
        """Solve D x - C y = t at fixed weight."""
        normalized_weight = _normalize_weight(target_weight)
        target_canonical = self.basis_builder.canonicalize(target_expr)
        target_vector = self.basis_builder.coordinates(
            target_canonical, normalized_weight
        )

        descendant_basis = self.descendants.span(sources, normalized_weight)
        c2_basis = self.c2_space.basis(normalized_weight)

        descendant_matrix = self._coordinate_matrix(descendant_basis, normalized_weight)
        c2_matrix = self._coordinate_matrix(c2_basis, normalized_weight)

        blocks = []
        if descendant_matrix.cols:
            blocks.append(descendant_matrix)
        if c2_matrix.cols:
            blocks.append(-c2_matrix)

        if not blocks:
            return None

        system_matrix = sp.Matrix.hstack(*blocks)
        if (
            system_matrix.rank()
            != sp.Matrix.hstack(system_matrix, target_vector).rank()
        ):
            return None

        solution_vector = self._solve_one_solution(system_matrix, target_vector)
        descendant_coeffs = [
            solution_vector[i, 0] for i in range(descendant_matrix.cols)
        ]
        c2_coeffs = [
            solution_vector[i, 0]
            for i in range(descendant_matrix.cols, solution_vector.rows)
        ]

        null_operator = self._linear_combination(descendant_basis, descendant_coeffs)
        c2_remainder = self._linear_combination(c2_basis, c2_coeffs)

        return {
            "weight": normalized_weight,
            "target": target_canonical,
            "target_vector": target_vector,
            "descendant_basis": descendant_basis,
            "descendant_matrix": descendant_matrix,
            "descendant_coeffs": descendant_coeffs,
            "c2_basis": c2_basis,
            "c2_matrix": c2_matrix,
            "c2_coeffs": c2_coeffs,
            "solution_vector": solution_vector,
            "null_operator": self.basis_builder.canonicalize(null_operator),
            "c2_remainder": self.basis_builder.canonicalize(c2_remainder),
        }

    def search_stress_tensor_nilpotency(
        self, n: int, sources: Iterable[Any], stress_tensor: Any = None
    ) -> Optional[dict[str, Any]]:
        """Search for a descendant congruent to T^(n) modulo C2."""
        tensor = stress_tensor if stress_tensor is not None else self.stress_tensor
        if tensor is None:
            raise ValueError("stress_tensor must be provided")

        target_expr = normal_product(*([tensor] * n))
        target_weight = _get_conformal_weight(target_expr)
        if target_weight is None:
            raise ValueError("Could not infer target weight for stress tensor product")

        result = self.search_from_sources(target_weight, sources, target_expr)
        if result is not None:
            result["n"] = n
        return result

    def _coordinate_matrix(self, expressions: Iterable[Any], weight: Any) -> sp.Matrix:
        columns = [
            self.basis_builder.coordinates(expr, weight=weight) for expr in expressions
        ]
        if not columns:
            basis_dim = len(self.basis_builder.basis(weight))
            return sp.zeros(basis_dim, 0)
        return sp.Matrix.hstack(*columns)

    def _linear_combination(self, basis: Iterable[Any], coeffs: Iterable[Any]) -> Any:
        terms = []
        for coeff, expr in zip(coeffs, basis):
            coeff = sp.sympify(coeff)
            if coeff != 0:
                terms.append(coeff * expr)
        if not terms:
            return Zero
        return sp.Add(*terms)

    def _solve_one_solution(self, matrix: sp.Matrix, rhs: sp.Matrix) -> sp.Matrix:
        solution = matrix.gauss_jordan_solve(rhs)
        if isinstance(solution, tuple):
            vector = solution[0]
            parameters = solution[1]
            if parameters.shape[0] > 0:
                substitutions = {param: 0 for param in parameters}
                vector = vector.xreplace(substitutions)
            return vector
        return solution
