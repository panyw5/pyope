"""Operator-space utilities for fixed-weight and sparse local-operator computations."""

from __future__ import annotations

import inspect
import shutil
from functools import cmp_to_key, lru_cache
from typing import Any, Iterable, MutableMapping, Optional, cast

import sympy as sp

from .backend import get_compute_backend
from .api import normal_product
from .constants import One, Zero
from .local_operator import (
    collect_operator_terms,
    extract_scalar_operator,
    get_operator_parity,
)
from .operators import BasisOperator, NormalOrderedOperator, Operator, d
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
        object.__setattr__(obj, "_realization", realization)
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
        generator = cast(
            RealizedGenerator,
            RealizedGenerator(name, realization=realization, **assumptions),
        )
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


def _ordered_unique_expressions(expressions: Iterable[Any]) -> list[Any]:
    return sorted(set(expressions), key=sp.srepr)


def _should_use_wolfram_precanonicalization() -> bool:
    return (
        get_compute_backend() == "wolfram" and shutil.which("wolframscript") is not None
    )


def _batch_wolfram_precanonicalize(expressions: list[Any]) -> list[Any]:
    if not expressions:
        return []

    from .wolfram_backend import chunk_exprs_for_wolfram, simplify_expr

    canonicalized: list[Any] = []
    for chunk in chunk_exprs_for_wolfram(expressions):
        canonicalized.extend(simplify_expr(chunk))
    return canonicalized


def _precanonicalize_expressions(expressions: Iterable[Any]) -> list[Any]:
    expr_list = list(expressions)
    if not expr_list or not _should_use_wolfram_precanonicalization():
        return expr_list

    return _batch_wolfram_precanonicalize(expr_list)


def _precanonicalize_expression(expr: Any) -> Any:
    return _precanonicalize_expressions([expr])[0]


def _canonicalization_mode() -> str:
    if _should_use_wolfram_precanonicalization():
        return "wolfram-pre"
    return get_compute_backend()


def _precanonicalized_expression_lookup(expressions: Iterable[Any]) -> dict[Any, Any]:
    expr_list = list(expressions)
    canonicalized = _precanonicalize_expressions(expr_list)
    return dict(zip(expr_list, canonicalized, strict=True))


def _terms_from_precanonical_expression(expr: Any) -> dict[Any, sp.Expr]:
    if expr == 0 or expr == Zero:
        return {}

    terms = {}
    for operator, coeff in collect_operator_terms(expr).items():
        coeff = sp.sympify(coeff)
        if coeff == 0:
            continue

        monomial = One if operator == 1 else operator
        terms[monomial] = sp.sympify(terms.get(monomial, 0) + coeff)

    return {monomial: coeff for monomial, coeff in terms.items() if coeff != 0}


def _coordinates_from_precanonical_expression(
    expr: Any,
    basis_builder: "LocalOperatorBasis",
    weight: Any,
    max_occurence: Any = None,
) -> sp.Matrix:
    normalized_max_occurence = basis_builder._resolve_max_occurence(max_occurence)
    target_weight = _normalize_weight(weight)

    if expr == 0 or expr == Zero:
        basis = basis_builder.list(
            target_weight, max_occurence=normalized_max_occurence
        )
        return sp.zeros(len(basis), 1)

    expr_weight = _get_conformal_weight(expr)
    if expr_weight is not None and not _weights_equal(expr_weight, target_weight):
        raise ValueError("Expression weight does not match requested basis weight")

    basis = basis_builder.list(target_weight, max_occurence=normalized_max_occurence)
    index = {op: i for i, op in enumerate(basis)}
    vector = sp.zeros(len(basis), 1)

    for operator, coeff in collect_operator_terms(expr).items():
        if operator == 1:
            operator = One
        if operator not in index:
            raise ValueError(
                f"Operator {operator} is not in the canonical basis at weight {target_weight}"
            )
        vector[index[operator], 0] += sp.sympify(coeff)

    return vector


def _independent_from_vectors(
    expressions: Iterable[Any], vector_getter: Any
) -> list[Any]:
    ordered = _ordered_unique_expressions(expressions)
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


class SparseLinearContext:
    """On-demand sparse linear algebra over canonical monomial supports."""

    def __init__(self, canonicalizer: Any):
        self.canonicalizer = canonicalizer
        self._eliminator = _SparseIndependentEliminator()

    def sparse_terms(self, expr: Any) -> dict[Any, sp.Expr]:
        sparse_terms = getattr(self.canonicalizer, "sparse_terms", None)
        if callable(sparse_terms):
            return cast(dict[Any, sp.Expr], sparse_terms(expr))
        return _canonical_terms(expr, self.canonicalizer.canonicalize)

    def reduce_vector(
        self, sparse_terms: dict[Any, sp.Expr]
    ) -> tuple[dict[Any, sp.Expr], dict[Any, sp.Expr]]:
        return self._eliminator.reduce_with_coefficients(sparse_terms)

    def insert_expr(self, expr: Any) -> bool:
        return self._eliminator.insert(self.sparse_terms(expr))

    def independent_subset(self, expressions: Iterable[Any]) -> list[Any]:
        ordered = _ordered_unique_expressions(expressions)
        eliminator = _SparseIndependentEliminator()
        independent = []

        for expr in ordered:
            if eliminator.insert(self.sparse_terms(expr)):
                independent.append(expr)

        return independent

    def zero_relations(self, expressions: Iterable[Any]) -> list[dict[str, Any]]:
        return _zero_relations_from_vectors(expressions, self.sparse_terms)


def _canonical_terms(expr: Any, canonicalizer: Any) -> dict[Any, sp.Expr]:
    """Return a sparse canonical expansion of an operator expression."""
    canonical = canonicalizer(expr)  # OPE computation entry
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


def _compare_operators_for_sort(left: Operator, right: Operator) -> int:
    """Adapter from compare_operators() to Python's cmp-style sort."""
    comparison = ope_registry.compare_operators(left, right)
    if comparison > 0:
        return -1
    if comparison < 0:
        return 1
    return 0


def _build_canonical_basis_monomial(factors: Iterable[Operator]) -> Any:
    """Build a canonical right-nested monomial without simplification."""
    factor_tuple = tuple(factors)
    if not factor_tuple:
        return One
    if len(factor_tuple) == 1:
        return factor_tuple[0]

    result: Operator = factor_tuple[-1]
    for factor in reversed(factor_tuple[:-1]):
        result = NormalOrderedOperator(factor, result)
    return result


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


_realize_cache: dict[Any, Any] = {}


def clear_realize_cache() -> None:
    """Clear the realization memoization cache."""
    _realize_cache.clear()


def _realize_expr(expr: Any) -> Any:
    """Recursively expand realized generators into their underlying expressions.

    Optimizations over naive recursive expansion:
    - **P0 – incremental simplification**: after each ``normal_product`` call the
      result is immediately combined (like terms merged) so that intermediate
      expression sizes stay small instead of snowballing.
    - **P1 – RealizedGenerator caching**: each ``RealizedGenerator`` is expanded
      at most once; subsequent encounters reuse the cached result.
    """
    # P1: memoize RealizedGenerator expansions
    if isinstance(expr, RealizedGenerator):
        cached = _realize_cache.get(expr)
        if cached is not None:
            return cached
        result = _realize_expr(expr.realization)
        _realize_cache[expr] = result
        return result

    if expr is Zero or expr is One:
        return expr

    if isinstance(expr, sp.Add):
        # P0: combine like terms after expanding each summand
        expanded = sp.Add(*[_realize_expr(arg) for arg in expr.args])
        return _combine_like_terms_preserving_metadata(expanded)

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
            # P0: combine like terms immediately after normal ordering
            result = normal_product(_realize_expr(left), _realize_expr(right))
            return _combine_like_terms_preserving_metadata(result)
        return expr

    return expr


def realize(expr: Any) -> Any:
    """Expand realized generators in an operator expression."""
    return expr.realize()


def realize_and_simplify(expr: Any) -> Any:
    """Expand realized generators and canonicalize the resulting expression."""
    realized = realize(expr)
    precanonicalized = _precanonicalize_expression(realized)
    if _should_use_wolfram_precanonicalization():
        return _combine_like_terms_preserving_metadata(precanonicalized)
    return _combine_like_terms_preserving_metadata(simplify(precanonicalized))


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

    context = SparseLinearContext(basis_builder)
    use_wolfram_precanonicalization = _should_use_wolfram_precanonicalization()
    processed_lookup: dict[Any, Any] = {}
    if use_wolfram_precanonicalization:
        processed_lookup = _precanonicalized_expression_lookup(
            _ordered_unique_expressions(expressions)
        )

    def vector_getter(expr: Any) -> dict[Any, sp.Expr]:
        _validate_expression_weight(expr, weight)
        if not use_wolfram_precanonicalization:
            return context.sparse_terms(expr)
        return _terms_from_precanonical_expression(processed_lookup[expr])

    return _independent_from_vectors(expressions, vector_getter)


def list_zero_relations(
    expressions: Iterable[Any],
    basis_builder: "LocalOperatorBasis",
    weight: Any = None,
    max_occurence: Any = None,
) -> list[dict[str, Any]]:
    """Return a basis of linear relations via local canonical expansions."""

    context = SparseLinearContext(basis_builder)
    use_wolfram_precanonicalization = _should_use_wolfram_precanonicalization()
    processed_lookup: dict[Any, Any] = {}
    if use_wolfram_precanonicalization:
        processed_lookup = _precanonicalized_expression_lookup(expressions)

    def vector_getter(expr: Any) -> dict[Any, sp.Expr]:
        _validate_expression_weight(expr, weight)
        if not use_wolfram_precanonicalization:
            return context.sparse_terms(expr)
        return _terms_from_precanonical_expression(processed_lookup[expr])

    return _zero_relations_from_vectors(expressions, vector_getter)


def independent_under_realization(
    expressions: Iterable[Any],
    free_field_basis: "LocalOperatorBasis",
    weight: Any = None,
    max_occurence: Any = None,
) -> list[Any]:
    """Return expressions that stay linearly independent after realization."""

    ordered = _ordered_unique_expressions(expressions)
    use_wolfram_precanonicalization = _should_use_wolfram_precanonicalization()
    processed_lookup: dict[Any, Any] = {}

    if use_wolfram_precanonicalization and ordered:
        realized_exprs = [realize(expr) for expr in ordered]
        canonicalized = _batch_wolfram_precanonicalize(realized_exprs)
        processed_lookup = dict(zip(ordered, canonicalized, strict=True))

    def vector_getter(expr: Any) -> sp.Matrix:
        target_weight = weight
        if target_weight is None:
            target_weight = _get_conformal_weight(expr)
            if target_weight is None:
                raise ValueError(
                    "Could not infer conformal weight under realization; provide weight"
                )

        if not use_wolfram_precanonicalization:
            return realized_coordinates(
                expr,
                free_field_basis,
                weight=target_weight,
                max_occurence=max_occurence,
            )

        return _coordinates_from_precanonical_expression(
            processed_lookup[expr],
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
        stress_tensor: Any = None,
        gradings: Any = None,
        max_occurence: Any = None,
    ):
        raw_generators = tuple(generators)
        self.stress_tensor = stress_tensor
        self.gradings = gradings
        self._basis = self
        self.max_occurence = _normalize_max_occurence(max_occurence)

        if not raw_generators:
            raise ValueError("LocalOperatorBasis requires at least one generator")

        has_nonpositive_generator = False
        for generator in raw_generators:
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

        self.generators = tuple(
            sorted(raw_generators, key=cmp_to_key(_compare_operators_for_sort))
        )

    def _resolve_max_occurence(self, value: Any = None) -> int | None:
        normalized = _normalize_max_occurence(value)
        if normalized is None:
            return self.max_occurence
        return normalized

    def canonicalize(self, expr: Any) -> Any:
        """Canonicalize an operator expression for basis comparisons."""
        simplified = self._canonicalize_cached(
            expr,
            ope_registry.version,
            _canonicalization_mode(),
        )
        combined = _combine_like_terms_preserving_metadata(simplified)
        if combined == 0:
            return Zero
        return combined

    @staticmethod
    @lru_cache(maxsize=None)
    def _canonicalize_cached(
        expr: Any, registry_version: int, canonicalization_mode: str
    ) -> Any:
        del registry_version
        processed = expr
        if canonicalization_mode == "wolfram-pre":
            processed = _precanonicalize_expression(expr)
            return processed
        return simplify(processed)

    def sparse_terms(self, expr: Any) -> dict[Any, sp.Expr]:
        """Return sparse canonical coordinates without fixed-weight basis enumeration."""
        return _canonical_terms(expr, self.canonicalize)

    def basis(self, weight: Any, max_occurence: Any = None) -> list[Any]:
        """Alias for list(); retained for interface compatibility."""
        return self.list(weight, max_occurence=max_occurence)

    def sector_of(self, expr: Any) -> dict[str, Any]:
        """Return the conformal-weight/parity sector of an expression."""
        weight = _get_conformal_weight(expr)
        parity = get_operator_parity(expr)
        return {"weight": weight, "parity": parity}

    def nested_stress_tensor(self, n: int) -> Any:
        """Return the canonicalized n-fold normal product of the stress tensor."""
        if self.stress_tensor is None:
            raise ValueError("stress_tensor must be provided")
        if n <= 0:
            raise ValueError("n must be positive")
        return self.canonicalize(normal_product(*([self.stress_tensor] * n)))

    def enumerate_candidates(self, weight: Any) -> list[Any]:
        """Enumerate raw canonical candidate expressions of fixed weight."""
        normalized_weight = _normalize_weight(weight)
        basis_terms = self._basis_terms(
            normalized_weight, self._resolve_max_occurence()
        )
        return list(basis_terms)

    def list(self, weight: Any, max_occurence: Any = None) -> list[Any]:
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
            basis = self.list(weight, max_occurence=normalized_max_occurence)
            return sp.zeros(len(basis), 1)

        expr_weight = _get_conformal_weight(canonical)
        if expr_weight is None and weight is None:
            raise ValueError(
                "Could not infer conformal weight; please provide weight explicitly"
            )

        target_weight = expr_weight if weight is None else _normalize_weight(weight)
        if expr_weight is not None and not _weights_equal(expr_weight, target_weight):
            raise ValueError("Expression weight does not match requested basis weight")

        basis = self.list(target_weight, max_occurence=normalized_max_occurence)
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
        nonpositive_bosonic_counts = tuple(0 for _ in atomic_operators)
        basis_terms: list[Any] = []
        seen_terms: set[Any] = set()

        for factors in self._factor_tuples(
            weight,
            0,
            atomic_operators,
            max_occurence=max_occurence,
            nonpositive_bosonic_counts=nonpositive_bosonic_counts,
        ):
            monomial = _build_canonical_basis_monomial(factors)
            operator_weight = _get_conformal_weight(monomial)
            if operator_weight is None or not _weights_equal(operator_weight, weight):
                continue
            if weight == 0 and _is_negative_single_letter_fermion_sector(monomial):
                continue
            if monomial in seen_terms:
                continue
            seen_terms.add(monomial)
            basis_terms.append(monomial)

        return tuple(basis_terms)

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

        atoms.sort(key=lambda item: cmp_to_key(_compare_operators_for_sort)(item[0]))
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


class LocalOperatorCanonicalizer(LocalOperatorBasis):
    """Backward-compatible alias for LocalOperatorBasis.

    New code should use LocalOperatorBasis directly.
    """

    pass


def _c2_weight_step(generators: Any) -> sp.Expr:
    """Return the smallest weight unit implied by the given generators.

    For integer-weight generators this is 1; for half-integer generators it is
    1/2, etc.  The step controls how finely we sweep h_a when enumerating C2
    generators of the form NO(d(a), phi).
    """
    denoms: list[int] = []
    for g in generators:
        w = _get_conformal_weight(g)
        if w is None:
            continue
        r = sp.Rational(w)
        denoms.append(int(r.q))
    if not denoms:
        return sp.Integer(1)
    lcm_val = denoms[0]
    for denom in denoms[1:]:
        lcm_val = int(sp.lcm(lcm_val, denom))
    return sp.Rational(1, lcm_val)


class C2Space:
    """Construct the fixed-weight C2 subspace."""

    def __init__(self, basis_builder: LocalOperatorBasis):
        self.basis_builder = basis_builder
        self._compat_reducer = None

    def reducer(self):
        """Return a compatibility reducer backed by the new sparse C2 API."""
        if self._compat_reducer is None:
            from .c2 import GenericC2Reducer

            self._compat_reducer = GenericC2Reducer(cast(Any, self.basis_builder))
        return self._compat_reducer

    def generators(self, weight: Any) -> list[Any]:
        """Generate canonical C2 candidates of the form NO(d(a), phi).

        C2(V)_w = span{:(∂a)φ: | a ∈ V[h_a], φ ∈ V[h_φ], h_a + h_φ + 1 = w}.
        We iterate over ALL basis elements a (not just primary generators) so that
        derived elements like ∂^n T are correctly included.
        """
        normalized_weight = _normalize_weight(weight)
        generated = set()

        step = _c2_weight_step(self.basis_builder.generators)
        h_a = step
        while sp.simplify(normalized_weight - h_a - 1) >= 0:
            phi_weight = sp.simplify(normalized_weight - h_a - 1)
            a_list = self.basis_builder.list(h_a)
            phi_list = self.basis_builder.list(phi_weight)

            for a in a_list:
                for phi in phi_list:
                    candidate = self.basis_builder.canonicalize(
                        normal_product(d(a), phi)
                    )
                    if candidate == Zero:
                        continue

                    candidate_weight = _get_conformal_weight(candidate)
                    if candidate_weight is None or not _weights_equal(
                        candidate_weight, normalized_weight
                    ):
                        continue

                    generated.add(candidate)

            h_a = sp.simplify(h_a + step)

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

        return self.is_zero_mod_c2(canonical, weight=target_weight)

    def quotient_normal_form(self, expr: Any, weight: Any = None) -> Any:
        """Compatibility helper exposing reducer-style quotient normal forms."""
        canonical = self.basis_builder.canonicalize(expr)
        if canonical == Zero:
            return Zero
        if weight is not None:
            expr_weight = _get_conformal_weight(canonical)
            target_weight = _normalize_weight(weight)
            if expr_weight is not None and not _weights_equal(
                expr_weight, target_weight
            ):
                raise ValueError("Expression weight does not match requested C2 weight")
        return self.reducer().quotient_normal_form(canonical)

    def is_zero_mod_c2(self, expr: Any, weight: Any = None) -> bool:
        """Compatibility helper exposing reducer-style C2 triviality checks."""
        return self.quotient_normal_form(expr, weight=weight) == Zero

    def solve_c2_witness(self, expr: Any, weight: Any = None) -> Any:
        """Compatibility helper exposing reducer witness diagnostics."""
        canonical = self.basis_builder.canonicalize(expr)
        if weight is not None:
            expr_weight = _get_conformal_weight(canonical)
            target_weight = _normalize_weight(weight)
            if expr_weight is not None and not _weights_equal(
                expr_weight, target_weight
            ):
                raise ValueError("Expression weight does not match requested C2 weight")
        return self.reducer().solve_c2_witness(canonical)


class LegacyC2NullSearcher:
    """Search for descendant operators equivalent to a target modulo C2.

    Note: This is the legacy interface. For new code, prefer
    ``null_search.C2NullSearcher`` (exported as ``QuotientC2NullSearcher``).
    """

    def __init__(
        self,
        basis_builder: LocalOperatorBasis,
        descendants: Optional[Any] = None,
        c2_space: Optional[C2Space] = None,
        stress_tensor: Any = None,
        c2_reducer: Any = None,
    ):
        if descendants is None:
            from .descendants import DescendantSpace as ExternalDescendantSpace

            descendants = ExternalDescendantSpace(cast(Any, basis_builder))
        self.basis_builder = basis_builder
        self.descendants = descendants
        self.c2_space = c2_space or C2Space(basis_builder)
        self.stress_tensor = stress_tensor
        self.c2_reducer = c2_reducer or self.c2_space.reducer()

    def quotient_precheck(self, target_expr: Any) -> Any:
        """Bridge legacy searcher to the new reducer-based precheck API."""
        from .null_search import C2NullSearcher as QuotientC2NullSearcher

        canonicalizer = LocalOperatorBasis(
            self.basis_builder.generators,
            stress_tensor=self.stress_tensor,
            gradings=getattr(self.basis_builder, "gradings", None),
            max_occurence=self.basis_builder.max_occurence,
        )
        searcher = QuotientC2NullSearcher(
            canonicalizer=canonicalizer,
            linear_context=SparseLinearContext(canonicalizer),
            descendants=self.descendants,
            c2_reducer=self.c2_reducer,
        )
        return searcher.quotient_precheck(target_expr)

    def search_from_sources(
        self,
        target_weight: Any,
        sources: Iterable[Any],
        target_expr: Any,
    ) -> Optional[dict[str, Any]]:
        """Solve descendant lift via the new sparse quotient core."""
        from .null_search import C2NullSearcher as QuotientC2NullSearcher

        normalized_weight = _normalize_weight(target_weight)
        target_canonical = self.basis_builder.canonicalize(target_expr)
        target_vector = self.basis_builder.coordinates(
            target_canonical, normalized_weight
        )

        canonicalizer = LocalOperatorBasis(
            self.basis_builder.generators,
            stress_tensor=self.stress_tensor,
            gradings=getattr(self.basis_builder, "gradings", None),
            max_occurence=self.basis_builder.max_occurence,
        )
        searcher = QuotientC2NullSearcher(
            canonicalizer=canonicalizer,
            linear_context=SparseLinearContext(canonicalizer),
            descendants=self.descendants,
            c2_reducer=self.c2_reducer,
        )
        result = searcher.search_from_sources(
            normalized_weight,
            sources,
            target_canonical,
        )
        if result is None or result.status == "obstructed":
            return None
        return self._legacy_result_payload(result, override_weight=normalized_weight)

    def search_stress_tensor_nilpotency(
        self, n: int, sources: Iterable[Any], stress_tensor: Any = None
    ) -> Optional[dict[str, Any]]:
        """Search for a descendant congruent to T^(n) modulo C2."""
        tensor = stress_tensor if stress_tensor is not None else self.stress_tensor
        if tensor is None:
            raise ValueError("stress_tensor must be provided")

        from .null_search import C2NullSearcher as QuotientC2NullSearcher

        canonicalizer = LocalOperatorBasis(
            self.basis_builder.generators,
            stress_tensor=tensor,
            gradings=getattr(self.basis_builder, "gradings", None),
            max_occurence=self.basis_builder.max_occurence,
        )
        searcher = QuotientC2NullSearcher(
            canonicalizer=canonicalizer,
            linear_context=SparseLinearContext(canonicalizer),
            descendants=self.descendants,
            singular_constraints=self._make_singular_analyzer(
                tensor,
            ),
            c2_reducer=self.c2_reducer,
        )
        result = searcher.search_stress_tensor_nilpotency(n, sources)
        if result is None or result.status == "obstructed":
            return None
        return self._legacy_result_payload(result, override_n=n)

    def _make_singular_analyzer(self, tensor: Any) -> Any:
        from .singularity import (
            SingularVectorAnalyzer as ExternalSingularVectorAnalyzer,
        )

        return ExternalSingularVectorAnalyzer(
            cast(Any, self.basis_builder),
            generators=self.basis_builder.generators,
            stress_tensor=tensor,
        )

    def _legacy_result_payload(
        self,
        result: Any,
        override_weight: Any = None,
        override_n: int | None = None,
    ) -> dict[str, Any]:
        payload = result.legacy_payload(self.basis_builder)
        if override_weight is not None:
            payload["weight"] = override_weight
            payload["target_vector"] = self.basis_builder.coordinates(
                payload["target"],
                override_weight,
            )
            payload["descendant_matrix"] = self._coordinate_matrix(
                payload["descendant_basis"],
                override_weight,
            )
            payload["c2_matrix"] = self._coordinate_matrix(
                payload["c2_basis"],
                override_weight,
            )
        if override_n is not None:
            payload["n"] = override_n
        return payload

    def _coordinate_matrix(self, expressions: Iterable[Any], weight: Any) -> sp.Matrix:
        columns = [
            self.basis_builder.coordinates(expr, weight=weight) for expr in expressions
        ]
        if not columns:
            basis_dim = len(self.basis_builder.list(weight))
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


# Backward-compatible alias
C2NullSearcher = LegacyC2NullSearcher
