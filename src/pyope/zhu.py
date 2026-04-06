"""Research-layer APIs for Zhu's associative algebra $A(V)=V/O(V)$."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any

import sympy as sp

from .api import NO, bracket
from .constants import One, Zero
from .local_operator import collect_operator_terms
from .operator_spaces import (
    LocalOperatorBasis,
    SparseLinearContext,
    _c2_weight_step,
    _get_conformal_weight,
    _normalize_weight,
)
from .operators import d


def _normalized_homogeneous_weight(expr: Any) -> sp.Expr:
    weight = _get_conformal_weight(expr)
    if weight is None:
        raise ValueError("Could not infer conformal weight for Zhu product")
    return _normalize_weight(weight)


def _homogeneous_components(
    canonicalizer: LocalOperatorBasis, expr: Any
) -> dict[sp.Expr, Any]:
    canonical = canonicalizer.canonicalize(expr)
    if canonical == Zero:
        return {}

    grouped: dict[sp.Expr, list[Any]] = {}
    for operator, coeff in collect_operator_terms(canonical).items():
        monomial = One if operator == 1 else operator
        weight = _get_conformal_weight(monomial)
        if weight is None:
            raise ValueError("Could not infer conformal weight for Zhu decomposition")
        grouped.setdefault(_normalize_weight(weight), []).append(
            sp.sympify(coeff) * monomial
        )

    return {
        weight: canonicalizer.canonicalize(sp.Add(*terms))
        for weight, terms in grouped.items()
    }


def _homogeneous_component(
    canonicalizer: LocalOperatorBasis, expr: Any, target_weight: Any
) -> Any:
    components = _homogeneous_components(canonicalizer, expr)
    return canonicalizer.canonicalize(
        components.get(_normalize_weight(target_weight), Zero)
    )


def zhu_star_product(left: Any, right: Any) -> Any:
    """Return Zhu's star product $a*b$ for homogeneous ``left``.

    Formula:
        a * b = sum_{i >= 0} binom(wt(a), i) a_{i-1} b

    Using pyope conventions, ``a_{-1} b = NO(a, b)`` and for ``i >= 1`` we use
    the singular OPE coefficient ``{ab}_{i}`` accessed via ``bracket(a, b, i)``.
    """

    weight = _normalized_homogeneous_weight(left)
    result = NO(left, right)
    ope = None
    from .api import OPE

    ope = OPE(left, right)
    for pole in range(1, ope.max_pole + 1):
        coeff = sp.binomial(weight, pole)
        if coeff == 0:
            continue
        term = bracket(left, right, pole)
        if term == Zero:
            continue
        result = result + coeff * term
    return result


def zhu_circle_product(left: Any, right: Any) -> Any:
    """Return Zhu's circle product generating $O(V)$.

    Formula:
        a circle b = sum_{i >= 0} binom(wt(a), i) a_{i-2} b

    We realize the ``i=0`` term as ``a_{-2} b = NO(d(a), b)``.
    For ``i >= 1``, ``a_{i-2} b`` is the ``(i-1)``-st OPE pole, with
    ``a_{-1} b = NO(a, b)`` when ``i=1``.
    """

    weight = _normalized_homogeneous_weight(left)
    result = NO(d(left), right)
    ope = None
    from .api import OPE

    ope = OPE(left, right)
    for index in range(1, ope.max_pole + 2):
        coeff = sp.binomial(weight, index)
        if coeff == 0:
            continue
        term = bracket(left, right, index - 1)
        if term == Zero:
            continue
        result = result + coeff * term
    return result


@dataclass
class ZhuReductionWitness:
    expr: Any
    remainder: Any
    generators: list[Any]
    coefficients: list[sp.Expr]
    ov_part: Any
    target_weight: sp.Expr | None


class AbstractZhuReducer(ABC):
    """Common interface for quotient reduction modulo Zhu's $O(V)$."""

    @abstractmethod
    def quotient_normal_form(self, expr: Any) -> Any:
        raise NotImplementedError

    def is_zero_mod_ov(self, expr: Any) -> bool:
        return bool(self.quotient_normal_form(expr) == Zero)

    @abstractmethod
    def solve_zhu_witness(self, expr: Any) -> ZhuReductionWitness:
        raise NotImplementedError


class GenericZhuReducer(AbstractZhuReducer):
    """Sparse reducer for Zhu's algebra quotient modulo $O(V)$."""

    def __init__(
        self,
        canonicalizer: LocalOperatorBasis,
        linear_context: SparseLinearContext | None = None,
        generator_provider: Any = None,
    ):
        self.canonicalizer = canonicalizer
        self.linear_context = linear_context or SparseLinearContext(canonicalizer)
        self.generator_provider = generator_provider

    def candidate_generators_for_expr(self, expr: Any) -> list[Any]:
        canonical_expr = self.canonicalizer.canonicalize(expr)
        if canonical_expr == Zero:
            return []

        target_weight = _get_conformal_weight(canonical_expr)
        if target_weight is None:
            raise ValueError("Could not infer conformal weight for Zhu reduction")
        target_weight = _normalize_weight(target_weight)

        return self.candidate_generators_for_weight(target_weight)

    def candidate_generators_for_weight(self, target_weight: Any) -> list[Any]:
        target_weight = _normalize_weight(target_weight)

        candidates: set[Any] = set()
        step = _c2_weight_step(self.canonicalizer.generators)

        if self.generator_provider is not None:
            left_basis = self.generator_provider
        else:
            left_basis = []
            h_a = step
            while sp.simplify(h_a) <= target_weight:
                left_basis.extend(self.canonicalizer.list(h_a))
                h_a = sp.simplify(h_a + step)

        for left in left_basis:
            left_weight = _get_conformal_weight(left)
            if left_weight is None:
                continue
            right_weight = sp.simplify(target_weight - left_weight)
            if right_weight < 0:
                continue
            for right in self.canonicalizer.list(right_weight):
                candidate = _homogeneous_component(
                    self.canonicalizer,
                    zhu_circle_product(left, right),
                    target_weight,
                )
                if candidate == Zero:
                    continue
                candidates.add(candidate)

        return self.linear_context.independent_subset(sorted(candidates, key=sp.srepr))

    def quotient_normal_form(self, expr: Any) -> Any:
        return self.solve_zhu_witness(expr).remainder

    def solve_zhu_witness(self, expr: Any) -> ZhuReductionWitness:
        canonical_expr = self.canonicalizer.canonicalize(expr)
        if canonical_expr == Zero:
            return ZhuReductionWitness(
                expr=expr,
                remainder=Zero,
                generators=[],
                coefficients=[],
                ov_part=Zero,
                target_weight=None,
            )

        components = _homogeneous_components(self.canonicalizer, canonical_expr)
        if len(components) == 1:
            target_weight = next(iter(components))
            return self._solve_homogeneous_zhu_witness(
                components[target_weight], target_weight
            )

        witnesses = [
            self._solve_homogeneous_zhu_witness(component, weight)
            for weight, component in sorted(
                components.items(), key=lambda item: sp.srepr(item[0])
            )
        ]

        return ZhuReductionWitness(
            expr=expr,
            remainder=self.canonicalizer.canonicalize(
                sp.Add(*[witness.remainder for witness in witnesses])
            ),
            generators=[
                generator for witness in witnesses for generator in witness.generators
            ],
            coefficients=[
                coeff for witness in witnesses for coeff in witness.coefficients
            ],
            ov_part=self.canonicalizer.canonicalize(
                sp.Add(*[witness.ov_part for witness in witnesses])
            ),
            target_weight=None,
        )

    def _solve_homogeneous_zhu_witness(
        self, canonical_expr: Any, target_weight: Any
    ) -> ZhuReductionWitness:
        target_weight = _normalize_weight(target_weight)
        candidates = self.candidate_generators_for_weight(target_weight)
        expr_terms = self.linear_context.sparse_terms(canonical_expr)

        eliminator = SparseLinearContext(self.canonicalizer)
        independent_candidates: list[Any] = []
        candidate_vectors: list[dict[Any, sp.Expr]] = []
        for candidate in candidates:
            vector = eliminator.sparse_terms(candidate)
            reduced, _ = eliminator.reduce_vector(vector)
            if not reduced:
                continue
            eliminator._eliminator.insert_reduced(reduced)
            independent_candidates.append(candidate)
            candidate_vectors.append(vector)

        reduced_terms, coeff_map = eliminator.reduce_vector(expr_terms)
        coeffs = []
        for pivot, candidate_vec in zip(
            eliminator._eliminator.pivot_order, candidate_vectors
        ):
            raw_coeff = sp.sympify(coeff_map.get(pivot, 0))
            pivot_val_in_candidate = sp.sympify(candidate_vec.get(pivot, 0))
            if pivot_val_in_candidate != 0:
                coeffs.append(raw_coeff / pivot_val_in_candidate)
            else:
                coeffs.append(raw_coeff)

        ov_part_terms = []
        for coeff, candidate in zip(coeffs, independent_candidates):
            if coeff != 0:
                ov_part_terms.append(coeff * candidate)
        ov_part = self.canonicalizer.canonicalize(
            sp.Add(*ov_part_terms) if ov_part_terms else Zero
        )

        remainder_terms = [
            coeff * monomial for monomial, coeff in reduced_terms.items()
        ]
        remainder = self.canonicalizer.canonicalize(
            sp.Add(*remainder_terms) if remainder_terms else Zero
        )

        return ZhuReductionWitness(
            expr=canonical_expr,
            remainder=remainder,
            generators=independent_candidates,
            coefficients=coeffs,
            ov_part=ov_part,
            target_weight=target_weight,
        )


class ZhuSpace:
    """Convenience facade for fixed-weight computations modulo Zhu's $O(V)$."""

    def __init__(self, basis_builder: LocalOperatorBasis):
        self.basis_builder = basis_builder
        self._compat_reducer: GenericZhuReducer | None = None

    def reducer(self) -> GenericZhuReducer:
        if self._compat_reducer is None:
            self._compat_reducer = GenericZhuReducer(self.basis_builder)
        return self._compat_reducer

    def generators(self, weight: Any) -> list[Any]:
        normalized_weight = _normalize_weight(weight)
        return self.reducer().candidate_generators_for_weight(normalized_weight)

    def basis(self, weight: Any) -> list[Any]:
        normalized_weight = _normalize_weight(weight)
        return self.basis_builder.list_independent_ops(
            self.generators(normalized_weight), weight=normalized_weight
        )

    def contains(self, expr: Any, weight: Any = None) -> bool:
        canonical = self.basis_builder.canonicalize(expr)
        if canonical == Zero:
            return True

        if weight is None:
            return bool(self.reducer().is_zero_mod_ov(canonical))

        target_weight = _normalize_weight(weight)
        component = _homogeneous_component(self.basis_builder, canonical, target_weight)
        return self.is_zero_mod_ov(component, weight=target_weight)

    def quotient_normal_form(self, expr: Any, weight: Any = None) -> Any:
        canonical = self.basis_builder.canonicalize(expr)
        if canonical == Zero:
            return Zero
        if weight is not None:
            canonical = _homogeneous_component(
                self.basis_builder, canonical, _normalize_weight(weight)
            )
        return self.reducer().quotient_normal_form(canonical)

    def is_zero_mod_ov(self, expr: Any, weight: Any = None) -> bool:
        return bool(self.quotient_normal_form(expr, weight=weight) == Zero)

    def solve_zhu_witness(self, expr: Any, weight: Any = None) -> ZhuReductionWitness:
        canonical = self.basis_builder.canonicalize(expr)
        if weight is not None:
            canonical = _homogeneous_component(
                self.basis_builder, canonical, _normalize_weight(weight)
            )
        return self.reducer().solve_zhu_witness(canonical)

    def multiply(self, left: Any, right: Any, weight: Any = None) -> Any:
        return self.quotient_normal_form(zhu_star_product(left, right), weight=weight)
