"""Sparse $C_2$ quotient reducers for local-operator expressions."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp

from .api import normal_product
from .constants import Zero
from .operator_spaces import (
    LocalOperatorBasis,
    SparseLinearContext,
    _c2_weight_step,
    _get_conformal_weight,
    _normalize_weight,
    _weights_equal,
)
from .operators import d


@dataclass
class C2ReductionWitness:
    expr: Any
    remainder: Any
    generators: list[Any]
    coefficients: list[sp.Expr]
    c2_part: Any
    target_weight: sp.Expr | None


class AbstractC2Reducer(ABC):
    """Common interface for quotient reduction modulo $C_2$."""

    @abstractmethod
    def quotient_normal_form(self, expr: Any) -> Any:
        raise NotImplementedError

    def is_zero_mod_c2(self, expr: Any) -> bool:
        return self.quotient_normal_form(expr) == Zero

    @abstractmethod
    def solve_c2_witness(self, expr: Any) -> C2ReductionWitness:
        raise NotImplementedError


class GenericC2Reducer(AbstractC2Reducer):
    """Realization-free sparse reducer using on-demand `NO(d(a), phi)` spans."""

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
            raise ValueError("Could not infer conformal weight for C2 reduction")
        target_weight = _normalize_weight(target_weight)

        candidates: set[Any] = set()
        # C2(V)_w = span{:(∂a)φ: | a ∈ V[h_a], φ ∈ V, h_a + h_φ + 1 = w}.
        # Iterate over ALL basis elements a, not just primary generators, so that
        # derived elements like ∂^n g are correctly included.
        step = _c2_weight_step(self.canonicalizer.generators)
        if self.generator_provider is not None:
            # When an explicit generator_provider is set, fall back to the primary loop
            provider_generators = self.generator_provider
            for generator in provider_generators:
                generator_weight = _get_conformal_weight(generator)
                if generator_weight is None:
                    continue
                phi_weight = sp.simplify(target_weight - generator_weight - 1)
                if phi_weight < 0:
                    continue
                for phi in self.canonicalizer.list(phi_weight):
                    candidate = self.canonicalizer.canonicalize(
                        normal_product(d(generator), phi)
                    )
                    if candidate == Zero:
                        continue
                    candidate_weight = _get_conformal_weight(candidate)
                    if candidate_weight is None or not _weights_equal(
                        candidate_weight, target_weight
                    ):
                        continue
                    candidates.add(candidate)
        else:
            h_a = step
            while sp.simplify(target_weight - h_a - 1) >= 0:
                phi_weight = sp.simplify(target_weight - h_a - 1)
                a_list = self.canonicalizer.list(h_a)
                phi_list = self.canonicalizer.list(phi_weight)
                for a in a_list:
                    for phi in phi_list:
                        candidate = self.canonicalizer.canonicalize(
                            normal_product(d(a), phi)
                        )
                        if candidate == Zero:
                            continue
                        candidate_weight = _get_conformal_weight(candidate)
                        if candidate_weight is None or not _weights_equal(
                            candidate_weight, target_weight
                        ):
                            continue
                        candidates.add(candidate)
                h_a = sp.simplify(h_a + step)

        return self.linear_context.independent_subset(sorted(candidates, key=sp.srepr))

    def quotient_normal_form(self, expr: Any) -> Any:
        return self.solve_c2_witness(expr).remainder

    def solve_c2_witness(self, expr: Any) -> C2ReductionWitness:
        canonical_expr = self.canonicalizer.canonicalize(expr)
        if canonical_expr == Zero:
            return C2ReductionWitness(
                expr=expr,
                remainder=Zero,
                generators=[],
                coefficients=[],
                c2_part=Zero,
                target_weight=None,
            )

        target_weight = _get_conformal_weight(canonical_expr)
        if target_weight is None:
            raise ValueError("Could not infer conformal weight for C2 reduction")
        target_weight = _normalize_weight(target_weight)

        candidates = self.candidate_generators_for_expr(canonical_expr)
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
        # The eliminator normalizes pivot vectors (leading coeff = 1). coeff_map[pivot]
        # is the elimination factor for the normalized pivot, not for the actual
        # candidate (which may have a different leading coefficient at that pivot).
        # Correct alpha_i = coeff_map[pivot_i] / candidate_vectors[i][pivot_i].
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

        c2_part_terms = []
        for coeff, candidate in zip(coeffs, independent_candidates):
            if coeff != 0:
                c2_part_terms.append(coeff * candidate)
        c2_part = self.canonicalizer.canonicalize(
            sp.Add(*c2_part_terms) if c2_part_terms else Zero
        )

        remainder_terms = [
            coeff * monomial for monomial, coeff in reduced_terms.items()
        ]
        remainder = self.canonicalizer.canonicalize(
            sp.Add(*remainder_terms) if remainder_terms else Zero
        )

        return C2ReductionWitness(
            expr=expr,
            remainder=remainder,
            generators=independent_candidates,
            coefficients=coeffs,
            c2_part=c2_part,
            target_weight=target_weight,
        )
