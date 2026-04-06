"""Sparse-first descendant generation utilities."""

from __future__ import annotations

from typing import Any, Iterable

import sympy as sp

from .api import normal_product
from .constants import Zero
from .operators import d
from .operator_spaces import (
    LocalOperatorBasis,
    _get_conformal_weight,
    _normalize_weight,
    _weights_equal,
)


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

        queue = [canonical_source]
        seen = {canonical_source}
        generated = []

        while queue:
            current = queue.pop(0)
            current_weight = _get_conformal_weight(current)
            if current_weight is None or current_weight > normalized_target:
                continue

            if _weights_equal(current_weight, normalized_target):
                generated.append(current)
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

        return sorted(set(generated), key=sp.srepr)

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
        queue: list[Any],
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
