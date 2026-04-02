"""Realization-aware $C_2$ reducers for explicit free-field embeddings."""

from __future__ import annotations

from typing import Any

from .c2 import AbstractC2Reducer, C2ReductionWitness, GenericC2Reducer
from .realizations import (
    DerivativeKillingRealizationBackend,
    IdentityRealizationBackend,
    RealizationBackend,
)


class FreeFieldC2Reducer(AbstractC2Reducer):
    """Backend-aware reducer shell for future free-field acceleration."""

    def __init__(
        self,
        canonicalizer: Any,
        realization_backend: RealizationBackend | None = None,
        fallback_reducer: AbstractC2Reducer | None = None,
    ):
        self.canonicalizer = canonicalizer
        self.realization_backend = realization_backend or IdentityRealizationBackend(
            canonicalizer
        )
        self.fallback_reducer = fallback_reducer or GenericC2Reducer(canonicalizer)

    def quotient_normal_form(self, expr: Any) -> Any:
        realized = self.realization_backend.realize(expr)
        candidate = self.realization_backend.quotient_normal_form(realized)
        if candidate == realized:
            return self.fallback_reducer.quotient_normal_form(expr)
        return self.canonicalizer.canonicalize(candidate)

    def solve_c2_witness(self, expr: Any) -> C2ReductionWitness:
        remainder = self.quotient_normal_form(expr)
        if (
            remainder == expr
            or self.realization_backend.__class__ is IdentityRealizationBackend
        ):
            return self.fallback_reducer.solve_c2_witness(expr)
        return C2ReductionWitness(
            expr=expr,
            remainder=self.canonicalizer.canonicalize(remainder),
            generators=[],
            coefficients=[],
            c2_part=self.canonicalizer.canonicalize(expr - remainder),
            target_weight=None,
        )


class DerivativeKillingFreeFieldC2Reducer(FreeFieldC2Reducer):
    """Concrete first-step free-field reducer using derivative-killing quotient."""

    def __init__(
        self, canonicalizer: Any, fallback_reducer: AbstractC2Reducer | None = None
    ):
        super().__init__(
            canonicalizer,
            realization_backend=DerivativeKillingRealizationBackend(canonicalizer),
            fallback_reducer=fallback_reducer,
        )
