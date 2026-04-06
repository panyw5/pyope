from typing import Any, Dict, Iterable, List, Optional

import sympy as sp

from .operators import BasisOperator, Operator

class RealizedGenerator(BasisOperator):
    @property
    def realization(self) -> Any: ...

def make_realized(expressions: Any, **assumptions: Any) -> list[RealizedGenerator]: ...
def realize(expr: Any) -> Any: ...
def realize_and_simplify(expr: Any) -> Any: ...
def clear_realize_cache() -> None: ...
def realized_coordinates(
    expr: Any, free_field_basis: "LocalOperatorBasis", weight: Any = None
) -> sp.Matrix: ...
def independent_under_realization(
    expressions: Iterable[Any],
    free_field_basis: "LocalOperatorBasis",
    weight: Any = None,
) -> list[Any]: ...
def list_independent_ops(
    expressions: Iterable[Any],
    basis: "LocalOperatorBasis",
    weight: Any = None,
) -> list[Any]: ...
def list_zero_relations(
    expressions: Iterable[Any],
    basis: "LocalOperatorBasis",
    weight: Any = None,
) -> list[dict[str, Any]]: ...

class SparseLinearContext:
    def __init__(self) -> None: ...

class LocalOperatorBasis:
    generators: tuple[Operator, ...]
    stress_tensor: Any
    gradings: Any
    max_weight: Optional[sp.Expr]
    def __init__(
        self,
        generators: Iterable[Operator],
        stress_tensor: Any = None,
        gradings: Any = None,
        max_weight: Any = None,
        max_occurence: Any = None,
    ) -> None: ...
    def canonicalize(self, expr: Any) -> Any: ...
    def sparse_terms(self, expr: Any) -> dict[Any, sp.Expr]: ...
    def basis(self, weight: Any, max_occurence: Any = None) -> list[Any]: ...
    def enumerate_candidates(self, weight: Any) -> list[Any]: ...
    def list(self, weight: Any, max_occurence: Any = None) -> list[Any]: ...
    def sector_of(self, expr: Any) -> dict[str, Any]: ...
    def nested_stress_tensor(self, n: int) -> Any: ...
    def coordinates(
        self, expr: Any, weight: Any = None, max_occurence: Any = None
    ) -> sp.Matrix: ...
    def realized_coordinates(
        self,
        expr: Any,
        free_field_basis: "LocalOperatorBasis",
        weight: Any = None,
        max_occurence: Any = None,
    ) -> sp.Matrix: ...
    def independent_under_realization(
        self,
        expressions: Iterable[Any],
        free_field_basis: "LocalOperatorBasis",
        weight: Any = None,
        max_occurence: Any = None,
    ) -> list[Any]: ...
    def list_independent_ops(
        self,
        expressions: Iterable[Any],
        weight: Any = None,
        max_occurence: Any = None,
    ) -> list[Any]: ...
    def list_zero_relations(
        self,
        expressions: Iterable[Any],
        weight: Any = None,
        max_occurence: Any = None,
    ) -> list[dict[str, Any]]: ...

class LocalOperatorCanonicalizer(LocalOperatorBasis): ...

class C2Space:
    basis_builder: LocalOperatorBasis
    def __init__(self, basis_builder: LocalOperatorBasis) -> None: ...
    def generators(self, weight: Any) -> List[Any]: ...
    def basis(self, weight: Any) -> List[Any]: ...
    def contains(self, expr: Any, weight: Any = None) -> bool: ...

class C2NullSearcher:
    basis_builder: LocalOperatorBasis
    def __init__(
        self,
        basis_builder: LocalOperatorBasis,
        stress_tensor: Any = None,
    ) -> None: ...
    def search_from_sources(
        self, target_weight: Any, sources: Iterable[Any], target_expr: Any
    ) -> Optional[Dict[str, Any]]: ...
