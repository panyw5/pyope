from typing import Any, Dict, Iterable, List, Optional

import sympy as sp

from .operators import BasisOperator, Operator

class RealizedGenerator(BasisOperator):
    @property
    def realization(self) -> Any: ...

def make_realized(expressions: Any, **assumptions: Any) -> list[RealizedGenerator]: ...
def realize(expr: Any) -> Any: ...
def realize_and_simplify(expr: Any) -> Any: ...
def realized_coordinates(
    expr: Any, free_field_basis: "LocalOperatorBasis", weight: Any = None
) -> sp.Matrix: ...
def independent_under_realization(
    expressions: Iterable[Any],
    free_field_basis: "LocalOperatorBasis",
    weight: Any = None,
) -> list[Any]: ...

class LocalOperatorBasis:
    generators: tuple[Operator, ...]
    max_weight: Optional[sp.Expr]
    def __init__(
        self, generators: Iterable[Operator], max_weight: Any = None
    ) -> None: ...
    def canonicalize(self, expr: Any) -> Any: ...
    def enumerate_candidates(self, weight: Any) -> list[Any]: ...
    def basis(self, weight: Any) -> list[Any]: ...
    def coordinates(self, expr: Any, weight: Any = None) -> sp.Matrix: ...
    def realized_coordinates(
        self, expr: Any, free_field_basis: "LocalOperatorBasis", weight: Any = None
    ) -> sp.Matrix: ...
    def independent_under_realization(
        self,
        expressions: Iterable[Any],
        free_field_basis: "LocalOperatorBasis",
        weight: Any = None,
    ) -> list[Any]: ...

class DescendantSpace:
    basis_builder: LocalOperatorBasis
    def __init__(self, basis_builder: LocalOperatorBasis) -> None: ...
    def generate(self, source: Any, target_weight: Any) -> List[Any]: ...
    def basis(self, source: Any, target_weight: Any) -> List[Any]: ...
    def span(self, sources: Iterable[Any], target_weight: Any) -> List[Any]: ...

class SingularVectorAnalyzer:
    basis_builder: LocalOperatorBasis
    generators: tuple[Any, ...]
    stress_tensor: Any
    def __init__(
        self,
        basis_builder: LocalOperatorBasis,
        generators: Optional[Iterable[Any]] = None,
        stress_tensor: Any = None,
    ) -> None: ...
    def positive_mode_constraints(
        self, expr: Any, generators: Optional[Iterable[Any]] = None
    ) -> Dict[Any, Dict[int, Any]]: ...
    def is_singular(
        self, expr: Any, generators: Optional[Iterable[Any]] = None
    ) -> bool: ...
    def find_singular_vectors(
        self,
        weight: Any,
        ansatz: Optional[Iterable[Any]] = None,
        generators: Optional[Iterable[Any]] = None,
    ) -> List[Dict[str, Any]]: ...

class C2Space:
    basis_builder: LocalOperatorBasis
    def __init__(self, basis_builder: LocalOperatorBasis) -> None: ...
    def generators(self, weight: Any) -> List[Any]: ...
    def basis(self, weight: Any) -> List[Any]: ...
    def contains(self, expr: Any, weight: Any = None) -> bool: ...

class C2NullSearcher:
    basis_builder: LocalOperatorBasis
    descendants: DescendantSpace
    c2_space: C2Space
    stress_tensor: Any
    def __init__(
        self,
        basis_builder: LocalOperatorBasis,
        descendants: Optional[DescendantSpace] = None,
        c2_space: Optional[C2Space] = None,
        stress_tensor: Any = None,
    ) -> None: ...
    def search_from_sources(
        self, target_weight: Any, sources: Iterable[Any], target_expr: Any
    ) -> Optional[Dict[str, Any]]: ...
    def search_stress_tensor_nilpotency(
        self, n: int, sources: Iterable[Any], stress_tensor: Any = None
    ) -> Optional[Dict[str, Any]]: ...
