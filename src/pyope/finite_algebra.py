"""Finite-dimensional associative algebra helpers for research workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp


def _normalize_basis_key(name: Any) -> str:
    return str(name)


@dataclass(frozen=True)
class AlgebraElement:
    """Coordinate vector in a fixed ordered basis."""

    algebra: FiniteDimensionalAlgebra
    coordinates: tuple[sp.Expr, ...]

    def as_expr(self) -> sp.Expr:
        terms = []
        for coeff, basis_name in zip(self.coordinates, self.algebra.basis_names):
            if coeff != 0:
                terms.append(coeff * sp.Symbol(_normalize_basis_key(basis_name)))
        return sp.Add(*terms) if terms else sp.Integer(0)

    def is_zero(self) -> bool:
        return all(sp.simplify(coord) == 0 for coord in self.coordinates)

    def equals(self, other: Any) -> bool:
        other_element = self.algebra.coerce_element(other)
        return all(
            sp.simplify(left - right) == 0
            for left, right in zip(self.coordinates, other_element.coordinates)
        )

    def is_scalar_multiple_of(self, other: Any) -> tuple[bool, sp.Expr | None]:
        other_element = self.algebra.coerce_element(other)
        scalar: sp.Expr | None = None
        for left, right in zip(self.coordinates, other_element.coordinates):
            left = sp.sympify(left)
            right = sp.sympify(right)
            if sp.simplify(right) == 0:
                if sp.simplify(left) != 0:
                    return False, None
                continue
            candidate = sp.simplify(left / right)
            if scalar is None:
                scalar = candidate
                continue
            if sp.simplify(candidate - scalar) != 0:
                return False, None

        if scalar is None:
            if self.is_zero() and other_element.is_zero():
                return True, sp.Integer(0)
            return False, None
        return True, scalar

    def __repr__(self) -> str:
        return f"AlgebraElement({sp.sstr(self.as_expr())})"

    def __add__(self, other: Any) -> AlgebraElement:
        other_element = self.algebra.coerce_element(other)
        return self.algebra.element(
            [
                left + right
                for left, right in zip(self.coordinates, other_element.coordinates)
            ]
        )

    def __sub__(self, other: Any) -> AlgebraElement:
        other_element = self.algebra.coerce_element(other)
        return self.algebra.element(
            [
                left - right
                for left, right in zip(self.coordinates, other_element.coordinates)
            ]
        )

    def __rmul__(self, scalar: Any) -> AlgebraElement:
        scalar = sp.sympify(scalar)
        return self.algebra.element([scalar * coord for coord in self.coordinates])

    def __mul__(self, other: Any) -> AlgebraElement:
        return self.algebra.multiply(self, other)


@dataclass(frozen=True)
class MultiplicationTableEntry:
    left: Any
    right: Any
    product: AlgebraElement


@dataclass(frozen=True)
class OneDimensionalRepresentation:
    values: dict[Any, sp.Expr]

    def evaluate(self, element: AlgebraElement) -> sp.Expr:
        total = sp.Integer(0)
        for coeff, basis_name in zip(element.coordinates, element.algebra.basis_names):
            total += sp.sympify(coeff) * sp.sympify(self.values[basis_name])
        return sp.simplify(total)


@dataclass(frozen=True)
class IrreducibleRepresentationClassification:
    dimension: int | None
    irreducible_representations: list[Any]
    method: str
    notes: str


@dataclass(frozen=True)
class AlgebraValidationIssue:
    kind: str
    message: str
    left: Any | None = None
    middle: Any | None = None
    right: Any | None = None
    expected: AlgebraElement | None = None
    actual: AlgebraElement | None = None


@dataclass(frozen=True)
class AssociativityValidationResult:
    is_associative: bool
    issues: list[AlgebraValidationIssue]


@dataclass(frozen=True)
class IdentityValidationResult:
    is_identity_valid: bool
    issues: list[AlgebraValidationIssue]


class FiniteDimensionalAlgebra:
    """Finite-dimensional associative algebra described by structure constants."""

    def __init__(
        self,
        basis_names: Iterable[Any],
        structure_constants: dict[tuple[Any, Any], dict[Any, Any]],
        identity: Any | None = None,
    ):
        ordered_basis = tuple(basis_names)
        if not ordered_basis:
            raise ValueError("basis_names must not be empty")

        self.basis_names = ordered_basis
        self.index = {name: idx for idx, name in enumerate(self.basis_names)}
        if len(self.index) != len(self.basis_names):
            raise ValueError("basis_names must be unique")

        self.identity = identity
        self.structure_constants: dict[tuple[Any, Any], tuple[sp.Expr, ...]] = {}
        zero = tuple(sp.Integer(0) for _ in self.basis_names)
        for left in self.basis_names:
            for right in self.basis_names:
                raw = structure_constants.get((left, right), {})
                coords = [sp.Integer(0) for _ in self.basis_names]
                for basis_name, coeff in raw.items():
                    if basis_name not in self.index:
                        raise ValueError(
                            "Unknown basis element in structure constants: "
                            f"{basis_name}"
                        )
                    coords[self.index[basis_name]] = sp.sympify(coeff)
                self.structure_constants[(left, right)] = tuple(coords) if raw else zero

    @property
    def dimension(self) -> int:
        return len(self.basis_names)

    def basis_element(self, name: Any) -> AlgebraElement:
        if name not in self.index:
            raise ValueError(f"Unknown basis element: {name}")
        coords = [sp.Integer(0) for _ in self.basis_names]
        coords[self.index[name]] = sp.Integer(1)
        return self.element(coords)

    def basis(self) -> list[AlgebraElement]:
        return [self.basis_element(name) for name in self.basis_names]

    def element(self, coordinates: Iterable[Any]) -> AlgebraElement:
        coords = tuple(sp.sympify(value) for value in coordinates)
        if len(coords) != self.dimension:
            raise ValueError("Coordinate vector has wrong length")
        return AlgebraElement(self, coords)

    def coerce_element(self, value: Any) -> AlgebraElement:
        if isinstance(value, AlgebraElement):
            if value.algebra is not self:
                raise ValueError("Element belongs to a different algebra")
            return value
        if isinstance(value, (list, tuple)):
            return self.element(value)
        if value in self.index:
            return self.basis_element(value)
        raise TypeError(f"Cannot coerce {type(value)} to AlgebraElement")

    def multiply(self, left: Any, right: Any) -> AlgebraElement:
        left_element = self.coerce_element(left)
        right_element = self.coerce_element(right)
        coords = [sp.Integer(0) for _ in self.basis_names]
        for i, left_coeff in enumerate(left_element.coordinates):
            if left_coeff == 0:
                continue
            left_name = self.basis_names[i]
            for j, right_coeff in enumerate(right_element.coordinates):
                if right_coeff == 0:
                    continue
                right_name = self.basis_names[j]
                product = self.structure_constants[(left_name, right_name)]
                for k, basis_coeff in enumerate(product):
                    if basis_coeff != 0:
                        coords[k] += left_coeff * right_coeff * basis_coeff
        return self.element(coords)

    def multiplication_table(self) -> list[MultiplicationTableEntry]:
        entries = []
        for left in self.basis_names:
            for right in self.basis_names:
                entries.append(
                    MultiplicationTableEntry(
                        left=left,
                        right=right,
                        product=self.multiply(left, right),
                    )
                )
        return entries

    def multiplication_tensor(self) -> sp.MutableDenseNDimArray:
        tensor = sp.MutableDenseNDimArray.zeros(
            self.dimension, self.dimension, self.dimension
        )
        for i, left in enumerate(self.basis_names):
            for j, right in enumerate(self.basis_names):
                product = self.structure_constants[(left, right)]
                for k, coeff in enumerate(product):
                    tensor[i, j, k] = sp.sympify(coeff)
        return tensor

    def left_multiplication_operator(self, basis_name: Any) -> sp.Matrix:
        return self.left_regular_matrix(basis_name)

    def right_multiplication_operator(self, basis_name: Any) -> sp.Matrix:
        return self.right_regular_matrix(basis_name)

    def left_regular_matrix(self, element: Any) -> sp.Matrix:
        target = self.coerce_element(element)
        columns = []
        for basis_name in self.basis_names:
            product = self.multiply(target, basis_name)
            columns.append(sp.Matrix(product.coordinates))
        return sp.Matrix.hstack(*columns)

    def right_regular_matrix(self, element: Any) -> sp.Matrix:
        target = self.coerce_element(element)
        columns = []
        for basis_name in self.basis_names:
            product = self.multiply(basis_name, target)
            columns.append(sp.Matrix(product.coordinates))
        return sp.Matrix.hstack(*columns)

    def validate_associativity(self) -> AssociativityValidationResult:
        issues = []
        for left in self.basis_names:
            for middle in self.basis_names:
                for right in self.basis_names:
                    lhs = self.multiply(self.multiply(left, middle), right)
                    rhs = self.multiply(left, self.multiply(middle, right))
                    if not lhs.equals(rhs):
                        issues.append(
                            AlgebraValidationIssue(
                                kind="associativity",
                                message=(
                                    "Associativity fails for "
                                    f"({left} * {middle}) * {right} = {lhs.as_expr()} "
                                    "but "
                                    f"{left} * ({middle} * {right}) = {rhs.as_expr()}"
                                ),
                                left=left,
                                middle=middle,
                                right=right,
                                expected=rhs,
                                actual=lhs,
                            )
                        )
        return AssociativityValidationResult(
            is_associative=not issues,
            issues=issues,
        )

    def validate_identity(self) -> IdentityValidationResult:
        if self.identity is None:
            raise ValueError("Cannot validate identity when no identity is declared")
        if self.identity not in self.index:
            raise ValueError(f"Declared identity is not in basis: {self.identity}")

        issues = []
        identity_element = self.basis_element(self.identity)
        for basis_name in self.basis_names:
            basis_element = self.basis_element(basis_name)
            left_product = self.multiply(self.identity, basis_name)
            if not left_product.equals(basis_element):
                issues.append(
                    AlgebraValidationIssue(
                        kind="left-identity",
                        message=(
                            f"Declared identity {self.identity} fails "
                            f"on the left of {basis_name}: got "
                            f"{left_product.as_expr()}, expected "
                            f"{basis_element.as_expr()}"
                        ),
                        left=self.identity,
                        right=basis_name,
                        expected=basis_element,
                        actual=left_product,
                    )
                )
            right_product = self.multiply(basis_name, self.identity)
            if not right_product.equals(basis_element):
                issues.append(
                    AlgebraValidationIssue(
                        kind="right-identity",
                        message=(
                            f"Declared identity {self.identity} fails "
                            f"on the right of {basis_name}: got "
                            f"{right_product.as_expr()}, expected "
                            f"{basis_element.as_expr()}"
                        ),
                        left=basis_name,
                        right=self.identity,
                        expected=basis_element,
                        actual=right_product,
                    )
                )

        if not identity_element.equals(self.identity):
            raise AssertionError("Internal identity coercion failed")

        return IdentityValidationResult(
            is_identity_valid=not issues,
            issues=issues,
        )

    def solve_one_dimensional_representations(
        self, symbol_map: dict[Any, sp.Symbol] | None = None
    ) -> list[OneDimensionalRepresentation]:
        symbols = symbol_map or {
            basis_name: sp.Symbol(f"chi_{_normalize_basis_key(basis_name)}")
            for basis_name in self.basis_names
        }
        equations = []
        for left in self.basis_names:
            for right in self.basis_names:
                lhs = symbols[left] * symbols[right]
                product = self.multiply(left, right)
                rhs = sp.Integer(0)
                for coeff, basis_name in zip(product.coordinates, self.basis_names):
                    rhs += sp.sympify(coeff) * symbols[basis_name]
                equations.append(sp.Eq(lhs, sp.expand(rhs)))

        if self.identity is not None:
            equations.append(sp.Eq(symbols[self.identity], 1))

        unknowns = [symbols[basis_name] for basis_name in self.basis_names]
        solutions = sp.solve(equations, unknowns, dict=True)
        representations = []
        for solution in solutions:
            values = {
                basis_name: sp.simplify(solution[symbols[basis_name]])
                for basis_name in self.basis_names
            }
            if all(value == 0 for value in values.values()):
                continue
            representations.append(OneDimensionalRepresentation(values))
        return representations

    def classify_irreducible_representations(
        self,
    ) -> IrreducibleRepresentationClassification:
        one_dimensional = self.solve_one_dimensional_representations()
        return IrreducibleRepresentationClassification(
            dimension=None,
            irreducible_representations=one_dimensional,
            method="one-dimensional search",
            notes=(
                "Current implementation classifies one-dimensional "
                "irreducible representations "
                "by solving the character equations from the structure constants. "
                "Higher-dimensional irreducible module classification is "
                "not implemented yet."
            ),
        )


def build_finite_dimensional_algebra(
    basis_names: Iterable[Any],
    structure_constants: dict[tuple[Any, Any], dict[Any, Any]],
    identity: Any | None = None,
) -> FiniteDimensionalAlgebra:
    return FiniteDimensionalAlgebra(
        basis_names=basis_names,
        structure_constants=structure_constants,
        identity=identity,
    )
