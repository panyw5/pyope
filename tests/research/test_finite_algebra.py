import sympy as sp

from pyope import (
    FiniteDimensionalAlgebra,
    build_finite_dimensional_algebra,
)


def test_finite_dimensional_algebra_multiplies_from_structure_constants():
    algebra = build_finite_dimensional_algebra(
        basis_names=["1", "x"],
        structure_constants={
            ("1", "1"): {"1": 1},
            ("1", "x"): {"x": 1},
            ("x", "1"): {"x": 1},
            ("x", "x"): {"1": 1, "x": 1},
        },
        identity="1",
    )

    assert isinstance(algebra, FiniteDimensionalAlgebra)
    assert algebra.multiply("1", "x").coordinates == (0, 1)
    assert algebra.multiply("x", "x").coordinates == (1, 1)


def test_finite_algebra_exposes_multiplication_table_and_regular_matrix():
    algebra = build_finite_dimensional_algebra(
        basis_names=["1", "x"],
        structure_constants={
            ("1", "1"): {"1": 1},
            ("1", "x"): {"x": 1},
            ("x", "1"): {"x": 1},
            ("x", "x"): {"1": 1, "x": 1},
        },
        identity="1",
    )

    table = algebra.multiplication_table()
    table_dict = {
        (entry.left, entry.right): entry.product.coordinates for entry in table
    }

    assert table_dict[("x", "x")] == (1, 1)

    left_matrix = algebra.left_regular_matrix("x")
    assert left_matrix == sp.Matrix([[0, 1], [1, 1]])

    right_matrix = algebra.right_regular_matrix("x")
    assert right_matrix == sp.Matrix([[0, 1], [1, 1]])

    tensor = algebra.multiplication_tensor()
    assert tensor[1, 1, 0] == 1
    assert tensor[1, 1, 1] == 1

    assert algebra.left_multiplication_operator("x") == left_matrix
    assert algebra.right_multiplication_operator("x") == right_matrix


def test_algebra_element_helpers_detect_zero_equality_and_scalar_multiple():
    algebra = build_finite_dimensional_algebra(
        basis_names=["1", "x"],
        structure_constants={
            ("1", "1"): {"1": 1},
            ("1", "x"): {"x": 1},
            ("x", "1"): {"x": 1},
            ("x", "x"): {"1": 1, "x": 1},
        },
        identity="1",
    )

    zero = algebra.element([0, 0])
    x = algebra.basis_element("x")
    two_x = algebra.element([0, 2])

    assert zero.is_zero() is True
    assert x.equals([0, 1]) is True
    assert x.equals("x") is True
    assert x.is_scalar_multiple_of(two_x) == (True, sp.Rational(1, 2))
    assert x.is_scalar_multiple_of([1, 1]) == (False, None)


def test_validate_associativity_accepts_associative_example():
    algebra = build_finite_dimensional_algebra(
        basis_names=["1", "eps"],
        structure_constants={
            ("1", "1"): {"1": 1},
            ("1", "eps"): {"eps": 1},
            ("eps", "1"): {"eps": 1},
            ("eps", "eps"): {},
        },
        identity="1",
    )

    result = algebra.validate_associativity()

    assert result.is_associative is True
    assert result.issues == []


def test_validate_associativity_reports_nonassociative_example():
    algebra = build_finite_dimensional_algebra(
        basis_names=["a", "b"],
        structure_constants={
            ("a", "a"): {"a": 1},
            ("a", "b"): {"a": 1},
            ("b", "a"): {},
            ("b", "b"): {"b": 1},
        },
    )

    result = algebra.validate_associativity()

    assert result.is_associative is False
    assert result.issues
    assert result.issues[0].kind == "associativity"
    assert result.issues[0].actual is not None
    assert result.issues[0].expected is not None


def test_validate_identity_accepts_declared_identity():
    algebra = build_finite_dimensional_algebra(
        basis_names=["1", "x"],
        structure_constants={
            ("1", "1"): {"1": 1},
            ("1", "x"): {"x": 1},
            ("x", "1"): {"x": 1},
            ("x", "x"): {"1": 1, "x": 1},
        },
        identity="1",
    )

    result = algebra.validate_identity()

    assert result.is_identity_valid is True
    assert result.issues == []


def test_validate_identity_reports_invalid_declared_identity():
    algebra = build_finite_dimensional_algebra(
        basis_names=["e", "x"],
        structure_constants={
            ("e", "e"): {"e": 1},
            ("e", "x"): {"x": 1},
            ("x", "e"): {},
            ("x", "x"): {},
        },
        identity="e",
    )

    result = algebra.validate_identity()

    assert result.is_identity_valid is False
    assert [issue.kind for issue in result.issues] == ["right-identity"]


def test_validate_identity_requires_declared_identity():
    algebra = build_finite_dimensional_algebra(
        basis_names=["e", "f"],
        structure_constants={
            ("e", "e"): {"e": 1},
            ("e", "f"): {},
            ("f", "e"): {},
            ("f", "f"): {"f": 1},
        },
    )

    try:
        algebra.validate_identity()
    except ValueError as exc:
        assert "no identity is declared" in str(exc)
    else:
        raise AssertionError("Expected validate_identity to reject missing identity")


def test_dual_numbers_have_one_one_dimensional_irreducible_representation():
    algebra = build_finite_dimensional_algebra(
        basis_names=["1", "eps"],
        structure_constants={
            ("1", "1"): {"1": 1},
            ("1", "eps"): {"eps": 1},
            ("eps", "1"): {"eps": 1},
            ("eps", "eps"): {},
        },
        identity="1",
    )

    reps = algebra.solve_one_dimensional_representations()

    assert len(reps) == 1
    assert reps[0].values == {"1": 1, "eps": 0}


def test_semisimple_two_idempotent_algebra_has_two_one_dimensional_representations():
    algebra = build_finite_dimensional_algebra(
        basis_names=["e", "f"],
        structure_constants={
            ("e", "e"): {"e": 1},
            ("e", "f"): {},
            ("f", "e"): {},
            ("f", "f"): {"f": 1},
        },
    )

    reps = algebra.solve_one_dimensional_representations()
    value_sets = [rep.values for rep in reps]

    assert len(reps) == 2
    assert {"e": 0, "f": 1} in value_sets
    assert {"e": 1, "f": 0} in value_sets


def test_irreducible_representation_classification_reports_current_scope():
    algebra = build_finite_dimensional_algebra(
        basis_names=["1", "eps"],
        structure_constants={
            ("1", "1"): {"1": 1},
            ("1", "eps"): {"eps": 1},
            ("eps", "1"): {"eps": 1},
            ("eps", "eps"): {},
        },
        identity="1",
    )

    classification = algebra.classify_irreducible_representations()

    assert classification.method == "one-dimensional search"
    assert len(classification.irreducible_representations) == 1
    assert "Higher-dimensional" in classification.notes
