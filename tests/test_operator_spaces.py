import pytest
import sympy as sp

from pyope import (
    BasisOperator,
    Bosonic,
    independent_under_realization,
    LocalOperatorBasis,
    list_zero_relations,
    make_realized,
    NO,
    One,
    RealizedGenerator,
    d,
    realize,
    realize_and_simplify,
    realized_coordinates,
)


def test_local_operator_basis_enumerates_virasoro_weight_four_basis():
    T = BasisOperator("T_basis_w4", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T], max_weight=4)

    weight_four_basis = basis.basis(4)

    assert d(T, 2) in weight_four_basis
    assert NO(T, T) in weight_four_basis
    assert len(weight_four_basis) == 2


def test_local_operator_basis_canonicalizes_and_extracts_coordinates():
    T = BasisOperator("T_basis_coord", conformal_weight=2)
    J = BasisOperator("J_basis_coord", conformal_weight=1)
    Bosonic(T, J)

    basis = LocalOperatorBasis([T, J], max_weight=3)
    expr = 2 * NO(T, J) + 3 * NO(T, J) + d(T)

    canonical = basis.canonicalize(expr)
    coordinates = basis.coordinates(expr, weight=3)
    ordered_basis = basis.basis(3)
    index = {op: i for i, op in enumerate(ordered_basis)}

    assert canonical == 5 * NO(T, J) + d(T)
    assert coordinates[index[NO(T, J)], 0] == 5
    assert coordinates[index[d(T)], 0] == 1


def test_local_operator_basis_weight_zero_is_vacuum():
    T = BasisOperator("T_basis_vac", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T])

    assert basis.basis(0) == [One]


def test_local_operator_basis_rejects_weight_mismatch():
    T = BasisOperator("T_basis_mismatch", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T], max_weight=4)

    try:
        basis.coordinates(NO(T, T), weight=3)
    except ValueError as exc:
        assert "weight" in str(exc).lower()
    else:
        raise AssertionError("Expected coordinates() to reject mismatched weights")


def test_realized_generator_can_be_used_in_local_operator_basis():
    b = BasisOperator("b_basis_realized", fermionic=True, conformal_weight=2)
    c = BasisOperator("c_basis_realized", fermionic=True, conformal_weight=-1)
    Bosonic(b, c)

    T = RealizedGenerator(
        "T_realized",
        realization=-2 * NO(b, d(c)) - NO(d(b), c),
        conformal_weight=2,
    )
    Bosonic(T)

    basis = LocalOperatorBasis([T], max_weight=4)
    weight_four_basis = basis.basis(4)

    assert d(T, 2) in weight_four_basis
    assert NO(T, T) in weight_four_basis


def test_realized_generator_infers_weight_and_realization():
    J = BasisOperator("J_basis_realized_infer", conformal_weight=1)
    Bosonic(J)

    composite = NO(J, J)
    W = RealizedGenerator("W_realized", realization=composite)

    assert W.conformal_weight == 2
    assert W.realization == composite


def test_make_realized_preserves_input_order_and_infers_names():
    b = BasisOperator("b_make_realized", fermionic=True, conformal_weight=2)
    c = BasisOperator("c_make_realized", fermionic=True, conformal_weight=-1)
    beta = BasisOperator("beta_make_realized", conformal_weight=1.5)
    gamma = BasisOperator("gamma_make_realized", conformal_weight=-0.5)
    Bosonic(b, c, beta, gamma)

    J0 = NO(b, c) + 2 * NO(beta, gamma)
    Jminus = (
        -NO(beta, NO(gamma, gamma)) - NO(gamma, NO(b, c)) + sp.Rational(3, 2) * d(gamma)
    )

    J0, Jminus = make_realized([J0, Jminus])

    assert isinstance(J0, RealizedGenerator)
    assert isinstance(Jminus, RealizedGenerator)
    assert J0.name == "J0"
    assert Jminus.name == "Jminus"
    assert J0.realization == NO(b, c) + 2 * NO(beta, gamma)
    assert Jminus.realization == (
        -NO(beta, NO(gamma, gamma)) - NO(gamma, NO(b, c)) + sp.Rational(3, 2) * d(gamma)
    )


def test_make_realized_requires_expression_iterable():
    J = BasisOperator("J_make_realized_single", conformal_weight=1)
    Bosonic(J)

    W = NO(J, J)

    with pytest.raises(TypeError):
        make_realized(W)


def test_make_realized_returns_generators_in_input_order():
    J = BasisOperator("J_make_realized_order", conformal_weight=1)
    Bosonic(J)

    W = NO(J, J)
    dW = d(J)

    first, second = make_realized([W, dW])

    assert first.name == "W"
    assert second.name == "dW"


def test_make_realized_prefers_latest_alias_name():
    beta = BasisOperator("beta_make_realized_alias", conformal_weight=1)
    Bosonic(beta)

    Jplus = beta

    (Jplus_realized,) = make_realized([Jplus])

    assert Jplus_realized.name == "Jplus"
    assert Jplus_realized.realization == beta


def test_realize_expands_realized_generator_recursively():
    J = BasisOperator("J_realize_base", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize", realization=NO(J, J))

    expr = NO(W, d(W))
    realized = realize(expr)

    assert realized == NO(NO(J, J), d(NO(J, J)))


def test_operator_realize_method_expands_realized_generator_recursively():
    J = BasisOperator("J_realize_method_base", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize_method", realization=NO(J, J))

    expr = NO(W, d(W))

    realize_method = getattr(expr, "realize")

    assert realize_method() == NO(NO(J, J), d(NO(J, J)))


def test_linear_combination_realize_method_expands_realized_generator_recursively():
    J = BasisOperator("J_realize_add_base", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize_add", realization=NO(J, J))

    expr = 2 * W + d(W)

    realize_method = getattr(expr, "realize")

    assert realize_method() == 2 * NO(J, J) + d(NO(J, J))


def test_realize_and_simplify_produces_free_field_expression():
    J = BasisOperator("J_realize_simplify", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize_simplify", realization=NO(J, J))

    realized = realize_and_simplify(d(W))

    assert realized == 2 * NO(d(J), J)


def test_list_zero_relations_finds_dependent_linear_combination():
    J = BasisOperator("J_zero_relation", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J], max_weight=2)
    jj = NO(J, J)
    expressions = [jj, 2 * jj, d(J)]

    relations = list_zero_relations(expressions, basis, weight=2)

    assert len(relations) == 1
    relation = relations[0]

    assert relation["operators"] == expressions
    assert relation["coefficients"] == [-2, 1, 0]
    assert relation["terms"] == [(jj, -2), (2 * jj, 1)]
    assert basis.canonicalize(relation["relation"]) == 0


def test_local_operator_basis_list_zero_relations_returns_same_kernel_basis():
    J = BasisOperator("J_zero_relation_method", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J], max_weight=2)
    relations = basis.list_zero_relations([NO(J, J), 2 * NO(J, J), d(J)], weight=2)

    assert len(relations) == 1
    assert relations[0]["coefficients"] == [-2, 1, 0]


def test_independent_under_realization_filters_dependent_abstract_basis():
    J = BasisOperator("J_realize_indep", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize_indep", realization=NO(J, J))
    free_field_basis = LocalOperatorBasis([J], max_weight=2)

    expressions = [W, NO(J, J)]
    independent = independent_under_realization(
        expressions,
        free_field_basis=free_field_basis,
        weight=2,
    )

    assert len(independent) == 1
    assert independent[0] in expressions


def test_local_operator_basis_allows_nonpositive_fermionic_atoms_without_recursion():
    c = BasisOperator("c_basis_nonpositive", fermionic=True, conformal_weight=-1)

    basis = LocalOperatorBasis([c], max_weight=1)
    weight_one_basis = basis.basis(1)

    assert d(c, 2) in weight_one_basis
    assert len(weight_one_basis) > 0


def test_local_operator_basis_rejects_nonpositive_bosonic_atoms():
    phi = BasisOperator("phi_basis_nonpositive_boson", conformal_weight=0)
    Bosonic(phi)

    basis = LocalOperatorBasis([phi], max_weight=1)

    try:
        basis.basis(1)
    except ValueError as exc:
        assert "nonpositive bosonic atomic operators" in str(exc)
    else:
        raise AssertionError(
            "Expected LocalOperatorBasis to reject nonpositive bosonic atoms"
        )


def test_realized_coordinates_work_for_bc_free_field_ambient_basis():
    b = BasisOperator("b_basis_bc_coords", fermionic=True, conformal_weight=2)
    c = BasisOperator("c_basis_bc_coords", fermionic=True, conformal_weight=-1)
    Bosonic(b, c)

    T = RealizedGenerator(
        "T_basis_bc_coords",
        realization=-2 * NO(b, d(c)) - NO(d(b), c),
        conformal_weight=2,
    )

    abstract_basis = LocalOperatorBasis([T], max_weight=4)
    free_field_basis = LocalOperatorBasis([b, c], max_weight=4)

    expr = abstract_basis.basis(4)[0]
    coords = realized_coordinates(expr, free_field_basis, weight=4)

    assert coords.shape[1] == 1
    assert coords != sp.zeros(*coords.shape)
