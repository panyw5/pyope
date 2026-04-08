import pytest
import sympy as sp

import pyope.operator_spaces as operator_spaces_module
from pyope import (
    BasisOperator,
    Bosonic,
    compute_backend,
    independent_under_realization,
    LocalOperatorBasis,
    list_independent_ops,
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

    basis = LocalOperatorBasis([T])

    weight_four_basis = basis.list(4)

    assert d(T, 2) in weight_four_basis
    assert NO(T, T) in weight_four_basis
    assert len(weight_four_basis) == 2


def test_local_operator_basis_canonicalizes_and_extracts_coordinates():
    T = BasisOperator("T_basis_coord", conformal_weight=2)
    J = BasisOperator("J_basis_coord", conformal_weight=1)
    Bosonic(T, J)

    basis = LocalOperatorBasis([T, J])
    expr = 2 * NO(T, J) + 3 * NO(T, J) + d(T)

    canonical = basis.canonicalize(expr)
    coordinates = basis.coordinates(expr, weight=3)
    ordered_basis = basis.list(3)
    index = {op: i for i, op in enumerate(ordered_basis)}

    assert canonical == 5 * NO(T, J) + d(T)
    assert coordinates[index[NO(T, J)], 0] == 5
    assert coordinates[index[d(T)], 0] == 1


def test_local_operator_basis_canonicalize_cache_tracks_registry_changes():
    A = BasisOperator("A_basis_cache_version", conformal_weight=1)
    B = BasisOperator("B_basis_cache_version", conformal_weight=1)

    basis = LocalOperatorBasis([A, B])
    expr = NO(B, A)

    assert basis.canonicalize(expr) == NO(A, B)

    Bosonic(B, A)

    assert basis.canonicalize(expr) == NO(B, A)


def test_local_operator_basis_sorts_generators_before_enumeration():
    A = BasisOperator("A_basis_sorted", conformal_weight=1)
    B = BasisOperator("B_basis_sorted", conformal_weight=1)
    Bosonic(A, B)

    basis = LocalOperatorBasis([B, A])

    assert basis.generators == (A, B)
    assert basis.list(2) == [d(A), NO(A, A), NO(A, B), d(B), NO(B, B)]


def test_local_operator_basis_builds_canonical_right_nested_basis():
    A = BasisOperator("A_basis_nested", conformal_weight=1)
    B = BasisOperator("B_basis_nested", conformal_weight=1)
    C = BasisOperator("C_basis_nested", conformal_weight=1)
    Bosonic(A, B, C)

    basis = LocalOperatorBasis([C, A, B])

    assert NO(A, NO(B, C)) in basis.list(3)
    assert NO(C, NO(A, B)) not in basis.list(3)


def test_local_operator_basis_list_is_backend_independent():
    A = BasisOperator("A_basis_backend", conformal_weight=1)
    B = BasisOperator("B_basis_backend", conformal_weight=1)
    Bosonic(A, B)

    basis = LocalOperatorBasis([B, A])
    sympy_basis = basis.list(2)

    with compute_backend("wolfram"):
        wolfram_basis = basis.list(2)

    assert wolfram_basis == sympy_basis


def test_local_operator_basis_weight_zero_is_vacuum():
    T = BasisOperator("T_basis_vac", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T])

    assert basis.list(0) == [One]


def test_local_operator_basis_rejects_weight_mismatch():
    T = BasisOperator("T_basis_mismatch", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T])

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

    basis = LocalOperatorBasis([T])
    weight_four_basis = basis.list(4)

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


def test_realize_and_simplify_uses_wolfram_precanonicalization_when_available(
    monkeypatch,
):
    J = BasisOperator("J_realize_simplify_wolfram_pre", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize_simplify_wolfram_pre", realization=NO(J, J))

    monkeypatch.setattr(
        operator_spaces_module,
        "_should_use_wolfram_precanonicalization",
        lambda: True,
    )

    calls: list[list[object]] = []

    def fake_batch(exprs):
        calls.append(list(exprs))
        return [2 * NO(d(J), J)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_batch_wolfram_precanonicalize",
        fake_batch,
    )

    realized = realize_and_simplify(d(W))

    assert calls == [[d(NO(J, J))]]
    assert realized == 2 * NO(d(J), J)


def test_batch_wolfram_precanonicalize_uses_simplify_expr_not_eval(monkeypatch):
    A = BasisOperator("A_batch_precanonicalize_route", conformal_weight=1)
    B = BasisOperator("B_batch_precanonicalize_route", conformal_weight=1)
    Bosonic(A, B)

    simplify_calls: list[list[object]] = []
    eval_calls: list[object] = []

    monkeypatch.setattr(
        operator_spaces_module,
        "_should_use_wolfram_precanonicalization",
        lambda: True,
    )

    def fake_simplify_expr(value):
        simplify_calls.append(list(value))
        return [NO(A, B), NO(A, B)]

    def fake_evaluate_expr(value):
        eval_calls.append(value)
        raise AssertionError("evaluate_expr should not be used for precanonicalization")

    import pyope.wolfram_backend as wolfram_backend_module

    monkeypatch.setattr(wolfram_backend_module, "simplify_expr", fake_simplify_expr)
    monkeypatch.setattr(wolfram_backend_module, "evaluate_expr", fake_evaluate_expr)

    result = operator_spaces_module._batch_wolfram_precanonicalize([NO(B, A), NO(A, B)])

    assert simplify_calls == [[NO(B, A), NO(A, B)]]
    assert eval_calls == []
    assert result == [NO(A, B), NO(A, B)]


def test_local_operator_basis_canonicalize_uses_wolfram_precanonicalization_when_available(
    monkeypatch,
):
    A = BasisOperator("A_basis_canonicalize_wolfram_pre", conformal_weight=1)
    B = BasisOperator("B_basis_canonicalize_wolfram_pre", conformal_weight=1)
    Bosonic(A, B)

    basis = LocalOperatorBasis([A, B])

    monkeypatch.setattr(
        operator_spaces_module,
        "_should_use_wolfram_precanonicalization",
        lambda: True,
    )

    calls: list[list[object]] = []

    def fake_batch(exprs):
        calls.append(list(exprs))
        return [NO(A, B)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_batch_wolfram_precanonicalize",
        fake_batch,
    )

    import pyope.wolfram_backend as wolfram_backend_module

    monkeypatch.setattr(
        wolfram_backend_module,
        "evaluate_expr",
        lambda value: (_ for _ in ()).throw(
            AssertionError(
                "evaluate_expr should not be called after Wolfram precanonicalization"
            )
        ),
    )

    canonical = basis.canonicalize(NO(B, A))

    assert calls == [[NO(B, A)]]
    assert canonical == NO(A, B)


def test_local_operator_basis_coordinates_uses_wolfram_precanonicalization_when_available(
    monkeypatch,
):
    A = BasisOperator("A_basis_coordinates_wolfram_pre", conformal_weight=1)
    B = BasisOperator("B_basis_coordinates_wolfram_pre", conformal_weight=1)
    Bosonic(A, B)

    basis = LocalOperatorBasis([A, B])

    monkeypatch.setattr(
        operator_spaces_module,
        "_should_use_wolfram_precanonicalization",
        lambda: True,
    )

    calls: list[list[object]] = []

    def fake_batch(exprs):
        calls.append(list(exprs))
        return [NO(A, B)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_batch_wolfram_precanonicalize",
        fake_batch,
    )

    import pyope.wolfram_backend as wolfram_backend_module

    monkeypatch.setattr(
        wolfram_backend_module,
        "evaluate_expr",
        lambda value: (_ for _ in ()).throw(
            AssertionError(
                "evaluate_expr should not be called after Wolfram precanonicalization"
            )
        ),
    )

    coordinates = basis.coordinates(NO(B, A), weight=2)
    ordered_basis = basis.list(2)
    index = {op: i for i, op in enumerate(ordered_basis)}

    assert calls == [[NO(B, A)]]
    assert coordinates[index[NO(A, B)], 0] == 1


def test_list_independent_ops_filters_dependent_basis_elements():
    J = BasisOperator("J_list_indep", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    expressions = [NO(J, J), 2 * NO(J, J), d(J)]

    independent = list_independent_ops(expressions, basis, weight=2)

    assert len(independent) == 2
    assert d(J) in independent
    assert any(
        basis.canonicalize(expr) in {NO(J, J), 2 * NO(J, J)} for expr in independent
    )


def test_list_independent_ops_does_not_require_weight_for_homogeneous_inputs():
    J = BasisOperator("J_list_indep_no_weight", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    expressions = [NO(J, J), 2 * NO(J, J), d(J)]

    independent = list_independent_ops(expressions, basis)

    assert len(independent) == 2
    assert d(J) in independent


def test_list_independent_ops_avoids_basis_enumeration():
    J = BasisOperator("J_list_indep_local", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])

    def fail_basis(*args, **kwargs):
        raise AssertionError("basis() should not be called for local dependence checks")

    basis.basis = fail_basis  # type: ignore[method-assign]
    expressions = [NO(J, J), 2 * NO(J, J), d(J)]

    independent = basis.list_independent_ops(expressions)

    assert len(independent) == 2


def test_list_independent_ops_handles_incremental_dependence_chain():
    J = BasisOperator("J_list_indep_chain", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    jj = NO(J, J)
    dj = d(J)

    expressions = [jj, dj, jj + dj, 3 * jj - 2 * dj]

    independent = basis.list_independent_ops(expressions)

    assert len(independent) == 2
    assert all(expr in expressions for expr in independent)
    assert basis.list_zero_relations(independent) == []


def test_local_operator_basis_list_independent_ops_uses_default_max_occurence():
    phi = BasisOperator("phi_list_indep_zero", conformal_weight=0)
    Bosonic(phi)

    basis = LocalOperatorBasis([phi], max_occurence=1)
    expressions = [d(phi), NO(phi, d(phi)), 2 * NO(phi, d(phi))]

    independent = basis.list_independent_ops(expressions, weight=1)
    ordered_basis = basis.list(1)
    target_index = ordered_basis.index(basis.canonicalize(NO(phi, d(phi))))

    assert len(independent) == 2
    assert d(phi) in independent
    assert any(
        basis.coordinates(expr, weight=1)[target_index, 0] != 0 for expr in independent
    )


def test_list_zero_relations_finds_dependent_linear_combination():
    J = BasisOperator("J_zero_relation", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    jj = NO(J, J)
    expressions = [jj, 2 * jj, d(J)]

    relations = list_zero_relations(expressions, basis, weight=2)

    assert len(relations) == 1
    relation = relations[0]

    assert relation["operators"] == expressions
    assert relation["coefficients"] == [-2, 1, 0]
    assert relation["terms"] == [(jj, -2), (2 * jj, 1)]
    assert basis.canonicalize(relation["relation"]) == 0


def test_list_zero_relations_does_not_require_weight_for_homogeneous_inputs():
    J = BasisOperator("J_zero_relation_no_weight", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    jj = NO(J, J)
    expressions = [jj, 2 * jj, d(J)]

    relations = list_zero_relations(expressions, basis)

    assert len(relations) == 1
    assert relations[0]["coefficients"] == [-2, 1, 0]


def test_list_zero_relations_avoids_basis_enumeration():
    J = BasisOperator("J_zero_relation_local", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])

    def fail_basis(*args, **kwargs):
        raise AssertionError("basis() should not be called for local relation checks")

    basis.basis = fail_basis  # type: ignore[method-assign]
    jj = NO(J, J)
    relations = basis.list_zero_relations([jj, 2 * jj, d(J)])

    assert len(relations) == 1
    assert relations[0]["coefficients"] == [-2, 1, 0]


def test_list_zero_relations_handles_dependence_chain_with_compressed_coordinates():
    J = BasisOperator("J_zero_relation_chain", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    jj = NO(J, J)
    dj = d(J)
    expressions = [jj, dj, jj + dj, 3 * jj - 2 * dj]

    relations = basis.list_zero_relations(expressions)

    assert len(relations) == 2
    for relation in relations:
        assert relation["operators"] == expressions
        assert basis.canonicalize(relation["relation"]) == 0


def test_list_independent_ops_uses_wolfram_precanonicalization_when_available(
    monkeypatch,
):
    J = BasisOperator("J_list_indep_wolfram_pre", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    expressions = [2 * NO(J, J), NO(J, J), d(J)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_should_use_wolfram_precanonicalization",
        lambda: True,
    )

    calls: list[list[object]] = []

    def fake_batch(exprs):
        calls.append(list(exprs))
        return [NO(J, J), NO(J, J), d(J)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_batch_wolfram_precanonicalize",
        fake_batch,
    )

    independent = list_independent_ops(expressions, basis, weight=2)

    assert calls == [[d(J), 2 * NO(J, J), NO(J, J)]]
    assert independent == [d(J), NO(J, J)]


def test_list_zero_relations_uses_wolfram_precanonicalization_when_available(
    monkeypatch,
):
    J = BasisOperator("J_zero_relation_wolfram_pre", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    jj = NO(J, J)
    expressions = [jj, 2 * jj, d(J)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_should_use_wolfram_precanonicalization",
        lambda: True,
    )

    calls: list[list[object]] = []

    def fake_batch(exprs):
        calls.append(list(exprs))
        return [jj, jj, d(J)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_batch_wolfram_precanonicalize",
        fake_batch,
    )

    relations = list_zero_relations(expressions, basis, weight=2)

    assert calls == [[jj, 2 * jj, d(J)]]
    assert len(relations) == 1
    assert relations[0]["coefficients"] == [-1, 1, 0]


def test_independent_under_realization_uses_wolfram_precanonicalization_when_available(
    monkeypatch,
):
    J = BasisOperator("J_realize_indep_wolfram_pre", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize_indep_wolfram_pre", realization=NO(J, J))
    free_field_basis = LocalOperatorBasis([J])
    expressions = [W, NO(J, J)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_should_use_wolfram_precanonicalization",
        lambda: True,
    )

    calls: list[list[object]] = []

    def fake_batch(exprs):
        calls.append(list(exprs))
        return [NO(J, J), NO(J, J)]

    monkeypatch.setattr(
        operator_spaces_module,
        "_batch_wolfram_precanonicalize",
        fake_batch,
    )

    independent = independent_under_realization(
        expressions,
        free_field_basis=free_field_basis,
        weight=2,
    )

    assert calls == [[NO(J, J), NO(J, J)]]
    assert len(independent) == 1


def test_sparse_terms_distributes_symbolic_scalars_over_operator_sums():
    b = BasisOperator("b_sparse_expand", fermionic=True, conformal_weight=sp.Integer(2))
    c = BasisOperator(
        "c_sparse_expand", fermionic=True, conformal_weight=sp.Integer(-1)
    )
    Bosonic(b, c)

    canonicalizer = LocalOperatorBasis([b, c], max_occurence=4)
    x = NO(b, c)
    y = NO(b, d(c))
    a0, a1 = sp.symbols("a0:2")

    expr = 5 * x + y - a0 * (x + y) - a1 * (x - y)
    terms = canonicalizer.sparse_terms(expr)

    assert terms == {
        x: sp.expand(5 - a0 - a1),
        y: sp.expand(1 - a0 + a1),
    }


def test_local_operator_basis_list_zero_relations_uses_default_max_occurence():
    phi = BasisOperator("phi_zero_relation", conformal_weight=0)
    Bosonic(phi)

    basis = LocalOperatorBasis([phi], max_occurence=1)
    expr = NO(phi, d(phi))
    expressions = [d(phi), expr, 2 * expr]

    relations = basis.list_zero_relations(expressions, weight=1)

    assert len(relations) == 1
    relation = relations[0]

    assert relation["operators"] == expressions
    assert relation["coefficients"] == [0, -2, 1]
    assert basis.canonicalize(relation["relation"]) == 0


def test_independent_under_realization_filters_dependent_abstract_basis():
    J = BasisOperator("J_realize_indep", conformal_weight=1)
    Bosonic(J)

    W = RealizedGenerator("W_realize_indep", realization=NO(J, J))
    free_field_basis = LocalOperatorBasis([J])

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

    basis = LocalOperatorBasis([c], max_occurence=0)
    weight_one_basis = basis.list(1)

    assert d(c, 2) in weight_one_basis
    assert len(weight_one_basis) > 0


def test_local_operator_basis_requires_max_occurence_for_nonpositive_generators():
    phi = BasisOperator("phi_basis_nonpositive_boson", conformal_weight=0)
    Bosonic(phi)

    with pytest.raises(ValueError) as exc_info:
        LocalOperatorBasis([phi])

    assert "max_occurence" in str(exc_info.value)


def test_local_operator_basis_requires_max_occurence_for_negative_fermionic_generators():
    c = BasisOperator("c_basis_requires_cutoff", fermionic=True, conformal_weight=-1)

    with pytest.raises(ValueError) as exc_info:
        LocalOperatorBasis([c])

    assert "max_occurence" in str(exc_info.value)


def test_local_operator_basis_truncates_weight_zero_bosonic_atoms_with_max_occurence():
    phi = BasisOperator("phi_basis_zero_boson", conformal_weight=0)
    Bosonic(phi)

    basis = LocalOperatorBasis([phi], max_occurence=2)

    weight_zero_basis = basis.list(0, max_occurence=2)
    weight_one_basis = basis.list(1, max_occurence=1)

    assert set(weight_zero_basis) == {One, phi, NO(phi, phi)}
    assert len(weight_zero_basis) == 3
    assert d(phi) in weight_one_basis
    assert basis.canonicalize(NO(phi, d(phi))) in weight_one_basis
    assert len(weight_one_basis) == 2


def test_local_operator_basis_zero_max_occurence_keeps_only_vacuum_sector():
    phi = BasisOperator("phi_basis_zero_cutoff", conformal_weight=0)
    Bosonic(phi)

    basis = LocalOperatorBasis([phi], max_occurence=0)

    assert basis.list(0) == [One]
    assert basis.list(1) == [d(phi)]


def test_local_operator_basis_rejects_invalid_max_occurence():
    phi = BasisOperator("phi_basis_bad_cutoff", conformal_weight=0)
    Bosonic(phi)

    basis = LocalOperatorBasis([phi], max_occurence=0)

    with pytest.raises(ValueError):
        basis.list(0, max_occurence=-1)

    with pytest.raises(ValueError):
        basis.list(0, max_occurence=sp.Rational(1, 2))


def test_coordinates_supports_max_occurence_for_weight_zero_bosonic_atoms():
    phi = BasisOperator("phi_basis_coords_zero", conformal_weight=0)
    Bosonic(phi)

    basis = LocalOperatorBasis([phi], max_occurence=1)
    expr = basis.canonicalize(NO(phi, d(phi)))

    coords = basis.coordinates(expr, weight=1)
    ordered_basis = basis.list(1)
    index = {op: i for i, op in enumerate(ordered_basis)}

    assert coords[index[expr], 0] == 1


def test_realized_coordinates_supports_max_occurence_for_zero_weight_free_fields():
    beta = BasisOperator("beta_realized_zero_weight", conformal_weight=1)
    gamma = BasisOperator("gamma_realized_zero_weight", conformal_weight=0)
    Bosonic(beta, gamma)

    Jplus = RealizedGenerator("Jplus_realized_zero_weight", realization=beta)
    Bosonic(Jplus)

    abstract_basis = LocalOperatorBasis([Jplus])
    free_field_basis = LocalOperatorBasis([beta, gamma], max_occurence=1)

    expr = NO(Jplus, gamma)

    coords = abstract_basis.realized_coordinates(expr, free_field_basis, weight=1)
    ordered_basis = free_field_basis.list(1)
    canonical_realized = free_field_basis.canonicalize(NO(beta, gamma))
    index = {op: i for i, op in enumerate(ordered_basis)}

    assert coords[index[canonical_realized], 0] == 1


def test_realized_coordinates_work_for_bc_free_field_ambient_basis():
    b = BasisOperator("b_basis_bc_coords", fermionic=True, conformal_weight=2)
    c = BasisOperator("c_basis_bc_coords", fermionic=True, conformal_weight=-1)
    Bosonic(b, c)

    T = RealizedGenerator(
        "T_basis_bc_coords",
        realization=-2 * NO(b, d(c)) - NO(d(b), c),
        conformal_weight=2,
    )

    abstract_basis = LocalOperatorBasis([T])
    free_field_basis = LocalOperatorBasis([b, c], max_occurence=0)

    expr = abstract_basis.list(4)[0]
    coords = realized_coordinates(expr, free_field_basis, weight=4)

    assert coords.shape[1] == 1
    assert coords != sp.zeros(*coords.shape)


def test_local_operator_basis_allows_repeated_negative_bosonic_atoms_up_to_max_occurence():
    b = BasisOperator(
        "b_basis_repeated_negative_boson", fermionic=True, conformal_weight=2
    )
    c = BasisOperator(
        "c_basis_repeated_negative_boson", fermionic=True, conformal_weight=-1
    )
    beta = BasisOperator(
        "beta_basis_repeated_negative_boson", conformal_weight=sp.Rational(3, 2)
    )
    gamma = BasisOperator(
        "gamma_basis_repeated_negative_boson", conformal_weight=sp.Rational(-1, 2)
    )
    Bosonic(beta, gamma)

    basis = LocalOperatorBasis([b, c, beta, gamma], max_occurence=6)
    expr = NO(c, NO(d(beta, 3), NO(beta, NO(gamma, gamma))))

    canonical = basis.canonicalize(expr)

    assert canonical in basis.list(4)
    coords = basis.coordinates(expr, weight=4)
    ordered_basis = basis.list(4)
    index = {op: i for i, op in enumerate(ordered_basis)}
    assert coords[index[canonical], 0] == 1


def test_local_operator_basis_still_excludes_repeated_fermionic_single_letter_atoms():
    c = BasisOperator("c_basis_repeated_fermion", fermionic=True, conformal_weight=-1)

    basis = LocalOperatorBasis([c], max_occurence=3)

    assert basis.list(0) == [One]
    assert basis.canonicalize(NO(c, c)) == 0
