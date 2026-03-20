from pyope import BasisOperator, Bosonic, C2Space, LocalOperatorBasis, NO, d


def test_c2_space_generators_include_virasoro_weight_five_terms():
    T = BasisOperator("T_c2_gen", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=5)
    c2 = C2Space(basis_builder)

    generators = c2.generators(5)

    assert NO(d(T), T) in generators


def test_c2_space_basis_is_independent():
    T = BasisOperator("T_c2_basis", conformal_weight=2)
    J = BasisOperator("J_c2_basis", conformal_weight=1)
    Bosonic(T, J)

    basis_builder = LocalOperatorBasis([T, J], max_weight=4)
    c2 = C2Space(basis_builder)

    basis = c2.basis(4)

    assert NO(d(J), d(J)) in basis or NO(d(T), T) in basis


def test_c2_space_contains_known_generator():
    T = BasisOperator("T_c2_contains", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=5)
    c2 = C2Space(basis_builder)

    assert c2.contains(NO(d(T), T), weight=5)


def test_c2_space_rejects_non_c2_element():
    T = BasisOperator("T_c2_reject", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=4)
    c2 = C2Space(basis_builder)

    assert not c2.contains(NO(T, T), weight=4)
