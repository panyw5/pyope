from pyope import BasicOperator, Bosonic, DescendantSpace, LocalOperatorBasis, NO, d


def test_descendant_space_generates_virasoro_weight_four_descendants():
    T = BasicOperator("T_desc_gen", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T])
    descendants = DescendantSpace(basis_builder)

    generated = descendants.generate(T, 4)

    assert d(T, 2) in generated
    assert NO(T, T) in generated


def test_descendant_space_basis_returns_independent_span():
    T = BasicOperator("T_desc_basis", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T])
    descendants = DescendantSpace(basis_builder)

    basis = descendants.basis(T, 4)

    assert len(basis) == 2
    assert d(T, 2) in basis
    assert NO(T, T) in basis


def test_descendant_space_span_combines_multiple_sources():
    T = BasicOperator("T_desc_span", conformal_weight=2)
    J = BasicOperator("J_desc_span", conformal_weight=1)
    Bosonic(T, J)

    basis_builder = LocalOperatorBasis([T, J])
    descendants = DescendantSpace(basis_builder)

    span = descendants.span([J], 3)

    assert d(J, 2) in span
    assert NO(T, J) in span


def test_descendant_space_returns_empty_below_source_weight():
    T = BasicOperator("T_desc_empty", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T])
    descendants = DescendantSpace(basis_builder)

    assert descendants.generate(T, 1) == []
