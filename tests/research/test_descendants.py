from pyope import BasicOperator, Bosonic, DescendantSpace, LocalOperatorBasis, NO, d


def test_descendant_space_generates_virasoro_weight_four_descendants_new_module():
    T = BasicOperator("T_desc_gen_new", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T])
    descendants = DescendantSpace(basis_builder)

    generated = descendants.generate(T, 4)

    assert d(T, 2) in generated
    assert NO(T, T) in generated


def test_descendant_space_span_combines_multiple_sources_new_module():
    T = BasicOperator("T_desc_span_new", conformal_weight=2)
    J = BasicOperator("J_desc_span_new", conformal_weight=1)
    Bosonic(T, J)

    basis_builder = LocalOperatorBasis([T, J])
    descendants = DescendantSpace(basis_builder)

    span = descendants.span([J], 3)

    assert d(J, 2) in span
    assert NO(T, J) in span
