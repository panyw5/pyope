from pyope import (
    BasisOperator,
    Bosonic,
    LocalOperatorBasis,
    MakeOPE,
    NO,
    OPE,
    SingularVectorAnalyzer,
    d,
)


def test_singular_vector_analyzer_detects_primary_current():
    T = BasisOperator("T_sing_primary", conformal_weight=2)
    J = BasisOperator("J_sing_primary", conformal_weight=1)
    Bosonic(T, J)

    OPE[T, J] = MakeOPE([J, d(J)])

    basis_builder = LocalOperatorBasis([T, J])
    analyzer = SingularVectorAnalyzer(basis_builder, generators=[T], stress_tensor=T)

    constraints = analyzer.positive_mode_constraints(J)

    assert constraints[T] == {}
    assert analyzer.is_singular(J)


def test_singular_vector_analyzer_detects_nonsingular_stress_tensor_descendant():
    T = BasisOperator("T_sing_desc", conformal_weight=2)
    Bosonic(T)

    c = 0
    OPE[T, T] = MakeOPE([c, 0, 2 * T, d(T)])

    basis_builder = LocalOperatorBasis([T])
    analyzer = SingularVectorAnalyzer(basis_builder, generators=[T], stress_tensor=T)

    assert not analyzer.is_singular(d(T))


def test_singular_vector_analyzer_finds_primary_ansatz_solution():
    T = BasisOperator("T_sing_find", conformal_weight=2)
    J = BasisOperator("J_sing_find", conformal_weight=1)
    Bosonic(T, J)

    OPE[T, J] = MakeOPE([J, d(J)])
    OPE[T, T] = MakeOPE([0, 0, 2 * T, d(T)])

    basis_builder = LocalOperatorBasis([T, J])
    analyzer = SingularVectorAnalyzer(basis_builder, generators=[T], stress_tensor=T)

    solutions = analyzer.find_singular_vectors(2, ansatz=[T, NO(J, J)])

    assert any(solution["vector"] != 0 for solution in solutions)
