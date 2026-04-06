from sympy import Rational

from pyope import BasisOperator, Bosonic, C2Space, LocalOperatorBasis, NO, Zero, d


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


def test_c2_space_exposes_reducer_style_quotient_api():
    T = BasisOperator("T_c2_quotient_api", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=5)
    c2 = C2Space(basis_builder)

    assert c2.quotient_normal_form(NO(d(T), T), weight=5) == Zero
    assert c2.is_zero_mod_c2(NO(d(T), T), weight=5)

    witness = c2.solve_c2_witness(NO(d(T), T), weight=5)

    assert witness.remainder == Zero
    assert witness.c2_part == NO(d(T), T)


def test_c2_space_second_derivative_of_generator_is_in_c2():
    """d(d(T)) must be in C2 at weight 4.

    This requires sweeping *derived* basis elements (here ∂T at weight 3) as
    the 'a' factor in :(∂a)φ:.  The old code only used primary generators and
    therefore missed this case.
    """
    T = BasisOperator("T_c2_ddT", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=4)
    c2 = C2Space(basis_builder)

    # d(d(T)) = ∂²T has weight 4; it equals :(∂(∂T))·1: so it is in C2
    assert c2.contains(d(d(T)), weight=4)
    # NO(T,T) is NOT in C2[4]
    assert not c2.contains(NO(T, T), weight=4)


def test_c2_space_half_integer_weight_derived_element_is_in_c2():
    """For a weight-3/2 generator G, d(d(G)) must be in C2 at weight 7/2.

    This test verifies the _c2_weight_step logic (step = 1/2) for algebras
    with half-integer conformal weights.  d(G) has weight 5/2 and acts as 'a'
    in :(∂(∂G))·1:; the old primary-only loop only tried G (weight 3/2) and
    landed on φ-weight 1 where list() is empty, so it missed d(d(G)).
    """
    G = BasisOperator("G_c2_half", conformal_weight=Rational(3, 2))
    Bosonic(G)

    basis_builder = LocalOperatorBasis([G], max_weight=Rational(7, 2))
    c2 = C2Space(basis_builder)

    # d(d(G)) = ∂²G has weight 7/2; it equals :(∂(∂G))·1: so it is in C2
    assert c2.contains(d(d(G)), weight=Rational(7, 2))
    # d(G) has weight 5/2; it equals :(∂G)·1: so it also is in C2
    assert c2.contains(d(G), weight=Rational(5, 2))
