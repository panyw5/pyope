import sympy as sp

from pyope import (
    AbstractZhuReducer,
    BasicOperator,
    Bosonic,
    GenericZhuReducer,
    LocalOperatorBasis,
    NO,
    One,
    OPE,
    ZhuSpace,
    d,
    zhu_circle_product,
    zhu_star_product,
)
from pyope.operator_spaces import _get_conformal_weight


def test_zhu_star_product_matches_weight_one_current_formula():
    """Zhu star product for a weight-1 current uses NO + simple pole.

    Reference: Zhu product definition with a finite sum over OPE modes.
    For wt(J)=1, this is NO(J, b) plus the simple-pole term.
    """

    J = BasicOperator("J_zhu_star", conformal_weight=1)
    K = BasicOperator("K_zhu_star", conformal_weight=1)
    Bosonic(J, K)
    OPE[J, K] = OPE.make([3 * One])

    assert zhu_star_product(J, K) == NO(J, K) + 3 * One


def test_zhu_circle_product_uses_derivative_and_zero_pole_terms():
    """Zhu circle product realizes the $a_{-2}b$ term as NO(d(a), b).

    Reference: circle product is generated from the derivative term and the
    shifted OPE poles, with $a_{-2}b = NO(d(a), b)$.
    """

    J = BasicOperator("J_zhu_circle", conformal_weight=1)
    K = BasicOperator("K_zhu_circle", conformal_weight=1)
    Bosonic(J, K)
    OPE[J, K] = OPE.make([2 * One])

    assert zhu_circle_product(J, K) == NO(d(J), K) + NO(J, K)


def test_generic_zhu_reducer_kills_circle_generators_mod_ov():
    J = BasicOperator("J_zhu_reduce", conformal_weight=1)
    K = BasicOperator("K_zhu_reduce", conformal_weight=1)
    Bosonic(J, K)
    OPE[J, K] = OPE.make([One])

    basis = LocalOperatorBasis([J, K])
    reducer = GenericZhuReducer(basis)
    expr = zhu_circle_product(J, K)

    assert isinstance(reducer, AbstractZhuReducer)
    assert reducer.quotient_normal_form(expr) == 0
    witness = reducer.solve_zhu_witness(expr)
    assert witness.remainder == 0
    assert witness.ov_part == expr


def test_zhu_space_exposes_reducer_style_api_and_multiply():
    J = BasicOperator("J_zhu_space", conformal_weight=1)
    K = BasicOperator("K_zhu_space", conformal_weight=1)
    Bosonic(J, K)
    OPE[J, K] = OPE.make([5 * One])

    basis = LocalOperatorBasis([J, K])
    zhu = ZhuSpace(basis)
    circle_expr = zhu_circle_product(J, K)

    assert zhu.contains(circle_expr, weight=2)
    assert zhu.quotient_normal_form(circle_expr, weight=2) == 0
    assert zhu.is_zero_mod_ov(circle_expr, weight=2)
    witness = zhu.solve_zhu_witness(circle_expr, weight=2)
    assert witness.remainder == 0
    assert zhu.multiply(J, K) == 5
    assert zhu.multiply(J, K, weight=2) == 0


def test_zhu_space_generators_and_basis_are_fixed_weight():
    J = BasicOperator("J_zhu_generators", conformal_weight=1)
    K = BasicOperator("K_zhu_generators", conformal_weight=1)
    Bosonic(J, K)

    basis = LocalOperatorBasis([J, K])
    zhu = ZhuSpace(basis)

    generators = zhu.generators(2)
    assert 2 * NO(J, J) in generators or 2 * NO(J, K) in generators
    independent = zhu.basis(2)
    assert set(independent).issubset(set(generators))
    assert all(
        sp.simplify(_get_conformal_weight(expr) - 2) == 0 for expr in independent
    )
