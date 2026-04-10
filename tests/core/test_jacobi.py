"""测试 Jacobi 恒等式结果归一。"""

import sympy as sp

from pyope import (
    BasicOperator,
    Bosonic,
    MakeOPE,
    OPE,
    One,
    Zero,
    d,
    check_jacobi_identity,
    verify_jacobi_identity,
)
from pyope.jacobi import _normalize_jacobi_value
from pyope.ope_data import OPEData


def test_normalize_jacobi_value_zero_variants():
    T = BasicOperator("T")
    Bosonic(T)

    assert _normalize_jacobi_value(0) == Zero
    assert _normalize_jacobi_value(Zero) == Zero
    assert _normalize_jacobi_value(3 * Zero) == Zero
    assert _normalize_jacobi_value(OPEData({})) == Zero
    assert _normalize_jacobi_value(2 * T) == 2 * T


def test_check_jacobi_identity_returns_normalized_zeros():
    from pyope.registry import ope_registry

    ope_registry.clear()

    T = BasicOperator("T", conformal_weight=2)
    Bosonic(T)

    c = sp.Symbol("c")
    OPE[T, T] = MakeOPE([c / 2 * One, 0, 2 * T, d(T)])

    result = check_jacobi_identity(T, T, T)

    assert result
    assert all(value == Zero for row in result for value in row)
    assert all(value is Zero for row in result for value in row)


def test_verify_jacobi_identity_uses_same_normalization():
    from pyope.registry import ope_registry

    ope_registry.clear()

    T = BasicOperator("T", conformal_weight=2)
    Bosonic(T)

    c = sp.Symbol("c")
    OPE[T, T] = MakeOPE([c / 2 * One, 0, 2 * T, d(T)])

    assert verify_jacobi_identity(T, T, T) is True
