"""
测试 normal_product 函数
"""

import pytest
import sympy as sp

from pyope import (
    BasicOperator,
    NO,
    NO_product,
    NormalOrderedOperator,
    One,
    Zero,
    normal_product,
    d,
)


def test_normal_product_empty():
    """测试空参数列表"""
    result = normal_product()
    assert result == One


def test_normal_product_single():
    """测试单个算符"""
    T = BasicOperator("T", 2)
    result = normal_product(T)
    assert result == T


def test_normal_product_two_operators():
    """测试两个算符"""
    T = BasicOperator("T", 2)
    J = BasicOperator("J", 1)

    result = normal_product(T, J)
    expected = NO(T, J)

    assert result == expected


def test_no_binary_returns_normal_ordered_operator_instance():
    A = BasicOperator("A_no_binary_type", 1)
    B = BasicOperator("B_no_binary_type", 1)

    result = NO(A, B)

    assert isinstance(result, NormalOrderedOperator)


def test_normal_product_three_operators():
    """测试三个算符"""
    T = BasicOperator("T", 2)
    J = BasicOperator("J", 1)
    W = BasicOperator("W", 3)

    result = normal_product(T, J, W)
    expected = NO(T, NO(J, W))

    assert result == expected


def test_normal_product_four_operators():
    """测试四个算符"""
    T = BasicOperator("T", 2)
    J = BasicOperator("J", 1)
    W = BasicOperator("W", 3)
    L = BasicOperator("L", 1)

    result = normal_product(T, J, W, L)
    expected = NO(T, NO(J, NO(W, L)))

    assert result == expected


def test_normal_product_with_derivatives():
    """测试包含导数的算符"""
    T = BasicOperator("T", 2)
    J = BasicOperator("J", 1)

    result = normal_product(d(T), J, T)
    expected = NO(d(T), NO(J, T))

    assert result == expected


def test_normal_product_with_zero():
    """测试包含零算符"""
    T = BasicOperator("T", 2)

    result = normal_product(T, Zero)
    assert result == Zero

    result = normal_product(Zero, T)
    assert result == Zero


def test_normal_product_with_one():
    """测试包含单位算符"""
    T = BasicOperator("T", 2)
    J = BasicOperator("J", 1)

    result = normal_product(T, One, J)
    expected = NO(T, J)
    assert result == expected


def test_normal_product_with_scalars():
    """测试包含标量系数"""
    T = BasicOperator("T", 2)
    J = BasicOperator("J", 1)
    c = sp.Symbol("c")

    result = normal_product(c * T, J)
    expected = c * NO(T, J)

    assert result == expected


def test_normal_product_nested_structure():
    """验证嵌套结构的正确性"""
    A = BasicOperator("A", 1)
    B = BasicOperator("B", 1)
    C = BasicOperator("C", 1)
    D = BasicOperator("D", 1)

    result = normal_product(A, B, C, D)

    # 验证最外层是 NO(A, ...)
    assert isinstance(result, type(NO(A, B)))

    # 手动构建预期结构
    expected = NO(A, NO(B, NO(C, D)))
    assert result == expected


def test_normal_product_comparison_with_manual():
    """对比 normal_product 和手动构建的结果"""
    T = BasicOperator("T", 2)
    J = BasicOperator("J", 1)
    W = BasicOperator("W", 3)

    auto_result = normal_product(T, J, W)
    manual_result = NO(T, NO(J, W))

    assert auto_result == manual_result


def test_normal_product_many_operators():
    """测试多个算符（5个以上）"""
    ops = [BasicOperator(f"O{i}", 1) for i in range(6)]

    result = normal_product(*ops)

    # 手动构建预期结果
    expected = ops[-1]
    for op in reversed(ops[:-1]):
        expected = NO(op, expected)

    assert result == expected


def test_no_accepts_multiple_operators():
    A = BasicOperator("A_no_multi", 1)
    B = BasicOperator("B_no_multi", 1)
    C = BasicOperator("C_no_multi", 1)
    D = BasicOperator("D_no_multi", 1)

    result = NO(A, B, C, D)

    assert result == NO(A, NO(B, NO(C, D)))


def test_no_accepts_operator_list():
    X = BasicOperator("X_no_list", 1)
    Y = BasicOperator("Y_no_list", 1)
    Z = BasicOperator("Z_no_list", 1)
    U = BasicOperator("U_no_list", 1)

    result = NO([X, Y, Z, U])

    assert result == NO(X, NO(Y, NO(Z, U)))


def test_no_product_matches_backward_compatible_alias():
    T = BasicOperator("T_no_product", 2)
    J = BasicOperator("J_no_product", 1)
    W = BasicOperator("W_no_product", 3)

    assert NO_product(T, J, W) == normal_product(T, J, W)
    assert NO_product(T, J, W) == NO(T, J, W)


def test_no_empty_and_singleton_forms():
    T = BasicOperator("T_no_single", 2)

    assert NO() == One
    assert NO(T) == T
    assert NO([]) == One
    assert NO([T]) == T


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
