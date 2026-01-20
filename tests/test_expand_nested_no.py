"""
测试 expand_nested_no 函数

验证 OPEdefs.m 风格的嵌套正规序展开
"""

import pytest
from pyope import (
    BasisOperator,
    OPE,
    NO,
    simplify,
    expand_nested_no,
    One,
    Zero,
    d,
)
from pyope.api import MakeOPE


def test_kac_moody_nested_no_expansion():
    """
    测试 Kac-Moody 代数中的嵌套正规序展开

    验证 expand_nested_no(NO(J⁻, NO(J⁺, J⁰))) 给出 Mathematica 风格的输出：
    NO(J⁺, NO(J⁰, J⁻)) + 2*NO(J⁺, ∂J⁻) - NO(∂J⁰, J⁰)

    而不是 simplify 的完全规约结果：2*NO(J⁺, ∂J⁻)
    """
    J_plus = BasisOperator("J⁺", conformal_weight=1)
    J_zero = BasisOperator("J⁰", conformal_weight=1)
    J_minus = BasisOperator("J⁻", conformal_weight=1)

    k_val = 1

    OPE[J_plus, J_zero] = MakeOPE([-2 * J_plus])
    OPE[J_plus, J_minus] = MakeOPE([k_val * One, J_zero])
    OPE[J_zero, J_minus] = MakeOPE([-2 * J_minus])
    OPE[J_zero, J_zero] = MakeOPE([2 * k_val * One, 0])

    expr = NO(J_minus, NO(J_plus, J_zero))

    # 使用 expand_nested_no（OPEdefs 风格）
    expanded = expand_nested_no(expr)

    # 期望的结果（Mathematica 风格）
    # NO(J⁺, NO(J⁰, J⁻)) + 2*NO(J⁺, ∂J⁻) - NO(∂J⁰, J⁰)
    expected = (
        NO(J_plus, NO(J_zero, J_minus))
        + 2 * NO(J_plus, d(J_minus))
        - NO(d(J_zero), J_zero)
    )

    # 验证展开结果包含三重嵌套项
    assert str(expanded) != str(simplify(expr)), (
        "expand_nested_no 应该与 simplify 给出不同的结果"
    )

    # 验证展开结果与期望一致（通过简化差值）
    diff = simplify(expanded - expected)
    assert diff == Zero or diff == 0, f"展开结果不匹配，差值为: {diff}"


def test_expand_vs_simplify():
    """
    测试 expand_nested_no 与 simplify 的区别
    """
    J_plus = BasisOperator("J⁺", conformal_weight=1)
    J_zero = BasisOperator("J⁰", conformal_weight=1)
    J_minus = BasisOperator("J⁻", conformal_weight=1)

    k_val = 1

    OPE[J_plus, J_zero] = MakeOPE([-2 * J_plus])
    OPE[J_plus, J_minus] = MakeOPE([k_val * One, J_zero])
    OPE[J_zero, J_minus] = MakeOPE([-2 * J_minus])
    OPE[J_zero, J_zero] = MakeOPE([2 * k_val * One, 0])

    expr = NO(J_minus, NO(J_plus, J_zero))

    # expand_nested_no 保留三重嵌套项
    expanded = expand_nested_no(expr)

    # simplify 完全规约
    simplified = simplify(expr)

    # 验证 simplify 的结果是 2*NO(J⁺, ∂J⁻)
    expected_simplified = 2 * NO(J_plus, d(J_minus))
    diff_simplified = simplify(simplified - expected_simplified)
    assert diff_simplified == Zero or diff_simplified == 0, (
        f"simplify 结果不匹配，差值为: {diff_simplified}"
    )

    # 验证两者不同
    assert str(expanded) != str(simplified), (
        "expand_nested_no 和 simplify 应该给出不同的结果"
    )

    # 注意：expand_nested_no 保留三重嵌套项，而 simplify 会进一步规约
    # 所以它们在形式上不同，但在物理上等价（通过进一步简化可以证明）
    # 这里我们只验证它们确实不同
    print(f"expanded: {expanded}")
    print(f"simplified: {simplified}")


def test_expand_nested_no_depth():
    """
    测试 expand_nested_no 的深度控制
    """
    A = BasisOperator("A", conformal_weight=1)
    B = BasisOperator("B", conformal_weight=1)
    C = BasisOperator("C", conformal_weight=1)

    # 设置简单的 OPE（无收缩）
    OPE[A, B] = MakeOPE([])
    OPE[A, C] = MakeOPE([])
    OPE[B, C] = MakeOPE([])

    # 三重嵌套
    expr = NO(NO(A, B), C)

    # 深度 0：不展开
    result_0 = expand_nested_no(expr, max_depth=0)
    assert result_0 == expr, "深度 0 应该不展开"

    # 深度 1：展开一层
    result_1 = expand_nested_no(expr, max_depth=1)
    # 应该展开为 NO(A, NO(B, C))（因为没有 OPE 收缩）
    expected_1 = NO(A, NO(B, C))
    assert str(result_1) == str(expected_1), (
        f"深度 1 展开结果不匹配: {result_1} vs {expected_1}"
    )


def test_expand_nested_no_with_derivatives():
    """
    测试 expand_nested_no 对导数的处理
    """
    T = BasisOperator("T", conformal_weight=2)

    OPE[T, T] = MakeOPE([One / 2, 0, 2 * T, 0])

    # 包含导数的嵌套 NO
    expr = NO(T, NO(T, d(T)))

    # 默认不展开导数
    expanded = expand_nested_no(expr, expand_derivatives=False)

    # 验证结果中仍包含 d(T)
    assert "∂T" in str(expanded) or "d(T)" in str(expanded), (
        "expand_derivatives=False 应该保留导数形式"
    )


def test_expand_nested_no_linearity():
    """
    测试 expand_nested_no 的线性性
    """
    A = BasisOperator("A", conformal_weight=1)
    B = BasisOperator("B", conformal_weight=1)
    C = BasisOperator("C", conformal_weight=1)

    OPE[A, B] = MakeOPE([])
    OPE[A, C] = MakeOPE([])
    OPE[B, C] = MakeOPE([])

    expr1 = NO(NO(A, B), C)
    expr2 = NO(NO(B, A), C)

    # 展开和式
    result = expand_nested_no(expr1 + expr2)

    # 应该等于分别展开后相加
    expected = expand_nested_no(expr1) + expand_nested_no(expr2)

    diff = simplify(result - expected)
    assert diff == Zero or diff == 0, f"expand_nested_no 应该满足线性性，差值为: {diff}"


def test_expand_nested_no_scalar_multiplication():
    """
    测试 expand_nested_no 对标量乘法的处理
    """
    A = BasisOperator("A", conformal_weight=1)
    B = BasisOperator("B", conformal_weight=1)
    C = BasisOperator("C", conformal_weight=1)

    OPE[A, B] = MakeOPE([])
    OPE[A, C] = MakeOPE([])
    OPE[B, C] = MakeOPE([])

    expr = NO(NO(A, B), C)

    # 标量乘法
    result = expand_nested_no(3 * expr)
    expected = 3 * expand_nested_no(expr)

    diff = simplify(result - expected)
    assert diff == Zero or diff == 0, (
        f"expand_nested_no 应该保持标量乘法，差值为: {diff}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
