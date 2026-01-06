"""
OPEdefs.m 与 PyOPE 功能对比测试

本测试文件用于检测 OPEdefs.m (Mathematica) 与 PyOPE (Python) 之间的功能差异。
"""

import pytest
import sympy as sp
from sympy import Symbol

from pyope import (
    # 算符类
    Operator,
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    # 导数函数
    d,
    dn,
    # 局域算符
    LocalOperator,
    OperatorSum,
    OperatorProduct,
    is_local_operator,
    extract_scalar_operator,
    get_operator_parity,
    simplify_operator_expr,
    collect_operator_terms,
    # 常数算符
    ConstantOperator,
    One,
    Zero,
    Delta,
    # OPE 数据
    OPEData,
)


# ============================================================================
# 测试 1: 算符基本属性测试
# ============================================================================

def test_operator_basic_properties():
    """测试算符的基本属性"""
    # 创建玻色算符
    T = BasisOperator("T", bosonic=True)
    assert T.is_bosonic == True
    assert T.is_fermionic == False
    assert T.parity == 0
    
    # 创建费米算符
    psi = BasisOperator("psi", bosonic=False)
    assert psi.is_bosonic == False
    assert psi.is_fermionic == True
    assert psi.parity == 1
    
    # 测试索引算符
    i = Symbol('i')
    J = BasisOperator("J", bosonic=True, indexed=True)
    J_i = J[i]
    assert J_i.base_name == "J"
    assert J_i.indices == (i,)
    
    # 测试唯一性和哈希
    T2 = BasisOperator("T", bosonic=True)
    assert T == T2
    assert hash(T) == hash(T2)
    
    # 测试字符串表示
    assert "BasisOperator" in repr(T)


# ============================================================================
# 测试 2: 导数运算测试
# ============================================================================

def test_derivative_operations():
    """测试导数运算的线性性质"""
    T = BasisOperator("T", bosonic=True)
    
    # 测试一阶导数
    dT = d(T)
    assert isinstance(dT, DerivativeOperator)
    assert dT.base == T
    assert dT.order == 1
    
    # 测试二阶导数
    d2T = d(T, 2)
    assert d2T.order == 2
    
    # 测试线性性质：d(c*A) = c*d(A)
    assert d(2 * T) == 2 * d(T)
    
    # 测试线性性质：d(A+B) = d(A) + d(B)
    W = BasisOperator("W", bosonic=True)
    assert d(T + W) == d(T) + d(W)
    
    # 测试嵌套导数
    assert d(d(T)) == d(T, 2)
    
    # 测试 dn 函数
    d3T = dn(3, T)
    assert d3T.order == 3
    assert dn(2, d(T)) == d(T, 3)


# ============================================================================
# 测试 3: OPEData 基本操作测试
# ============================================================================

def test_ope_data_operations():
    """测试 OPEData 的基本操作"""
    T = BasisOperator("T", bosonic=True)
    
    # 从字典创建 OPEData
    poles = {4: sp.Rational(1, 2), 2: 2*T, 1: d(T)}
    ope = OPEData(poles)
    
    # 测试极点提取
    assert ope.pole(4) == sp.Rational(1, 2)
    assert ope.pole(2) == 2*T
    assert ope.pole(1) == d(T)
    assert ope.pole(3) == 0
    
    # 测试最高阶极点
    assert ope.max_pole == 4
    
    # 测试加减法
    ope2 = OPEData({2: T, 1: d(T)})
    ope_sum = ope + ope2
    assert ope_sum.pole(2) == 3*T
    
    # 测试数乘
    ope_scaled = 2 * ope
    assert ope_scaled.pole(2) == 4*T
    
    # 测试减法
    ope_diff = ope - ope2
    assert ope_diff.pole(2) == T
    
    # 测试从列表创建（Mathematica 风格）
    pole_list = [sp.Rational(1, 2), 0, 2*T, d(T)]
    ope_from_list = OPEData.from_list(pole_list)
    assert ope_from_list.pole(4) == sp.Rational(1, 2)
    assert ope_from_list.pole(3) == 0
    assert ope_from_list.pole(2) == 2*T
    assert ope_from_list.pole(1) == d(T)
    
    # 测试空 OPEData
    empty_ope = OPEData()
    assert empty_ope.is_zero()
    assert empty_ope.max_pole == 0


# ============================================================================
# 测试 4: 局域算符运算测试
# ============================================================================

def test_local_operator_operations():
    """测试局域算符的运算"""
    T = BasisOperator("T", bosonic=True)
    
    # 测试加法
    sum_op = T + T
    assert is_local_operator(sum_op)
    
    # 测试数乘
    scaled_op = 3 * T
    assert is_local_operator(scaled_op)
    
    # 测试提取标量和算符
    coeff, op = extract_scalar_operator(2 * T)
    assert coeff == 2
    assert op == T
    
    # 测试获取 parity
    assert get_operator_parity(T) == 0
    
    psi = BasisOperator("psi", bosonic=False)
    assert get_operator_parity(psi) == 1
    
    # 测试简化表达式
    expr = T + T
    simplified = simplify_operator_expr(expr)
    # sympy 会自动简化
    assert simplified == 2*T
    
    # 测试收集同类项
    expr = 2*T + 3*T
    terms = collect_operator_terms(expr)
    assert T in terms
    assert terms[T] == 5


# ============================================================================
# 测试 5: 常数算符测试
# ============================================================================

def test_constant_operators():
    """测试常数算符"""
    # 测试 One
    assert One.is_bosonic == True
    assert One.parity == 0
    
    # 测试 Zero
    assert Zero.is_bosonic == True
    assert Zero.parity == 0
    
    # 测试 Delta 函数
    assert Delta(1, 1) == 1
    assert Delta(1, 2) == 0
    i = Symbol('i')
    assert Delta(i, i) == 1
    j = Symbol('j')
    delta_ij = Delta(i, j)
    # 预期返回 KroneckerDelta 符号形式
    assert delta_ij == sp.KroneckerDelta(i, j)


# ============================================================================
# 测试 6: 未实现功能测试（预期失败）
# ============================================================================

def test_unimplemented_ope_calculation():
    """测试未实现的 OPE 计算功能（预期失败）"""
    T = BasisOperator("T", bosonic=True)
    
    # 这些功能尚未实现，应该抛出异常或返回错误
    try:
        # OPE 计算未实现
        from pyope import OPE
        result = OPE(T, T)
        pytest.fail("OPE function should not be implemented yet")
    except (ImportError, NameError, NotImplementedError):
        pass
    
    try:
        # bracket 计算未实现
        from pyope import bracket
        result = bracket(T, T, 1)
        pytest.fail("bracket function should not be implemented yet")
    except (ImportError, NameError, NotImplementedError):
        pass


# ============================================================================
# 测试 7: 导数莱布尼茨律测试（未完全实现）
# ============================================================================

def test_derivative_leibniz_rule():
    """测试导数的莱布尼茨律（未完全实现）"""
    T = BasisOperator("T", bosonic=True)
    
    # 当前实现只支持线性性质，不支持莱布尼茨律
    # 例如：d(NO(A,B)) 应该等于 NO(d(A),B) + NO(A,d(B))
    # 但这个功能尚未实现
    
    # 测试当前支持的线性性质
    assert d(T + T) == d(T) + d(T)
    assert d(2 * T) == 2 * d(T)
    
    # 测试导数的 parity
    dT = d(T)
    assert dT.parity == T.parity


# ============================================================================
# 测试 8: NormalOrderedOperator 基本测试
# ============================================================================

def test_normal_ordered_operator():
    """测试 NormalOrderedOperator 的基本功能"""
    T = BasisOperator("T", bosonic=True)
    W = BasisOperator("W", bosonic=True)
    
    # 创建正规序算符
    no_TW = NormalOrderedOperator(T, W)
    assert no_TW.left == T
    assert no_TW.right == W
    assert no_TW.factors == (T, W)
    
    # 测试 parity
    assert no_TW.parity == (T.parity + W.parity) % 2
    
    # 测试费米算符的 parity
    psi = BasisOperator("psi", bosonic=False)
    no_psi_psi = NormalOrderedOperator(psi, psi)
    assert no_psi_psi.parity == 0  # 1 + 1 = 2, 2 % 2 = 0
    
    # 测试相等性和哈希
    no_TW2 = NormalOrderedOperator(T, W)
    assert no_TW == no_TW2
    assert hash(no_TW) == hash(no_TW2)


# ============================================================================
# 测试 9: OPEData 简化测试
# ============================================================================

def test_ope_data_simplify():
    """测试 OPEData 的简化功能"""
    T = BasisOperator("T", bosonic=True)
    
    # 创建包含简化形式的 OPEData
    poles = {
        2: 2*T + 3*T,
        1: T - T
    }
    ope = OPEData(poles)
    
    # 简化
    simplified = ope.simplify()
    # sympy.simplify 应该自动简化表达式
    assert simplified.pole(2) == 5*T
    # T - T 应该被简化为 0
    assert simplified.pole(1) == 0 or simplified.pole(1) == sp.Integer(0)


# ============================================================================
# 测试 10: 复杂表达式测试
# ============================================================================

def test_complex_expressions():
    """测试复杂的算符表达式"""
    T = BasisOperator("T", bosonic=True)
    W = BasisOperator("W", bosonic=True)
    
    # 创建复杂表达式
    expr = 2*T + 3*W - T + 2*W
    assert is_local_operator(expr)
    
    # 提取标量和算符
    coeff, op = extract_scalar_operator(5*T)
    assert coeff == 5
    assert op == T
    
    # 测试导数的导数
    d2T = d(d(T))
    assert d2T.order == 2
    assert d2T.base == T


# ============================================================================
# 测试 11: 边界情况测试
# ============================================================================

def test_edge_cases():
    """测试边界情况"""
    T = BasisOperator("T", bosonic=True)
    
    # 测试零阶导数
    d0T = d(T, 0)
    # 当前实现可能不支持 0 阶导数，或者直接返回 T 本身
    
    # 测试负算符
    neg_op = -T
    assert is_local_operator(neg_op)
    
    # 测试零算符
    zero_op = 0 * T
    # sympy 会自动将其简化为 0
    
    # 测试空 OPEData
    empty_ope = OPEData()
    assert empty_ope.is_zero()
    assert empty_ope.pole(1) == 0


# ============================================================================
# 测试 12: 索引算符测试
# ============================================================================

def test_indexed_operators():
    """测试索引算符"""
    i = Symbol('i')
    j = Symbol('j')
    
    J = BasisOperator("J", bosonic=True, indexed=True)
    
    # 创建索引算符
    J_i = J[i]
    J_j = J[j]
    
    # 测试属性
    assert J_i.base_name == "J"
    assert J_i.indices == (i,)
    assert J_i.is_bosonic == True
    
    # 测试相等性
    J_i2 = J[i]
    assert J_i == J_i2
    
    # 测试不同索引的不相等性
    assert J_i != J_j
    
    # 测试非索引算符的索引操作
    T = BasisOperator("T", bosonic=True)
    try:
        T_i = T[i]
        pytest.fail("Non-indexed operator should not support indexing")
    except TypeError:
        pass


# ============================================================================
# 测试 13: 导数算符的 parity 测试
# ============================================================================

def test_derivative_operator_parity():
    """测试导数算符的 parity"""
    T = BasisOperator("T", bosonic=True)
    psi = BasisOperator("psi", bosonic=False)
    
    # 玻色算符的导数应该是玻色子
    dT = d(T)
    assert dT.parity == 0
    
    # 费米算符的导数应该是费米子
    dpsi = d(psi)
    assert dpsi.parity == 1
    
    # 高阶导数
    d2T = d(T, 2)
    assert d2T.parity == 0


# ============================================================================
# 测试 14: OPEData 相等性测试
# ============================================================================

def test_ope_data_equality():
    """测试 OPEData 的相等性"""
    T = BasisOperator("T", bosonic=True)
    
    poles1 = {2: 2*T, 1: d(T)}
    poles2 = {2: 2*T, 1: d(T)}
    
    ope1 = OPEData(poles1)
    ope2 = OPEData(poles2)
    
    assert ope1 == ope2
    
    # 测试不相等的情况
    poles3 = {2: 3*T, 1: d(T)}
    ope3 = OPEData(poles3)
    assert ope1 != ope3


# ============================================================================
# 测试 15: OPEData 字符串表示测试
# ============================================================================

def test_ope_data_string_representation():
    """测试 OPEData 的字符串表示"""
    T = BasisOperator("T", bosonic=True)
    
    # 空 OPEData
    empty_ope = OPEData()
    assert "OPEData" in str(empty_ope)
    
    # 非空 OPEData
    poles = {2: 2*T, 1: d(T)}
    ope = OPEData(poles)
    assert "OPEData" in str(ope)
    assert "2:" in str(ope) or "2" in str(ope)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
