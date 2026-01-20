"""
测试 NO(∂J⁰, J⁰) 不应该简化为 Zero 的 bug 修复

这个 bug 是由于 simplify 函数在应用 Thielemans 公式 (2.3.16) 时，
错误地使用 ope_AB.pole(0) 来获取 l=0 项，但 OPE 只包含奇异部分（l>=1），
导致主项 NO(A,B) 被丢失。
"""

import pytest
from pyope import BasisOperator, OPE, simplify, expand_nested_no, NO, d, One
from pyope.ope_data import OPEData
from pyope.registry import ope_registry


@pytest.fixture(autouse=True)
def reset_registry():
    """每个测试前重置 OPE 注册表"""
    ope_registry.clear()
    yield
    ope_registry.clear()


class TestNODerivativeBugFix:
    """测试 NO(∂J⁰, J⁰) bug 修复"""

    def test_no_derivative_not_zero(self):
        """测试 NO(∂J⁰, J⁰) 不应该简化为 Zero"""
        J_zero = BasisOperator("J⁰", conformal_weight=1)
        
        k_val = 1
        OPE[J_zero, J_zero] = OPEData({2: 2 * k_val * One})
        
        # NO(∂J⁰, J⁰) 不应该简化为 Zero
        expr = NO(d(J_zero), J_zero)
        result = simplify(expr)
        
        assert result != 0, "NO(∂J⁰, J⁰) 不应该简化为 Zero"
        # 根据新的排序规则，导数越多越靠左，所以应该保持为 NO(∂J⁰, J⁰)
        assert result == NO(d(J_zero), J_zero)

    def test_kac_moody_consistency(self):
        """测试 Kac-Moody 代数中 simplify 和 expand_nested_no 的一致性"""
        J_plus = BasisOperator("J⁺", conformal_weight=1)
        J_zero = BasisOperator("J⁰", conformal_weight=1)
        J_minus = BasisOperator("J⁻", conformal_weight=1)
        
        k_val = 1
        
        OPE[J_plus, J_zero] = OPEData({1: -2 * J_plus})
        OPE[J_plus, J_minus] = OPEData({2: k_val * One, 1: J_zero})
        OPE[J_zero, J_minus] = OPEData({1: -2 * J_minus})
        OPE[J_zero, J_zero] = OPEData({2: 2 * k_val * One})
        
        # 测试表达式
        expr = NO(J_minus, NO(J_plus, J_zero))
        
        # simplify 和 expand_nested_no 的结果应该在简化后相等
        result_simplify = simplify(expr)
        result_expand = expand_nested_no(expr)
        result_expand_simplified = simplify(result_expand)
        
        assert result_simplify == result_expand_simplified, \
            f"simplify 和 expand_nested_no 的结果应该相等\n" \
            f"simplify: {result_simplify}\n" \
            f"expand_nested_no 简化后: {result_expand_simplified}"

    def test_derivative_ordering(self):
        """测试导数算符的排序规则：导数越多越靠左"""
        T = BasisOperator("T", conformal_weight=2)
        
        # NO(T, ∂T) 应该重排为 NO(∂T, T)
        expr1 = NO(T, d(T))
        result1 = simplify(expr1)
        assert result1 == NO(d(T), T), "NO(T, ∂T) 应该重排为 NO(∂T, T)"
        
        # NO(∂T, T) 应该保持不变
        expr2 = NO(d(T), T)
        result2 = simplify(expr2)
        assert result2 == NO(d(T), T), "NO(∂T, T) 应该保持不变"
        
        # NO(T, ∂²T) 应该重排为 NO(∂²T, T)
        expr3 = NO(T, d(T, 2))
        result3 = simplify(expr3)
        assert result3 == NO(d(T, 2), T), "NO(T, ∂²T) 应该重排为 NO(∂²T, T)"
        
        # NO(∂T, ∂²T) 应该重排为 NO(∂²T, ∂T)
        expr4 = NO(d(T), d(T, 2))
        result4 = simplify(expr4)
        assert result4 == NO(d(T, 2), d(T)), "NO(∂T, ∂²T) 应该重排为 NO(∂²T, ∂T)"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
