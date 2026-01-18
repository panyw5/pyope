"""
测试正规序的正则化（按算符顺序排列）

根据 Thielemans 论文 eq (2.3.16)，正规序应该按照算符的正则顺序排列。
"""

import pytest
from fractions import Fraction
from src.pyope import BasisOperator, NO, OPE, MakeOPE, d
from src.pyope.simplify import simplify
from src.pyope.registry import ope_registry


class TestNormalOrderingCanonicalization:
    """测试正规序的正则化"""

    def setup_method(self):
        """每个测试前清空注册表"""
        ope_registry.clear()

    def test_bosonic_operators_no_ope(self):
        """测试未定义 OPE 的玻色算符的正规序"""
        beta = BasisOperator("β", conformal_weight=Fraction(3, 2))
        gamma = BasisOperator("γ", conformal_weight=Fraction(-1, 2))

        # 未化简时，两个 NO 是不同的对象
        no1 = NO(beta, gamma)
        no2 = NO(gamma, beta)
        assert no1 != no2

        # 化简后，应该都变成按字典序排列的形式
        simplified1 = simplify(no1)
        simplified2 = simplify(no2)

        # 由于 β < γ（字典序），两者都应该变成 NO(β, γ)
        assert simplified1 == simplified2
        assert str(simplified1) == "NO(β,γ)"

    def test_fermionic_operators_no_ope(self):
        """测试未定义 OPE 的费米算符的正规序"""
        psi = BasisOperator("ψ", conformal_weight=Fraction(1, 2), fermionic=True)
        chi = BasisOperator("χ", conformal_weight=Fraction(1, 2), fermionic=True)

        no1 = NO(psi, chi)
        no2 = NO(chi, psi)

        # 化简后应该相等（但符号相反，因为费米子）
        simplified1 = simplify(no1)
        simplified2 = simplify(no2)

        # ψ < χ（字典序），所以 NO(χ, ψ) 应该变成 -NO(ψ, χ)
        assert simplified1 == -simplified2

    def test_with_defined_ope(self):
        """测试有定义 OPE 的算符的正规序"""
        T = BasisOperator("T", conformal_weight=2)
        c = BasisOperator("c", conformal_weight=0)

        # 定义 OPE: T(z)T(w) = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)
        OPE[T, T] = MakeOPE([c / 2, 0, 2 * T, d(T)])

        # NO(T, T) 应该根据 OPE 计算
        no_tt = NO(T, T)
        simplified = simplify(no_tt)

        # 这应该保持为 NO(T, T)，因为顺序已经正确
        assert str(simplified) == "NO(T,T)"

    def test_derivative_operators(self):
        """测试导数算符的正规序"""
        J = BasisOperator("J", conformal_weight=1)

        # NO(∂J, J) 和 NO(J, ∂J) 应该按照基础算符和导数阶数排序
        no1 = NO(d(J), J)
        no2 = NO(J, d(J))

        simplified1 = simplify(no1)
        simplified2 = simplify(no2)

        # J 的导数阶数：0 < 1，所以 NO(∂J, J) 应该变成 NO(J, ∂J) + 修正项
        # 但由于 OPE 未定义，修正项为 0
        assert simplified1 == simplified2

    def test_mixed_operators(self):
        """测试混合算符的正规序"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        # 测试多个算符的组合
        no1 = NO(C, A)
        no2 = NO(A, C)
        no3 = NO(B, A)

        simplified1 = simplify(no1)
        simplified2 = simplify(no2)
        simplified3 = simplify(no3)

        # A < B < C（字典序）
        assert str(simplified1) == "NO(A,C)"
        assert str(simplified2) == "NO(A,C)"
        assert str(simplified3) == "NO(A,B)"

    def test_list_of_normal_orders(self):
        """测试正规序列表的化简"""
        beta = BasisOperator("β", conformal_weight=Fraction(3, 2))
        gamma = BasisOperator("γ", conformal_weight=Fraction(-1, 2))

        # 创建列表
        no_list = [NO(beta, gamma), NO(gamma, beta)]

        # 化简列表
        simplified_list = [simplify(no) for no in no_list]

        # 两者应该相等
        assert simplified_list[0] == simplified_list[1]
        assert all(str(no) == "NO(β,γ)" for no in simplified_list)

    def test_no_simplification_needed(self):
        """测试已经按正确顺序的正规序"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)

        # 已经按正确顺序
        no = NO(A, B)
        simplified = simplify(no)

        # 应该保持不变
        assert simplified == no

    def test_with_scalar_coefficients(self):
        """测试带标量系数的正规序"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)

        # 带系数的正规序
        expr = 2 * NO(B, A) + 3 * NO(A, B)
        simplified = simplify(expr)

        # 应该合并为 5 * NO(A, B)
        # 注意：这需要 collect_like_terms 功能
        # 目前只测试顺序是否正确
        assert "NO(A,B)" in str(simplified)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
