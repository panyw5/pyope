"""
测试算符的共形权重 (conformal weight) 属性

测试内容:
1. BasisOperator 的共形权重
2. DerivativeOperator 的共形权重
3. NormalOrderedOperator 的共形权重
4. 复杂算符组合的共形权重
"""

import pytest
from pyope import (
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    d, dn,
    NO,
    Bosonic,
    Fermionic,
)


class TestBasisOperatorConformalWeight:
    """测试基本算符的共形权重"""

    def test_basic_conformal_weight(self):
        """测试基本算符的共形权重设置和获取"""
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        assert T.conformal_weight == 2

        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        assert J.conformal_weight == 1

        psi = BasisOperator("psi", bosonic=False, conformal_weight=1.5)
        assert psi.conformal_weight == 1.5

    def test_fractional_conformal_weight(self):
        """测试分数共形权重（如超对称代数）"""
        G = BasisOperator("G", bosonic=False, conformal_weight=3/2)
        assert G.conformal_weight == 1.5

    def test_no_conformal_weight(self):
        """测试未指定共形权重的算符"""
        A = BasisOperator("A", bosonic=True)
        assert A.conformal_weight is None


class TestDerivativeOperatorConformalWeight:
    """测试导数算符的共形权重"""

    def test_first_derivative(self):
        """
        测试一阶导数算符的共形权重

        规则: weight(∂A) = weight(A) + 1
        """
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        dJ = d(J)

        assert dJ.conformal_weight == 2
        print(f"✓ ∂J 的共形权重: {dJ.conformal_weight} (期望: 2)")

    def test_second_derivative(self):
        """
        测试二阶导数算符的共形权重

        规则: weight(∂²A) = weight(A) + 2
        """
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        d2J = d(d(J))

        assert d2J.conformal_weight == 3
        print(f"✓ ∂²J 的共形权重: {d2J.conformal_weight} (期望: 3)")

    def test_high_order_derivative(self):
        """
        测试高阶导数算符的共形权重

        规则: weight(∂ⁿA) = weight(A) + n
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)

        # 三阶导数 (注意: dn 的参数顺序是 dn(order, operator))
        d3T = dn(3, T)
        assert d3T.conformal_weight == 5

        # 五阶导数
        d5T = dn(5, T)
        assert d5T.conformal_weight == 7

        print(f"✓ ∂³T 的共形权重: {d3T.conformal_weight} (期望: 5)")
        print(f"✓ ∂⁵T 的共形权重: {d5T.conformal_weight} (期望: 7)")

    def test_derivative_no_weight(self):
        """测试对未定义权重的算符求导"""
        A = BasisOperator("A", bosonic=True)
        dA = d(A)

        assert dA.conformal_weight is None


class TestNormalOrderedOperatorConformalWeight:
    """测试正规序算符的共形权重"""

    def test_basic_normal_ordered(self):
        """
        测试基本正规序的共形权重

        规则: weight(NO(A, B)) = weight(A) + weight(B)
        """
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        no_JJ = NO(J, J)

        assert no_JJ.conformal_weight == 2
        print(f"✓ NO(J, J) 的共形权重: {no_JJ.conformal_weight} (期望: 2)")

    def test_sugawara_tensor(self):
        """
        测试 Sugawara 张量的共形权重

        T_sugawara = NO(J, J) 应该有权重 2 (能量-动量张量)
        """
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        T_sugawara = NO(J, J)

        assert T_sugawara.conformal_weight == 2
        print(f"✓ Sugawara 张量的共形权重: {T_sugawara.conformal_weight}")

    def test_mixed_operators(self):
        """测试不同算符的正规序"""
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        J = BasisOperator("J", bosonic=True, conformal_weight=1)

        no_TJ = NO(T, J)
        assert no_TJ.conformal_weight == 3

        print(f"✓ NO(T, J) 的共形权重: {no_TJ.conformal_weight} (期望: 3)")

    def test_fermionic_normal_ordered(self):
        """测试费米算符的正规序"""
        psi = BasisOperator("psi", bosonic=False, conformal_weight=1.5)
        chi = BasisOperator("chi", bosonic=False, conformal_weight=0.5)
        Fermionic(psi, chi)

        no_psi_chi = NO(psi, chi)
        assert no_psi_chi.conformal_weight == 2.0

        print(f"✓ NO(psi, chi) 的共形权重: {no_psi_chi.conformal_weight} (期望: 2.0)")

    def test_normal_ordered_with_derivative(self):
        """
        测试包含导数算符的正规序

        例如: supercurrent G = NO(psi, ∂chi)
        """
        psi = BasisOperator("psi", bosonic=False, conformal_weight=1.5)
        chi = BasisOperator("chi", bosonic=False, conformal_weight=0.5)
        Fermionic(psi, chi)

        # G = NO(psi, ∂chi)
        G = NO(psi, d(chi))

        # weight(G) = weight(psi) + weight(∂chi) = 1.5 + 1.5 = 3.0
        assert G.conformal_weight == 3.0

        print(f"✓ Supercurrent G 的共形权重: {G.conformal_weight} (期望: 3.0)")

    def test_nested_normal_ordered(self):
        """测试嵌套的正规序算符"""
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(J)

        # NO(NO(J, J), J)
        no_JJ = NO(J, J)
        no_nested = NO(no_JJ, J)

        assert no_nested.conformal_weight == 3
        print(f"✓ NO(NO(J,J), J) 的共形权重: {no_nested.conformal_weight} (期望: 3)")

    def test_normal_ordered_no_weight(self):
        """测试未定义权重的算符的正规序"""
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True, conformal_weight=1)

        no_AB = NO(A, B)
        assert no_AB.conformal_weight is None


class TestComplexConformalWeight:
    """测试复杂算符组合的共形权重"""

    def test_w3_algebra(self):
        """
        测试 W₃ 代数中的高共形权重算符

        W 算符通常有权重 3
        """
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        assert W.conformal_weight == 3

        dW = d(W)
        assert dW.conformal_weight == 4

        print(f"✓ W 的共形权重: {W.conformal_weight}")
        print(f"✓ ∂W 的共形权重: {dW.conformal_weight}")

    def test_multiple_derivatives_in_normal_ordered(self):
        """测试正规序中包含多个导数"""
        J = BasisOperator("J", bosonic=True, conformal_weight=1)

        # NO(∂J, ∂J)
        no_dJ_dJ = NO(d(J), d(J))
        assert no_dJ_dJ.conformal_weight == 4

        # NO(J, ∂²J)
        no_J_d2J = NO(J, d(d(J)))
        assert no_J_d2J.conformal_weight == 4

        print(f"✓ NO(∂J, ∂J) 的共形权重: {no_dJ_dJ.conformal_weight}")
        print(f"✓ NO(J, ∂²J) 的共形权重: {no_J_d2J.conformal_weight}")

    def test_deeply_nested_operators(self):
        """测试深度嵌套的算符"""
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(J)

        # NO(NO(NO(J, J), J), J)
        no1 = NO(J, J)           # weight = 2
        no2 = NO(no1, J)         # weight = 3
        no3 = NO(no2, J)         # weight = 4

        assert no3.conformal_weight == 4
        print(f"✓ 深度嵌套算符的共形权重: {no3.conformal_weight}")

    def test_virasoro_primary_field(self):
        """
        测试 Virasoro primary field 的性质

        Primary field φ 应该满足: T(z)φ(w) ~ h*φ/(z-w)² + ...
        其中 h 是 φ 的共形权重
        """
        # 定义一个 primary field
        phi = BasisOperator("phi", bosonic=True, conformal_weight=2.5)

        assert phi.conformal_weight == 2.5

        # ∂φ 的共形权重应该是 h+1
        d_phi = d(phi)
        assert d_phi.conformal_weight == 3.5

        print(f"✓ Primary field φ 的共形权重: {phi.conformal_weight}")
        print(f"✓ ∂φ 的共形权重: {d_phi.conformal_weight}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
