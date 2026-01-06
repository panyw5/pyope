"""
测试左侧复合算符的 OPE 计算

本测试文件验证 OPE(NO(A,B), C) 的完整实现
"""

import pytest
import sympy as sp
from sympy import Symbol, simplify

from pyope import (
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    d,
    One,
    Zero,
    OPE,
    NO,
    bracket,
    MakeOPE,
    Bosonic,
    Fermionic,
    OPEData,
)


class TestCompositeLeftOPE:
    """测试左侧复合算符的 OPE 计算"""

    def setup_method(self):
        """设置测试环境"""
        # 创建测试用的算符
        self.T = BasisOperator("T", bosonic=True, conformal_weight=2)
        self.J = BasisOperator("J", bosonic=True, conformal_weight=1)
        self.W = BasisOperator("W", bosonic=True, conformal_weight=3)

        # 声明算符类型
        Bosonic(self.T, self.J, self.W)

        # 定义基本的 OPE
        self.c = Symbol('c')
        self.k = Symbol('k')

        # Virasoro OPE: T(z)T(w)
        OPE[self.T, self.T] = MakeOPE([
            self.c/2 * One,  # pole 4
            0,               # pole 3
            2 * self.T,      # pole 2
            d(self.T)        # pole 1
        ])

        # 电流代数 OPE: J(z)J(w)
        OPE[self.J, self.J] = MakeOPE([
            self.k * One,    # pole 2
            0                # pole 1
        ])

        # 混合 OPE: T(z)J(w)
        OPE[self.T, self.J] = MakeOPE([
            self.J,          # pole 2
            d(self.J)        # pole 1
        ])

    def test_composite_left_basic(self):
        """测试基本的左侧复合算符 OPE"""
        # 计算 OPE(NO(J,J), J)
        no_JJ = NO(self.J, self.J)
        ope_result = OPE(no_JJ, self.J)

        # 验证结果不为零
        assert not ope_result.is_zero()

        # 验证最高极点存在
        assert ope_result.max_pole > 0

    def test_composite_left_with_virasoro(self):
        """测试 Virasoro 算符的左侧复合 OPE"""
        # 计算 OPE(NO(T,J), J)
        no_TJ = NO(self.T, self.J)
        ope_result = OPE(no_TJ, self.J)

        # 验证结果
        assert not ope_result.is_zero()
        print(f"\nOPE(NO(T,J), J) 的最高极点: {ope_result.max_pole}")

        # 打印各个极点
        for q in range(ope_result.max_pole, 0, -1):
            pole_q = ope_result.pole(q)
            if pole_q != 0:
                print(f"  pole({q}): {pole_q}")

    def test_composite_left_zero_ope(self):
        """测试当基本 OPE 为零时的情况"""
        # 创建一个未定义 OPE 的算符
        X = BasisOperator("X", bosonic=True)
        Bosonic(X)

        # OPE(NO(X,X), X) 应该返回 NO(NO(X,X), X)
        no_XX = NO(X, X)
        ope_result = OPE(no_XX, X)

        # 验证结果
        assert not ope_result.is_zero()
        pole_0 = ope_result.pole(0)
        assert isinstance(pole_0, NormalOrderedOperator)

    def test_composite_left_linearity(self):
        """测试左侧复合算符 OPE 的线性性"""
        # OPE(NO(J,J), J+J) = OPE(NO(J,J), J) + OPE(NO(J,J), J)
        no_JJ = NO(self.J, self.J)

        ope1 = OPE(no_JJ, self.J + self.J)
        ope2 = OPE(no_JJ, self.J) + OPE(no_JJ, self.J)

        # 验证两者相等（通过比较各个极点）
        # 注意：不使用 simplify，因为它与自定义算符不兼容
        max_pole = max(ope1.max_pole, ope2.max_pole)
        for q in range(1, max_pole + 1):
            pole1 = ope1.pole(q)
            pole2 = ope2.pole(q)
            # 使用 sympy 的 expand 来展开表达式进行比较
            diff = sp.expand(pole1 - pole2)
            assert diff == 0, f"pole({q}): {pole1} != {pole2}"

    def test_composite_left_with_derivative(self):
        """测试包含导数算符的左侧复合 OPE"""
        # 计算 OPE(NO(d(T), J), J)
        dT = d(self.T)
        no_dTJ = NO(dT, self.J)
        ope_result = OPE(no_dTJ, self.J)

        # 验证结果不为零
        assert not ope_result.is_zero()
        print(f"\nOPE(NO(∂T,J), J) 的最高极点: {ope_result.max_pole}")

    def test_composite_left_parity(self):
        """测试费米算符的 parity 处理"""
        # 创建费米算符
        psi = BasisOperator("psi", bosonic=False, conformal_weight=1.5)
        chi = BasisOperator("chi", bosonic=False, conformal_weight=1.5)
        Fermionic(psi, chi)

        # 定义费米算符的 OPE
        OPE[psi, chi] = MakeOPE([
            One,  # pole 1
        ])

        # 计算 OPE(NO(psi, chi), psi)
        no_psi_chi = NO(psi, chi)
        ope_result = OPE(no_psi_chi, psi)

        # 验证结果
        assert not ope_result.is_zero()

    def test_composite_left_comparison_with_right(self):
        """比较左侧和右侧复合算符的 OPE"""
        # 对于对称的情况，左右应该有关系
        # OPE(NO(J,J), J) vs OPE(J, NO(J,J))

        no_JJ = NO(self.J, self.J)
        ope_left = OPE(no_JJ, self.J)
        ope_right = OPE(self.J, no_JJ)

        print(f"\nOPE(NO(J,J), J) 最高极点: {ope_left.max_pole}")
        print(f"OPE(J, NO(J,J)) 最高极点: {ope_right.max_pole}")

        # 两者应该都不为零
        assert not ope_left.is_zero()
        assert not ope_right.is_zero()

    def test_composite_left_nested(self):
        """测试嵌套的正规序算符"""
        # 创建 NO(NO(J,J), J)
        no_JJ = NO(self.J, self.J)
        no_nested = NO(no_JJ, self.J)

        # 计算 OPE(NO(NO(J,J), J), J)
        ope_result = OPE(no_nested, self.J)

        # 验证结果
        assert not ope_result.is_zero()
        print(f"\nOPE(NO(NO(J,J), J), J) 的最高极点: {ope_result.max_pole}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
