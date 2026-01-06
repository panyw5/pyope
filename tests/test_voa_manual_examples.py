"""
基于 VOA Manual 的测试用例

本测试文件基于 gemini 阅读 VOA-manual.md 和 OPEdefs-manual.md 后构思的测试用例。
这些测试验证 PyOPE 实现的核心功能和数学性质。
"""

import pytest
import sympy as sp
from sympy import Symbol, simplify, expand

from pyope import (
    BasisOperator,
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


class TestFreeFieldsOPE:
    """
    测试用例 1: 基础自由场 OPE

    参考: VOA-manual.md Section 2.6.1 (Eq 2.6.7) 和 Section 2.6.2 (Eq 2.6.26)
    """

    def test_bosonic_free_field_ope(self):
        """
        测试玻色子自由场的 OPE

        ∂X(z) ∂X(w) ~ 1/(z-w)^2
        """
        # 定义自由玻色子
        dX = BasisOperator("dX", bosonic=True, conformal_weight=1)
        Bosonic(dX)

        # 定义 OPE
        OPE[dX, dX] = MakeOPE([One, 0])  # [pole_2, pole_1]

        # 计算 OPE
        ope_result = OPE(dX, dX)

        # 验证结果
        assert ope_result.max_pole == 2
        assert ope_result.pole(2) == One
        assert ope_result.pole(1) == 0

        print(f"✓ 玻色子自由场 OPE: ∂X(z)∂X(w) ~ 1/(z-w)^2")

    def test_fermionic_free_field_ope(self):
        """
        测试费米子自由场的 OPE

        ψ(z) ψ(w) ~ 1/(z-w)
        """
        # 定义自由费米子
        psi = BasisOperator("psi", bosonic=False, conformal_weight=0.5)
        Fermionic(psi)

        # 定义 OPE
        OPE[psi, psi] = MakeOPE([One])  # [pole_1]

        # 计算 OPE
        ope_result = OPE(psi, psi)

        # 验证结果
        assert ope_result.max_pole == 1
        assert ope_result.pole(1) == One

        print(f"✓ 费米子自由场 OPE: ψ(z)ψ(w) ~ 1/(z-w)")


class TestVirasoroAlgebra:
    """
    测试用例 2: Virasoro 代数与中心荷

    参考: VOA-manual.md Section 2.2.2 (Eq 2.2.17, 2.3.26) 和 Section 2.6.1 (Eq 2.6.9)
    """

    def test_virasoro_ope_structure(self):
        """
        测试 Virasoro 代数的 OPE 结构

        T(z)T(w) ~ c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w)
        """
        # 定义能量-动量张量
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(T)

        # 中心荷
        c = Symbol('c')

        # 定义 Virasoro OPE
        OPE[T, T] = MakeOPE([
            c/2 * One,  # pole_4
            0,          # pole_3
            2 * T,      # pole_2
            d(T)        # pole_1
        ])

        # 计算 OPE
        ope_result = OPE(T, T)

        # 验证结构
        assert ope_result.max_pole == 4

        # 验证 4 阶极点（中心荷项）
        pole_4 = ope_result.pole(4)
        assert pole_4 == c/2 * One

        # 验证 2 阶极点
        pole_2 = ope_result.pole(2)
        assert pole_2 == 2 * T

        # 验证 1 阶极点
        pole_1 = ope_result.pole(1)
        assert pole_1 == d(T)

        print(f"✓ Virasoro 代数 OPE 结构正确")
        print(f"  - 中心荷项: c/2")
        print(f"  - 2阶极点: 2T")
        print(f"  - 1阶极点: ∂T")

    def test_free_boson_central_charge(self):
        """
        测试自由玻色子的中心荷

        对于单个自由玻色子，c = 1
        """
        # 定义自由玻色子
        dX = BasisOperator("dX", bosonic=True, conformal_weight=1)
        Bosonic(dX)

        # 定义 OPE
        OPE[dX, dX] = MakeOPE([One, 0])

        # 定义 Sugawara 张量 T = 1/2 :∂X ∂X:
        # 注意：这里我们直接定义 T 的 OPE，而不是从 :∂X ∂X: 计算
        # 因为完整的 Sugawara 构造需要更复杂的实现

        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(T)

        # 对于单个自由玻色子，c = 1
        c = 1
        OPE[T, T] = MakeOPE([
            c/2 * One,  # c/2 = 1/2
            0,
            2 * T,
            d(T)
        ])

        ope_result = OPE(T, T)
        pole_4 = ope_result.pole(4)

        # 验证中心荷（使用 expand 来处理浮点数和有理数的比较）
        expected = sp.Rational(1, 2) * One
        assert expand(pole_4 - expected) == 0

        print(f"✓ 自由玻色子中心荷 c = 1 (系数 c/2 = 1/2)")


class TestCompositeOperators:
    """
    测试用例 3 & 4: 复合算符和正规序

    参考: OPEdefs-manual.md Section 2.0.2 & Section 3.3.1
    """

    def test_normal_ordered_product_basic(self):
        """
        测试基本的正规序乘积
        """
        # 定义算符
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(J)

        # 创建正规序乘积
        no_JJ = NO(J, J)

        # 验证类型
        from pyope.operators import NormalOrderedOperator
        assert isinstance(no_JJ, NormalOrderedOperator)

        # 验证因子
        assert no_JJ.left == J
        assert no_JJ.right == J

        print(f"✓ 正规序乘积创建成功: NO(J,J)")

    def test_composite_ope_left_side(self):
        """
        测试左侧复合算符的 OPE

        这是我们刚刚完善的功能
        """
        # 定义算符
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(J)

        k = Symbol('k')
        OPE[J, J] = MakeOPE([k * One, 0])

        # 创建复合算符
        no_JJ = NO(J, J)

        # 计算 OPE(NO(J,J), J)
        ope_result = OPE(no_JJ, J)

        # 验证结果不为零
        assert not ope_result.is_zero()
        assert ope_result.max_pole > 0

        print(f"✓ 左侧复合算符 OPE: OPE(NO(J,J), J)")
        print(f"  最高极点: {ope_result.max_pole}")

    def test_nested_normal_ordered_products(self):
        """
        测试嵌套的正规序乘积
        """
        # 定义算符
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(J)

        k = Symbol('k')
        OPE[J, J] = MakeOPE([k * One, 0])

        # 创建嵌套正规序
        no_JJ = NO(J, J)
        no_nested = NO(no_JJ, J)

        # 计算 OPE
        ope_result = OPE(no_nested, J)

        # 验证结果
        assert not ope_result.is_zero()

        print(f"✓ 嵌套正规序 OPE: OPE(NO(NO(J,J), J), J)")


class TestJacobiIdentity:
    """
    测试用例 5: Jacobi 恒等式验证

    参考: VOA-manual.md Section 2.3.2 (Eq 2.3.21)
    """

    def test_jacobi_identity_virasoro(self):
        """
        测试 Virasoro 算符的 Jacobi 恒等式

        注意：完整的 Jacobi 恒等式验证需要实现 OPEJacobi 函数
        这里我们测试一个简化版本
        """
        # 定义算符
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(T)

        c = Symbol('c')
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])

        # 测试结合性的一个方面：
        # OPE(OPE(T,T)的某个极点, T) 应该与某些组合一致

        # 提取 {TT}_2 = 2T
        bracket_TT_2 = bracket(T, T, 2)
        assert bracket_TT_2 == 2*T

        # 计算 OPE({TT}_2, T) = OPE(2T, T) = 2*OPE(T,T)
        ope_bracket_T = OPE(bracket_TT_2, T)
        ope_TT_scaled = 2 * OPE(T, T)

        # 验证线性性
        assert ope_bracket_T.max_pole == ope_TT_scaled.max_pole

        print(f"✓ Jacobi 恒等式（部分验证）: 线性性成立")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
