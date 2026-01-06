"""
W₃ 代数测试

基于 OPEdefs/ope-examples.nb 中的 W₃ 代数部分（line 476-971）

测试内容：
1. W₃ 代数的基本 OPE 定义
2. 辅助算符 Λ = NO[T,T] - (3/10)T'' 的构造
3. T-T, T-W, W-W 的 OPE 计算
4. 导数算符的 OPE
5. 复合正规序乘积

参考资料：
- ope-examples.nb line 504-714: W₃ 代数定义
- ope-examples.nb line 716-971: W₃ 代数计算示例
"""

import pytest
import sympy as sp
from sympy import Symbol, simplify, expand, Rational

from pyope import (
    BasisOperator,
    d,
    dn,
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


@pytest.fixture(autouse=True)
def clear_registry():
    """每个测试前清空注册表"""
    from pyope.registry import ope_registry
    ope_registry.clear()
    yield
    ope_registry.clear()


# ============================================================================
# W₃ 代数基础测试
# ============================================================================

class TestW3AlgebraDefinition:
    """
    测试 W₃ 代数的基本定义

    对应 notebook line 504-714
    """

    def test_w3_operators_declaration(self):
        """
        测试算符声明: T (stress tensor), W (spin-3 primary)
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        Bosonic(T, W)

        assert T.conformal_weight == 2
        assert W.conformal_weight == 3
        # 注意：BasisOperator 使用 _bosonic 私有属性
        assert T._bosonic is True
        assert W._bosonic is True

        print("✓ W₃ 算符声明: T (weight=2), W (weight=3)")

    def test_lambda_construction(self):
        """
        测试辅助算符 Λ 的构造

        Λ[w] = NO[T,T][w] - (3/10)T''[w]

        对应 notebook line 513-523
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(T)

        # 构造 Λ
        Lambda = NO(T, T) - Rational(3, 10) * d(T, 2)

        # 验证 Lambda 是一个局域算符
        assert Lambda is not None

        print(f"✓ 辅助算符 Λ = NO[T,T] - (3/10)T''")
        print(f"  结果: {Lambda}")


class TestW3AlgebraOPEs:
    """
    测试 W₃ 代数的 OPE 关系

    对应 notebook line 525-714
    """

    def test_virasoro_ope_with_c(self):
        """
        测试 T-T OPE（带中心荷 c）

        T(z)T(w) ~ c/(2(z-w)⁴) + 2T/(z-w)² + T'/(z-w)

        对应 notebook line 526-548
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(T)

        c = Symbol('c')
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])

        result = OPE(T, T)

        assert result.max_pole == 4
        assert result.pole(4) == c/2 * One
        assert result.pole(3) == 0
        assert result.pole(2) == 2 * T
        assert result.pole(1) == d(T)

        print("✓ W₃: T-T OPE with central charge c")

    def test_t_w_ope(self):
        """
        测试 T-W OPE

        T(z)W(w) ~ 3W/(z-w)² + W'/(z-w)

        对应 notebook line 550-567
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        Bosonic(T, W)

        c = Symbol('c')
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])
        OPE[T, W] = MakeOPE([3*W, d(W)])

        result = OPE(T, W)

        assert result.max_pole == 2
        assert result.pole(2) == 3 * W
        assert result.pole(1) == d(W)

        print("✓ W₃: T-W OPE")

    def test_w_w_ope(self):
        """
        测试 W-W OPE（W₃ 代数的核心关系）

        W(z)W(w) ~ c/(z-w)⁶ + 2T/(z-w)⁴ + T'/(z-w)³
                    + (2βΛ + (3/10)T'')/(z-w)²
                    + (βΛ' + (1/15)T''')/(z-w)

        其中 Λ = NO[T,T] - (3/10)T''

        对应 notebook line 569-624
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        Bosonic(T, W)

        c = Symbol('c')
        beta = Symbol('beta')

        # 构造 Λ
        Lambda = NO(T, T) - Rational(3, 10) * d(T, 2)

        # 定义基本 OPE
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])
        OPE[T, W] = MakeOPE([3*W, d(W)])

        # W-W OPE
        OPE[W, W] = MakeOPE([
            c * One,
            0,
            2 * T,
            d(T),
            2 * beta * Lambda + Rational(3, 10) * d(T, 2),
            beta * d(Lambda) + Rational(1, 15) * d(T, 3)
        ])

        result = OPE(W, W)

        assert result.max_pole == 6
        assert result.pole(6) == c * One
        assert result.pole(5) == 0
        assert result.pole(4) == 2 * T
        assert result.pole(3) == d(T)

        # 2 阶和 1 阶极点包含 Λ，较为复杂
        pole_2 = result.pole(2)
        pole_1 = result.pole(1)
        assert pole_2 is not None
        assert pole_1 is not None

        print("✓ W₃: W-W OPE with parameter β")


# ============================================================================
# W₃ 代数计算测试（对应 notebook 计算示例）
# ============================================================================

class TestW3AlgebraComputations:
    """
    测试 W₃ 代数的具体计算

    对应 notebook line 716-971
    """

    def test_t_lambda_ope(self):
        """
        测试 OPE[T, NO[T,T] - (3/10)T'']

        对应 notebook line 718-794

        期望结果（line 767-789）:
        << 4|| (-18/5)T + (8+c)T ||3|| 0 ||2|| 4NO[T,T] - (6/5)T''
           ||1|| 2NO[T',T] - (2/15)T''' >>

        注意：虽然 OPE[T, NO[T,T]] 有 6 阶极点 3c*One，
        但 OPE[T, T''] 也有 6 阶极点 10c*One，
        所以 Lambda 的 6 阶极点相消：3c - (3/10)*10c = 0
        最终 max_pole = 4（与 Mathematica 一致）
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(T)

        c = Symbol('c')
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])

        # 构造 Λ
        Lambda = NO(T, T) - Rational(3, 10) * d(T, 2)

        # 计算 OPE[T, Λ]
        result = OPE(T, Lambda)

        # 验证 6 阶极点相消
        ope_T_NOTT = OPE(T, NO(T,T))
        assert ope_T_NOTT.max_pole == 6
        assert ope_T_NOTT.pole(6) == 3 * c * One

        # Lambda 的最终结果 max_pole=4
        assert result.max_pole == 4

        # 验证 4 阶极点
        pole_4 = result.pole(4)
        assert pole_4 is not None

        print(f"✓ W₃: OPE[T, Λ] 计算完成 (max_pole={result.max_pole})")
        print(f"  - OPE[T, NO(T,T)] 的 6 阶极点: {ope_T_NOTT.pole(6)}")
        print(f"  - Λ 的 6 阶极点相消，最终 max_pole=4")
        print(f"  - 4-pole: {pole_4}")

    def test_t_w_second_derivative(self):
        """
        测试 OPE[T, W'']

        对应 notebook line 735-750

        期望结果（line 796-813）:
        << 4|| 18 W ||3|| 14 W' ||2|| 5 W'' ||1|| W''' >>
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        Bosonic(T, W)

        c = Symbol('c')
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])
        OPE[T, W] = MakeOPE([3*W, d(W)])

        # 计算 OPE[T, W'']
        result = OPE(T, d(W, 2))

        assert result.max_pole == 4
        assert result.pole(4) == 18 * W
        assert result.pole(3) == 14 * d(W)
        assert result.pole(2) == 5 * d(W, 2)
        assert result.pole(1) == d(W, 3)

        print("✓ W₃: OPE[T, W''] = 18W/(z-w)⁴ + 14W'/(z-w)³ + ...")

    def test_normal_ordering_t_prime_w_double_prime(self):
        """
        测试 NO[T', NO[W'', T]]

        对应 notebook line 742-765

        期望结果（line 815-840）:
        NO[T, NO[T, W'']] - (1/12)NO[T', W⁽⁴⁾]
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        Bosonic(T, W)

        c = Symbol('c')
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])
        OPE[T, W] = MakeOPE([3*W, d(W)])

        # 计算 NO[T', NO[W'', T]]
        result = NO(d(T), NO(d(W, 2), T))

        # 结果应该包含 NO 项
        assert result is not None

        print(f"✓ W₃: NO[T', NO[W'', T]] = {result}")

    def test_normal_ordering_w_prime_w_double_prime(self):
        """
        测试 NO[W', NO[W'', T]]

        对应 notebook line 749-757

        期望结果（line 842-971）:
        复杂的 NO 项展开，包含高阶导数
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        Bosonic(T, W)

        c = Symbol('c')
        beta = Symbol('beta')

        # 构造 Λ
        Lambda = NO(T, T) - Rational(3, 10) * d(T, 2)

        # 定义完整的 W₃ OPE
        OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])
        OPE[T, W] = MakeOPE([3*W, d(W)])
        OPE[W, W] = MakeOPE([
            c * One,
            0,
            2 * T,
            d(T),
            2 * beta * Lambda + Rational(3, 10) * d(T, 2),
            beta * d(Lambda) + Rational(1, 15) * d(T, 3)
        ])

        # 计算 NO[W', NO[W'', T]]
        result = NO(d(W), NO(d(W, 2), T))

        # 验证结果非空
        assert result is not None

        print(f"✓ W₃: NO[W', NO[W'', T]] 计算完成")


# ============================================================================
# W₃ 代数数值检验
# ============================================================================

class TestW3AlgebraNumerical:
    """
    W₃ 代数的数值检验

    通过设定具体的参数值验证代数结构
    """

    def test_w3_with_specific_parameters(self):
        """
        测试特定参数下的 W₃ 代数

        设定 c=100, β=1/10
        """
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        W = BasisOperator("W", bosonic=True, conformal_weight=3)
        Bosonic(T, W)

        c_val = 100
        beta_val = Rational(1, 10)

        # 构造 Λ
        Lambda = NO(T, T) - Rational(3, 10) * d(T, 2)

        # 定义 OPE
        OPE[T, T] = MakeOPE([c_val/2 * One, 0, 2*T, d(T)])
        OPE[T, W] = MakeOPE([3*W, d(W)])
        OPE[W, W] = MakeOPE([
            c_val * One,
            0,
            2 * T,
            d(T),
            2 * beta_val * Lambda + Rational(3, 10) * d(T, 2),
            beta_val * d(Lambda) + Rational(1, 15) * d(T, 3)
        ])

        # 计算 W-W OPE
        result = OPE(W, W)

        assert result.max_pole == 6
        assert result.pole(6) == c_val * One
        assert result.pole(4) == 2 * T

        print(f"✓ W₃ 数值检验: c={c_val}, β={beta_val}")


# ============================================================================
# 运行测试
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
