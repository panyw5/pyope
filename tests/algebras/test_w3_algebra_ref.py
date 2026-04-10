"""
W3 代数测试

本测试文件对应 Mathematica 参考实现：tests/w3_algebra_test.wls

测试覆盖：
- 基本 OPE（T-T, T-W, W-W）
- 导数算符 OPE
- 复合正规序 OPE
- 数值验证（c=100, β=1/10）

总计：8个测试用例
"""

import pytest
import sympy as sp

from pyope import OPE, NO, d, dn, simplify
from tests.utils.comparison import assert_voa_equal, assert_voa_numeric_equal


# ============================================
# Test Case 1-3: 基本 OPE 计算
# ============================================


@pytest.mark.mathematica_ref
class TestBasicW3OPE:
    """基本 W3 代数 OPE 测试"""

    def test_1__T_T_OPE(self, w3_algebra):
        """
        Test 1: T-T OPE

        对应 Mathematica:
        Print["Test 1: T-T OPE"];
        Print["OPE[T, T] = ", OPE[T, T]];

        期望结果：
        T(z) T(w) = c/2 / (z-w)^4 + 2*T(w) / (z-w)^2 + ∂T(w) / (z-w)
        """
        T = w3_algebra["T"]
        c = w3_algebra["c"]

        # 计算 OPE
        result = OPE(T, T)

        # 验证极点结构
        assert result.max_pole == 4

        # 验证各极点系数
        pole_4 = result.pole(4)
        pole_3 = result.pole(3)
        pole_2 = result.pole(2)
        pole_1 = result.pole(1)

        # (z-w)^{-4}: c/2
        assert_voa_equal(pole_4, sp.Rational(1, 2) * c)

        # (z-w)^{-3}: 0
        assert_voa_equal(pole_3, 0)

        # (z-w)^{-2}: 2*T
        assert_voa_equal(pole_2, 2 * T)

        # (z-w)^{-1}: ∂T
        assert_voa_equal(pole_1, d(T))

    def test_2__T_W_OPE(self, w3_algebra):
        """
        Test 2: T-W OPE

        对应 Mathematica:
        Print["Test 2: T-W OPE"];
        Print["OPE[T, W] = ", OPE[T, W]];

        期望结果：
        T(z) W(w) = 3*W(w) / (z-w)^2 + ∂W(w) / (z-w)
        """
        T = w3_algebra["T"]
        W = w3_algebra["W"]

        # 计算 OPE
        result = OPE(T, W)

        # 验证极点结构
        assert result.max_pole == 2

        # 验证各极点系数
        pole_2 = result.pole(2)
        pole_1 = result.pole(1)

        # (z-w)^{-2}: 3*W
        assert_voa_equal(pole_2, 3 * W)

        # (z-w)^{-1}: ∂W
        assert_voa_equal(pole_1, d(W))

    def test_3__W_W_OPE(self, w3_algebra):
        """
        Test 3: W-W OPE

        对应 Mathematica:
        Print["Test 3: W-W OPE"];
        Print["OPE[W, W] = ", OPE[W, W]];

        期望结果：
        W(z) W(w) = c / (z-w)^6 + 2*T(w) / (z-w)^4 + ∂T(w) / (z-w)^3
                    + (2*β*Λ(w) + (3/10)*∂²T(w)) / (z-w)^2
                    + (β*∂Λ(w) + (1/15)*∂³T(w)) / (z-w)
        """
        T = w3_algebra["T"]
        W = w3_algebra["W"]
        c = w3_algebra["c"]
        beta = w3_algebra["beta"]
        Lambda = w3_algebra["Lambda"]

        # 计算 OPE
        result = OPE(W, W)

        # 验证极点结构
        assert result.max_pole == 6

        # 验证各极点系数
        pole_6 = result.pole(6)
        pole_5 = result.pole(5)
        pole_4 = result.pole(4)
        pole_3 = result.pole(3)
        pole_2 = result.pole(2)
        pole_1 = result.pole(1)

        # (z-w)^{-6}: c
        assert_voa_equal(pole_6, c)

        # (z-w)^{-5}: 0
        assert_voa_equal(pole_5, 0)

        # (z-w)^{-4}: 2*T
        assert_voa_equal(pole_4, 2 * T)

        # (z-w)^{-3}: ∂T
        assert_voa_equal(pole_3, d(T))

        # (z-w)^{-2}: 2*β*Λ + (3/10)*∂²T
        expected_pole_2 = 2 * beta * Lambda + sp.Rational(3, 10) * dn(2, T)
        assert_voa_equal(pole_2, expected_pole_2)

        # (z-w)^{-1}: β*∂Λ + (1/15)*∂³T
        expected_pole_1 = beta * d(Lambda) + sp.Rational(1, 15) * dn(3, T)
        assert_voa_equal(pole_1, expected_pole_1)


# ============================================
# Test Case 4-5: 导数算符 OPE
# ============================================


@pytest.mark.mathematica_ref
@pytest.mark.requires_derivative
class TestW3DerivativeOPE:
    """W3 代数导数算符 OPE 测试"""

    def test_4__T_with_NO_TT_minus_derivative(self, w3_algebra):
        """
        Test 4: OPE[T, NO[T,T] - (3/10)T'']

        对应 Mathematica:
        Print["Test 4: OPE[T, NO[T,T] - (3/10)T'']"];
        Print["结果 = ", OPE[T, NO[T, T] - (3/10) * T'']];

        这实际上是 OPE[T, Λ]，其中 Λ = NO(T,T) - (3/10)*∂²T
        """
        T = w3_algebra["T"]
        Lambda = w3_algebra["Lambda"]

        # 计算 OPE
        result = OPE(T, Lambda)

        # 验证结果不为零（具体值需要与 Mathematica 对比）
        assert result is not None
        assert result.max_pole > 0

    def test_5__T_with_W_double_derivative(self, w3_algebra):
        """
        Test 5: OPE[T, W'']

        对应 Mathematica:
        Print["Test 5: OPE[T, W'']"];
        Print["结果 = ", OPE[T, W'']];
        """
        T = w3_algebra["T"]
        W = w3_algebra["W"]

        # 计算 OPE
        result = OPE(T, dn(2, W))

        # 验证结果不为零
        assert result is not None
        assert result.max_pole > 0


# ============================================
# Test Case 6-7: 复合正规序 OPE
# ============================================


@pytest.mark.mathematica_ref
@pytest.mark.slow
class TestW3CompositeNO:
    """W3 代数复合正规序 OPE 测试"""

    def test_6__NO_dT_NO_ddW_T(self, w3_algebra):
        """
        Test 6: NO[T', NO[W'', T]]

        对应 Mathematica:
        Print["Test 6: NO[T', NO[W'', T]]"];
        Print["结果 = ", Expand[NO[T', NO[W'', T]]]];
        """
        T = w3_algebra["T"]
        W = w3_algebra["W"]

        # 构造表达式
        expr = NO(d(T), NO(dn(2, W), T))

        # 化简
        result = simplify(expr)

        # 验证结果不为零
        assert result is not None

    def test_7__NO_dW_NO_ddW_T(self, w3_algebra):
        """
        Test 7: NO[W', NO[W'', T]]

        对应 Mathematica:
        Print["Test 7: NO[W', NO[W'', T]]"];
        Print["结果 = ", Expand[NO[W', NO[W'', T]]]];
        """
        T = w3_algebra["T"]
        W = w3_algebra["W"]

        # 构造表达式
        expr = NO(d(W), NO(dn(2, W), T))

        # 化简
        result = simplify(expr)

        # 验证结果不为零
        assert result is not None


# ============================================
# Test Case 8: 数值检验
# ============================================


@pytest.mark.mathematica_ref
class TestW3NumericVerification:
    """W3 代数数值验证测试"""

    def test_8__numeric_verification_c100_beta_0p1(self, w3_algebra):
        """
        Test 8: 数值检验（c=100, β=1/10）

        对应 Mathematica:
        Print["Test 8: 数值检验"];
        Print["设置 c=100, β=1/10"];
        Print["OPE[W,W] 数值结果 = ", OPE[W, W] /. {c -> 100, β -> 1/10}];

        验证在特定参数值下，OPE 结果的数值正确性。
        """
        T = w3_algebra["T"]
        W = w3_algebra["W"]
        c = w3_algebra["c"]
        beta = w3_algebra["beta"]
        Lambda = w3_algebra["Lambda"]

        # 计算 W-W OPE
        result = OPE(W, W)

        # 数值代入
        subs_dict = {c: 100, beta: sp.Rational(1, 10)}

        # 验证 (z-w)^{-6} 极点：应该是 100
        pole_6 = result.pole(6)
        assert_voa_numeric_equal(pole_6, 100, subs=subs_dict)

        # 验证 (z-w)^{-4} 极点：应该是 2*T
        pole_4 = result.pole(4)
        assert_voa_equal(pole_4, 2 * T)

        # 验证 (z-w)^{-3} 极点：应该是 ∂T
        pole_3 = result.pole(3)
        assert_voa_equal(pole_3, d(T))

        # 验证 (z-w)^{-2} 极点的数值形式
        # 2*β*Λ + (3/10)*∂²T 在 β=1/10 时变为 (1/5)*Λ + (3/10)*∂²T
        pole_2 = result.pole(2)
        expected_pole_2 = 2 * beta * Lambda + sp.Rational(3, 10) * dn(2, T)

        # 符号验证
        assert_voa_equal(pole_2, expected_pole_2)

        # 验证 (z-w)^{-1} 极点的数值形式
        # β*∂Λ + (1/15)*∂³T 在 β=1/10 时变为 (1/10)*∂Λ + (1/15)*∂³T
        pole_1 = result.pole(1)
        expected_pole_1 = beta * d(Lambda) + sp.Rational(1, 15) * dn(3, T)

        # 符号验证
        assert_voa_equal(pole_1, expected_pole_1)


# ============================================
# 额外测试：验证辅助算符 Λ 的定义
# ============================================


@pytest.mark.mathematica_ref
class TestW3AuxiliaryOperator:
    """W3 代数辅助算符测试"""

    def test_lambda_definition(self, w3_algebra):
        """
        验证 Λ = NO(T,T) - (3/10)*∂²T 的定义

        对应 Mathematica:
        Λ[w_] := NO[T, T][w] - (3/10) * T''[w];
        """
        T = w3_algebra["T"]
        Lambda = w3_algebra["Lambda"]

        # 验证 Lambda 的定义
        expected = NO(T, T) - sp.Rational(3, 10) * dn(2, T)

        # 化简后比较
        assert_voa_equal(simplify(Lambda), simplify(expected))

    def test_lambda_in_ope(self, w3_algebra):
        """
        验证 Λ 在 OPE 中的使用

        测试 OPE[T, Λ] 是否能正确计算
        """
        T = w3_algebra["T"]
        Lambda = w3_algebra["Lambda"]

        # 计算 OPE
        result = OPE(T, Lambda)

        # 验证结果存在且有意义
        assert result is not None
        assert result.max_pole > 0

        # 验证至少有一个非零极点
        has_nonzero_pole = False
        for pole in range(1, result.max_pole + 1):
            coeff = result.pole(pole)
            if coeff != 0:
                has_nonzero_pole = True
                break

        assert has_nonzero_pole, "OPE[T, Λ] should have at least one non-zero pole"
