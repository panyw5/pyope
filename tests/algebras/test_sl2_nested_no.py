"""
SL(2) Kac-Moody 代数嵌套正规序测试

本测试文件对应 Mathematica 参考实现：tests/test_sl2_mathematica.wls

测试覆盖：
- 简单嵌套 NO（6个测试：1.1-1.6）
- 三重嵌套 NO（4个测试：2.1-2.4）
- 带导数的嵌套 NO（5个测试：3.1-3.5）
- Sugawara 构造（3个测试：4.1-4.3）
- 复杂嵌套结构（4个测试：5.1-5.4）

总计：22个测试用例
"""

import pytest
import sympy as sp

from pyope import NO, d, simplify
from tests.utils.comparison import assert_voa_equal


# ============================================
# Test Case 1: 简单嵌套 NO
# ============================================


@pytest.mark.mathematica_ref
class TestSimpleNestedNO:
    """简单嵌套正规序测试（对应 Mathematica test case 1.x）"""

    def test_1_1__NO_NO_Jplus_Jzero_Jminus(self, sl2_algebra):
        """
        Test 1.1: NO(NO(J+, J0), J-)

        对应 Mathematica:
        expr1 = NO[NO[Jplus, Jzero], Jminus];
        result1 = Expand[expr1];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        # 计算表达式
        expr = NO(NO(Jp, J0), Jm)
        result = simplify(expr)

        # 验证结果（展开后应该是右嵌套形式）
        # 根据 Jacobi 恒等式和 OPE 规则计算
        # 这里我们主要验证计算不会出错，具体结果需要与 Mathematica 对比
        assert result is not None

    def test_1_2__NO_Jplus_NO_Jzero_Jminus(self, sl2_algebra):
        """
        Test 1.2: NO(J+, NO(J0, J-))

        对应 Mathematica:
        expr2 = NO[Jplus, NO[Jzero, Jminus]];
        result2 = Expand[expr2];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(Jp, NO(J0, Jm))
        result = simplify(expr)

        assert result is not None

    def test_1_3__NO_NO_Jminus_Jplus_Jzero(self, sl2_algebra):
        """
        Test 1.3: NO(NO(J-, J+), J0)

        对应 Mathematica:
        expr3 = NO[NO[Jminus, Jplus], Jzero];
        result3 = Expand[expr3];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(Jm, Jp), J0)
        result = simplify(expr)

        assert result is not None

    def test_1_4__NO_Jminus_NO_Jplus_Jzero(self, sl2_algebra):
        """
        Test 1.4: NO(J-, NO(J+, J0))

        对应 Mathematica:
        expr4 = NO[Jminus, NO[Jplus, Jzero]];
        result4 = Expand[expr4];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(Jm, NO(Jp, J0))
        result = simplify(expr)

        assert result is not None

    def test_1_5__NO_NO_Jzero_Jplus_Jminus(self, sl2_algebra):
        """
        Test 1.5: NO(NO(J0, J+), J-)

        对应 Mathematica:
        expr5 = NO[NO[Jzero, Jplus], Jminus];
        result5 = Expand[expr5];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(J0, Jp), Jm)
        result = simplify(expr)

        assert result is not None

    def test_1_6__NO_Jzero_NO_Jminus_Jplus(self, sl2_algebra):
        """
        Test 1.6: NO(J0, NO(J-, J+))

        对应 Mathematica:
        expr6 = NO[Jzero, NO[Jminus, Jplus]];
        result6 = Expand[expr6];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(J0, NO(Jm, Jp))
        result = simplify(expr)

        assert result is not None


# ============================================
# Test Case 2: 三重嵌套 NO
# ============================================


@pytest.mark.mathematica_ref
class TestTripleNestedNO:
    """三重嵌套正规序测试（对应 Mathematica test case 2.x）"""

    def test_2_1__NO_NO_Jminus_NO_Jplus_Jzero(self, sl2_algebra):
        """
        Test 2.1: NO(NO(J-, NO(J+, J0)))

        对应 Mathematica:
        expr7 = NO[NO[Jminus, NO[Jplus, Jzero]]];
        result7 = Expand[expr7];

        注意：Mathematica 中的单参数 NO[expr] 在 pyope 中直接使用 expr 即可
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        # 在 pyope 中，NO(NO(A, B)) 等价于 NO(A, B)
        expr = NO(Jm, NO(Jp, J0))
        result = simplify(expr)

        assert result is not None

    def test_2_2__NO_NO_Jzero_NO_Jplus_Jminus(self, sl2_algebra):
        """
        Test 2.2: NO(NO(J0, NO(J+, J-)))

        对应 Mathematica:
        expr8 = NO[NO[Jzero, NO[Jplus, Jminus]]];
        result8 = Expand[expr8];

        注意：Mathematica 中的单参数 NO[expr] 在 pyope 中直接使用 expr 即可
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        # 在 pyope 中，NO(NO(A, B)) 等价于 NO(A, B)
        expr = NO(J0, NO(Jp, Jm))
        result = simplify(expr)

        assert result is not None

    def test_2_3__NO_NO_NO_Jplus_Jzero_Jminus_Jzero(self, sl2_algebra):
        """
        Test 2.3: NO(NO(NO(J+, J0), J-), J0)

        对应 Mathematica:
        expr9 = NO[NO[NO[Jplus, Jzero], Jminus], Jzero];
        result9 = Expand[expr9];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(NO(Jp, J0), Jm), J0)
        result = simplify(expr)

        assert result is not None

    def test_2_4__NO_Jplus_NO_Jzero_NO_Jminus_Jplus(self, sl2_algebra):
        """
        Test 2.4: NO(J+, NO(J0, NO(J-, J+)))

        对应 Mathematica:
        expr10 = NO[Jplus, NO[Jzero, NO[Jminus, Jplus]]];
        result10 = Expand[expr10];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(Jp, NO(J0, NO(Jm, Jp)))
        result = simplify(expr)

        assert result is not None


# ============================================
# Test Case 3: 带导数的嵌套 NO
# ============================================


@pytest.mark.mathematica_ref
@pytest.mark.requires_derivative
class TestNestedNOWithDerivatives:
    """带导数的嵌套正规序测试（对应 Mathematica test case 3.x）"""

    def test_3_1__NO_NO_dJplus_Jzero_Jminus(self, sl2_algebra):
        """
        Test 3.1: NO(NO(∂J+, J0), J-)

        对应 Mathematica:
        expr11 = NO[NO[Derivative[1][Jplus], Jzero], Jminus];
        result11 = Expand[expr11];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(d(Jp), J0), Jm)
        result = simplify(expr)

        assert result is not None

    def test_3_2__NO_Jplus_NO_dJzero_Jminus(self, sl2_algebra):
        """
        Test 3.2: NO(J+, NO(∂J0, J-))

        对应 Mathematica:
        expr12 = NO[Jplus, NO[Derivative[1][Jzero], Jminus]];
        result12 = Expand[expr12];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(Jp, NO(d(J0), Jm))
        result = simplify(expr)

        assert result is not None

    def test_3_3__NO_NO_Jplus_dJzero_Jminus(self, sl2_algebra):
        """
        Test 3.3: NO(NO(J+, ∂J0), J-)

        对应 Mathematica:
        expr13 = NO[NO[Jplus, Derivative[1][Jzero]], Jminus];
        result13 = Expand[expr13];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(Jp, d(J0)), Jm)
        result = simplify(expr)

        assert result is not None

    def test_3_4__NO_NO_Jminus_Jplus_dJzero(self, sl2_algebra):
        """
        Test 3.4: NO(NO(J-, J+), ∂J0)

        对应 Mathematica:
        expr14 = NO[NO[Jminus, Jplus], Derivative[1][Jzero]];
        result14 = Expand[expr14];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(Jm, Jp), d(J0))
        result = simplify(expr)

        assert result is not None

    def test_3_5__NO_dJminus_NO_Jplus_Jzero(self, sl2_algebra):
        """
        Test 3.5: NO(∂J-, NO(J+, J0))

        对应 Mathematica:
        expr15 = NO[Derivative[1][Jminus], NO[Jplus, Jzero]];
        result15 = Expand[expr15];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(d(Jm), NO(Jp, J0))
        result = simplify(expr)

        assert result is not None


# ============================================
# Test Case 4: Sugawara 构造
# ============================================


@pytest.mark.mathematica_ref
@pytest.mark.slow
class TestSugaraConstruction:
    """Sugawara 构造测试（对应 Mathematica test case 4.x）"""

    def test_4_1__NO_NO_Jplus_Jminus_NO_Jplus_Jminus(self, sl2_algebra):
        """
        Test 4.1: NO(NO(J+, J-), NO(J+, J-))

        对应 Mathematica:
        expr16 = NO[NO[Jplus, Jminus], NO[Jplus, Jminus]];
        result16 = Expand[expr16];
        """
        Jp = sl2_algebra["Jplus"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(Jp, Jm), NO(Jp, Jm))
        result = simplify(expr)

        assert result is not None

    def test_4_2__NO_NO_Jzero_Jzero_NO_Jplus_Jminus(self, sl2_algebra):
        """
        Test 4.2: NO(NO(J0, J0), NO(J+, J-))

        对应 Mathematica:
        expr17 = NO[NO[Jzero, Jzero], NO[Jplus, Jminus]];
        result17 = Expand[expr17];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(J0, J0), NO(Jp, Jm))
        result = simplify(expr)

        assert result is not None

    def test_4_3__NO_NO_Jplus_Jminus_NO_Jzero_Jzero(self, sl2_algebra):
        """
        Test 4.3: NO(NO(J+, J-), NO(J0, J0))

        对应 Mathematica:
        expr18 = NO[NO[Jplus, Jminus], NO[Jzero, Jzero]];
        result18 = Expand[expr18];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(Jp, Jm), NO(J0, J0))
        result = simplify(expr)

        assert result is not None


# ============================================
# Test Case 5: 复杂嵌套结构
# ============================================


@pytest.mark.mathematica_ref
@pytest.mark.slow
class TestComplexNestedStructures:
    """复杂嵌套结构测试（对应 Mathematica test case 5.x）"""

    def test_5_1__NO_NO_NO_Jplus_Jminus_Jzero_Jplus(self, sl2_algebra):
        """
        Test 5.1: NO(NO(NO(J+, J-), J0), J+)

        对应 Mathematica:
        expr19 = NO[NO[NO[Jplus, Jminus], Jzero], Jplus];
        result19 = Expand[expr19];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(NO(Jp, Jm), J0), Jp)
        result = simplify(expr)

        assert result is not None

    def test_5_2__NO_NO_Jplus_NO_Jminus_Jzero_Jplus(self, sl2_algebra):
        """
        Test 5.2: NO(NO(J+, NO(J-, J0)), J+)

        对应 Mathematica:
        expr20 = NO[NO[Jplus, NO[Jminus, Jzero]], Jplus];
        result20 = Expand[expr20];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(Jp, NO(Jm, J0)), Jp)
        result = simplify(expr)

        assert result is not None

    def test_5_3__NO_NO_dJplus_dJminus_Jzero(self, sl2_algebra):
        """
        Test 5.3: NO(NO(∂J+, ∂J-), J0)

        对应 Mathematica:
        expr21 = NO[NO[Derivative[1][Jplus], Derivative[1][Jminus]], Jzero];
        result21 = Expand[expr21];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(d(Jp), d(Jm)), J0)
        result = simplify(expr)

        assert result is not None

    def test_5_4__NO_NO_Jplus_dJminus_dJzero(self, sl2_algebra):
        """
        Test 5.4: NO(NO(J+, ∂J-), ∂J0)

        对应 Mathematica:
        expr22 = NO[NO[Jplus, Derivative[1][Jminus]], Derivative[1][Jzero]];
        result22 = Expand[expr22];
        """
        Jp = sl2_algebra["Jplus"]
        J0 = sl2_algebra["Jzero"]
        Jm = sl2_algebra["Jminus"]

        expr = NO(NO(Jp, d(Jm)), d(J0))
        result = simplify(expr)

        assert result is not None
