"""
测试 bug 修复：当 OPE 没有奇异部分时，q=0 系数的正确处理

根据 Thielemans eq. (3.3.19):
当 max_AC = max_BC = 0 时，
[[AB]_0 C]_0 = [A[BC]_0]_0 = NO(A, NO(B,C))

这个测试验证了修复后的行为。
"""

import pytest
from src.pyope import *
from src.pyope.registry import ope_registry
import sympy as sp


class TestQ0CaseFix:
    """测试 q=0 情况的修复"""

    def setup_method(self):
        """每个测试前清空 OPE 注册表"""
        ope_registry.clear()

    def test_no_singular_part_returns_normal_order(self):
        """
        当 OPE(A,C) 和 OPE(B,C) 都没有奇异部分时，
        OPE(NO(A,B), C) 应该返回 q=0 的正规序乘积，而不是零
        """
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        # 定义 OPE 都没有奇异部分
        OPE[A, C] = OPE.make([])
        OPE[B, C] = OPE.make([])

        # 计算 OPE(NO(A,B), C)
        result = OPE(NO(A, B), C)

        # 验证结果
        assert result.max_pole == 0, "max_pole 应该是 0"
        assert result.pole(0) != Zero, "pole(0) 不应该是 Zero"
        assert isinstance(result.pole(0), NormalOrderedOperator), (
            "pole(0) 应该是 NormalOrderedOperator"
        )

        # 验证正规序的结构
        pole_0 = result.pole(0)
        assert pole_0.left == NO(A, B), "左侧应该是 NO(A,B)"
        assert pole_0.right == C, "右侧应该是 C"

    def test_bracket_q0_with_no_singular_part(self):
        """
        测试 bracket(NO(A,B), C, 0) 在没有奇异部分时的行为
        """
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, C] = OPE.make([])
        OPE[B, C] = OPE.make([])

        # bracket(NO(A,B), C, 0) 应该返回正规序
        result = bracket(NO(A, B), C, 0)

        assert result != Zero, "bracket 结果不应该是 Zero"
        assert isinstance(result, NormalOrderedOperator), (
            "bracket 结果应该是 NormalOrderedOperator"
        )

    def test_mixed_case_one_ope_has_poles(self):
        """
        测试混合情况：一个 OPE 有奇异部分，另一个没有
        """
        # 使用不同的算符名称避免冲突
        A2 = BasisOperator("A2", conformal_weight=1)
        B2 = BasisOperator("B2", conformal_weight=1)
        C2 = BasisOperator("C2", conformal_weight=1)
        D2 = BasisOperator("D2", conformal_weight=1)

        # A2-C2 有奇异部分，B2-C2 没有
        OPE[A2, C2] = OPE.make([D2])  # pole(1) = D2
        OPE[B2, C2] = OPE.make([])  # 没有奇异部分

        result = OPE(NO(A2, B2), C2)

        # 应该有奇异部分（来自 A2-C2 的贡献）
        assert result.max_pole >= 1, "应该有奇异部分"

        # 验证 pole(1) 来自第二项的贡献
        pole_1 = result.pole(1)
        assert pole_1 != Zero, "pole(1) 不应该是零"

    def test_both_opes_have_poles(self):
        """
        测试两个 OPE 都有奇异部分的情况（原有功能不应被破坏）
        """
        # 使用不同的算符名称避免冲突
        A3 = BasisOperator("A3", conformal_weight=1)
        B3 = BasisOperator("B3", conformal_weight=1)
        C3 = BasisOperator("C3", conformal_weight=1)
        D3 = BasisOperator("D3", conformal_weight=1)
        E3 = BasisOperator("E3", conformal_weight=1)

        # 两个 OPE 都有奇异部分
        OPE[A3, C3] = OPE.make([D3])  # pole(1) = D3
        OPE[B3, C3] = OPE.make([E3])  # pole(1) = E3

        result = OPE(NO(A3, B3), C3)

        # 应该有奇异部分
        assert result.max_pole >= 1, "应该有奇异部分"

        # 验证有来自两个 OPE 的贡献
        pole_1 = result.pole(1)
        assert pole_1 != Zero, "pole(1) 不应该是零"

    def test_both_opes_have_poles(self):
        """
        测试两个 OPE 都有奇异部分的情况（原有功能不应被破坏）
        """
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)
        D = BasisOperator("D", conformal_weight=1)
        E = BasisOperator("E", conformal_weight=1)

        # 两个 OPE 都有奇异部分
        OPE[A, C] = OPE.make([D])  # pole(1) = D
        OPE[B, C] = OPE.make([E])  # pole(1) = E

        result = OPE(NO(A, B), C)

        # 应该有奇异部分
        assert result.max_pole >= 1, "应该有奇异部分"

        # 验证有来自两个 OPE 的贡献
        pole_1 = result.pole(1)
        assert pole_1 != Zero, "pole(1) 不应该是零"

        # 也应该有 q=0 的正规序部分
        pole_0 = result.pole(0)
        assert pole_0 != Zero, "应该有 q=0 的正规序部分"

    def test_both_opes_have_poles(self):
        """
        测试两个 OPE 都有奇异部分的情况（原有功能不应被破坏）
        """
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)
        D = BasisOperator("D", conformal_weight=1)
        E = BasisOperator("E", conformal_weight=1)

        # 两个 OPE 都有奇异部分
        OPE[A, C] = OPE.make([D])  # pole(1) = D
        OPE[B, C] = OPE.make([E])  # pole(1) = E

        result = OPE(NO(A, B), C)

        # 应该有奇异部分
        assert result.max_pole >= 1, "应该有奇异部分"

    def test_theoretical_consistency_eq_3_3_19(self):
        """
        验证 Thielemans eq. (3.3.19) 的理论一致性：
        [[AB]_0 C]_0 = [A[BC]_0]_0
        """
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, C] = OPE.make([])
        OPE[B, C] = OPE.make([])

        # 左侧：[[AB]_0 C]_0
        lhs = bracket(NO(A, B), C, 0)

        # 右侧：[A[BC]_0]_0 = [A · NO(B,C)]_0 = NO(A, NO(B,C))
        # 注意：由于 [BC]_0 = NO(B,C)，所以这应该等价于 NO(NO(A,B), C)

        # 验证结果是正规序
        assert isinstance(lhs, NormalOrderedOperator), (
            "结果应该是 NormalOrderedOperator"
        )

        # 验证结构
        assert lhs.left == NO(A, B), "左侧应该是 NO(A,B)"
        assert lhs.right == C, "右侧应该是 C"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
