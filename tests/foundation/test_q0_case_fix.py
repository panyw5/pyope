"""
测试 OPEData 契约：OPE 只捕获奇异部分 (q >= 1)

根据系统设计契约：
- OPE(A, B) 返回奇异部分（q >= 1）
- 正规序（q = 0）通过 NO() 函数和 bracket(A, B, 0) 获取
- OPEData 只存储极点（poles），不包含正则部分

参考：simplify.py:512 注释：
"OPE(A,B) 只包含奇异部分（l>=1），0 阶由 NO 定义"
"""

import pytest
import sympy as sp

from pyope import *
from pyope.registry import ope_registry


class TestOPEDataContract:
    """测试 OPEData 的设计契约：只存储奇异部分"""

    def setup_method(self):
        """每个测试前清空 OPE 注册表"""
        ope_registry.clear()

    def test_no_singular_part_returns_zero(self):
        """
        当 OPE(A,C) 和 OPE(B,C) 都没有奇异部分时，
        OPE(NO(A,B), C) 应该返回空 OPEData（奇异部分为零）
        """
        A = BasicOperator("A", conformal_weight=1)
        B = BasicOperator("B", conformal_weight=1)
        C = BasicOperator("C", conformal_weight=1)

        # 定义 OPE 都没有奇异部分
        OPE[A, C] = OPE.make([])
        OPE[B, C] = OPE.make([])

        # 计算 OPE(NO(A,B), C)
        result = OPE(NO(A, B), C)

        # 验证结果：奇异部分为空
        assert result.is_zero(), "没有奇异部分时，OPEData 应该是零"
        assert result.max_pole == 0, "max_pole 应该是 0"

        # pole(0) 应该返回 Zero（因为 OPEData 不存储 q=0）
        assert result.pole(0) == Zero, "pole(0) 应该返回 Zero"

    def test_bracket_q0_returns_normal_order(self):
        """
        测试 bracket(NO(A,B), C, 0) 在没有奇异部分时的行为
        bracket(A, B, 0) 应该返回正规序，而不是从 OPEData.pole(0) 获取
        """
        A = BasicOperator("A", conformal_weight=1)
        B = BasicOperator("B", conformal_weight=1)
        C = BasicOperator("C", conformal_weight=1)

        OPE[A, C] = OPE.make([])
        OPE[B, C] = OPE.make([])

        # bracket(NO(A,B), C, 0) 应该返回正规序
        result = bracket(NO(A, B), C, 0)

        # 验证：结果是正规序，不是从 OPEData 获取
        assert result != Zero, "bracket 结果不应该是 Zero"
        assert isinstance(result, NormalOrderedOperator), (
            "bracket 结果应该是 NormalOrderedOperator"
        )

        # 验证结构：bracket(NO(A,B), C, 0) = NO(NO(A,B), C)
        assert result.left == NO(A, B), "左侧应该是 NO(A,B)"
        assert result.right == C, "右侧应该是 C"

    def test_mixed_case_one_ope_has_poles(self):
        """
        测试混合情况：一个 OPE 有奇异部分，另一个没有
        """
        # 使用不同的算符名称避免冲突
        A = BasicOperator("A", conformal_weight=1)
        B = BasicOperator("B", conformal_weight=1)
        C = BasicOperator("C", conformal_weight=1)
        D = BasicOperator("D", conformal_weight=1)

        # A-C 有奇异部分，B-C 没有
        OPE[A, C] = OPE.make([D])  # pole(1) = D
        OPE[B, C] = OPE.make([])  # 没有奇异部分

        result = OPE(NO(A, B), C)

        # 应该有奇异部分（来自 A-C 的贡献）
        assert result.max_pole >= 1, "应该有奇异部分"
        assert not result.is_zero(), "不应该为零"

        # 验证 pole(1) 来自 A-C 的贡献
        pole_1 = result.pole(1)
        assert pole_1 != Zero, "pole(1) 不应该是零"

        # 但 pole(0) 仍然是 Zero（OPEData 不存储 q=0）
        assert result.pole(0) == Zero, "pole(0) 应该是 Zero"

    def test_both_opes_have_poles(self):
        """
        测试两个 OPE 都有奇异部分的情况
        """
        A = BasicOperator("A", conformal_weight=1)
        B = BasicOperator("B", conformal_weight=1)
        C = BasicOperator("C", conformal_weight=1)
        D = BasicOperator("D", conformal_weight=1)
        E = BasicOperator("E", conformal_weight=1)

        # 两个 OPE 都有奇异部分
        OPE[A, C] = OPE.make([D])  # pole(1) = D
        OPE[B, C] = OPE.make([E])  # pole(1) = E

        result = OPE(NO(A, B), C)

        # 应该有奇异部分
        assert result.max_pole >= 1, "应该有奇异部分"
        assert not result.is_zero(), "不应该为零"

        # 验证有来自两个 OPE 的贡献
        pole_1 = result.pole(1)
        assert pole_1 != Zero, "pole(1) 不应该是零"

        # 但 pole(0) 仍然是 Zero（OPEData 不存储 q=0）
        assert result.pole(0) == Zero, "pole(0) 应该是 Zero"

    def test_theoretical_consistency_bracket_vs_ope(self):
        """
        验证理论一致性：
        1. OPE(NO(A,B), C) 返回奇异部分（这里为空）
        2. bracket(NO(A,B), C, 0) 返回正规序
        两者是分离的，不冲突
        """
        A = BasicOperator("A", conformal_weight=1)
        B = BasicOperator("B", conformal_weight=1)
        C = BasicOperator("C", conformal_weight=1)

        OPE[A, C] = OPE.make([])
        OPE[B, C] = OPE.make([])

        # OPE 返回奇异部分（为空）
        ope_result = OPE(NO(A, B), C)
        assert ope_result.is_zero(), "OPE 应该返回空的奇异部分"

        # bracket 返回正规序
        bracket_result = bracket(NO(A, B), C, 0)
        assert bracket_result != Zero, "bracket 应该返回正规序"
        assert isinstance(bracket_result, NormalOrderedOperator), (
            "bracket 结果应该是 NormalOrderedOperator"
        )

        # 验证结构
        assert bracket_result.left == NO(A, B), "左侧应该是 NO(A,B)"
        assert bracket_result.right == C, "右侧应该是 C"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
