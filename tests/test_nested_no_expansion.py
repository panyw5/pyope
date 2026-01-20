"""
测试嵌套 NO 的展开行为

验证 simplify 函数的 preserve_nested_structure 参数：
- 默认行为（False）：自动展开左嵌套为右嵌套
- preserve_nested_structure=True：仅在需要排序时展开
"""

import pytest
from pyope import BasisOperator, OPE, simplify, NO
from pyope.operators import NormalOrderedOperator


class TestNestedNOExpansion:
    """测试嵌套 NO 的展开行为"""

    def setup_method(self):
        """设置 Kac-Moody 代数"""
        self.J_plus = BasisOperator("J⁺", conformal_weight=1)
        self.J_zero = BasisOperator("J⁰", conformal_weight=1)
        self.J_minus = BasisOperator("J⁻", conformal_weight=1)

        k_val = 1
        OPE[self.J_plus, self.J_zero] = OPE.make([-2 * self.J_plus])
        OPE[self.J_plus, self.J_minus] = OPE.make([k_val, self.J_zero])
        OPE[self.J_zero, self.J_minus] = OPE.make([-2 * self.J_minus])
        OPE[self.J_zero, self.J_zero] = OPE.make([2 * k_val, 0])

    def test_default_expands_left_nesting(self):
        """默认行为：展开左嵌套 NO(NO(A,B),C) 为右嵌套"""
        expr = NO(NO(self.J_plus, self.J_zero), self.J_minus)
        result = simplify(expr)

        # 结果不应该是左嵌套形式
        assert not isinstance(result, NormalOrderedOperator) or not isinstance(
            result.left, NormalOrderedOperator
        ), "默认应该展开左嵌套"

        # 结果应该包含右嵌套形式 NO(J⁺, NO(J⁰, J⁻))
        result_str = str(result)
        assert "NO(J⁺,NO(J⁰,J⁻))" in result_str, (
            f"结果应包含右嵌套形式，实际: {result_str}"
        )

    def test_preserve_structure_keeps_left_nesting(self):
        """preserve_nested_structure=True：保持左嵌套结构"""
        expr = NO(NO(self.J_plus, self.J_zero), self.J_minus)
        result = simplify(expr, preserve_nested_structure=True)

        # 结果应该保持左嵌套形式
        assert isinstance(result, NormalOrderedOperator), (
            "结果应该是 NormalOrderedOperator"
        )
        assert isinstance(result.left, NormalOrderedOperator), "左侧应该保持嵌套结构"
        assert result.left.left == self.J_plus
        assert result.left.right == self.J_zero
        assert result.right == self.J_minus

    def test_both_modes_expand_when_reordering_needed(self):
        """两种模式在需要重排序时都应该展开"""
        expr = NO(NO(self.J_zero, self.J_minus), self.J_plus)

        # 默认模式
        result_default = simplify(expr)
        # 保持结构模式
        result_preserve = simplify(expr, preserve_nested_structure=True)

        # 两者都应该展开（因为需要重排序）
        assert str(result_default) == str(result_preserve), (
            "需要重排序时，两种模式应该产生相同结果"
        )

        # 结果应该包含右嵌套形式
        result_str = str(result_default)
        assert "NO(J⁰,NO(J⁻,J⁺))" in result_str

    def test_default_expands_without_ope(self):
        """默认行为：即使没有 OPE 也展开左嵌套"""
        # 创建没有定义 OPE 的算符
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        expr = NO(NO(A, B), C)
        result = simplify(expr)

        # 应该转换为右嵌套形式（即使没有 OPE）
        assert isinstance(result, NormalOrderedOperator)
        # 检查是否是右嵌套形式 NO(A, NO(B, C))
        if isinstance(result.right, NormalOrderedOperator):
            assert result.left == A
            assert result.right.left == B
            assert result.right.right == C

    def test_preserve_structure_without_ope(self):
        """preserve_nested_structure=True：没有 OPE 时保持结构"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        expr = NO(NO(A, B), C)
        result = simplify(expr, preserve_nested_structure=True)

        # 应该保持左嵌套形式
        assert isinstance(result, NormalOrderedOperator)
        assert isinstance(result.left, NormalOrderedOperator)
        assert result.left.left == A
        assert result.left.right == B
        assert result.right == C

    def test_expansion_produces_correction_terms(self):
        """展开应该产生校正项（来自 Jacobi 恒等式）"""
        expr = NO(NO(self.J_plus, self.J_zero), self.J_minus)
        result = simplify(expr)

        result_str = str(result)

        # 应该包含导数项（校正项）
        assert "∂" in result_str, "应该包含导数校正项"

        # 应该包含多个项（不只是简单的重排）
        assert "+" in result_str or "-" in result_str, "应该包含多个项（主项 + 校正项）"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
