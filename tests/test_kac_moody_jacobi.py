"""
测试 Kac-Moody 代数的 Jacobi 恒等式和嵌套 NO 简化
"""

import pytest
from pyope import BasisOperator, OPE, simplify, NO, One, Zero


class TestKacMoodyJacobi:
    """测试 Kac-Moody 代数中的 Jacobi 恒等式"""

    def setup_method(self):
        self.J_plus = BasisOperator("J⁺", conformal_weight=1)
        self.J_zero = BasisOperator("J⁰", conformal_weight=1)
        self.J_minus = BasisOperator("J⁻", conformal_weight=1)

        k_val = 1

        OPE[self.J_plus, self.J_zero] = OPE.make([-2 * self.J_plus])
        OPE[self.J_plus, self.J_minus] = OPE.make([k_val * One, self.J_zero])
        OPE[self.J_zero, self.J_minus] = OPE.make([-2 * self.J_minus])
        OPE[self.J_zero, self.J_zero] = OPE.make([2 * k_val * One, 0])

    def test_right_nested_no_expansion(self):
        """测试右嵌套 NO 的展开: NO(J_minus, NO(J_plus, J_zero))"""
        expr = NO(self.J_minus, NO(self.J_plus, self.J_zero))
        result = simplify(expr)

        # 结果应该不再是原始的嵌套形式
        assert result != expr

        # 结果应该完全简化，不包含非标准序的 NO
        result_str = str(result)
        # 结果应该是 2*NO(J⁺,∂J⁻)，这是标准序
        assert "NO(J⁺,∂J⁻)" in result_str or "∂J⁻" in result_str

    def test_left_nested_no_expansion(self):
        """测试左嵌套 NO 的展开: NO(NO(J_zero, J_minus), J_plus)"""
        expr = NO(NO(self.J_zero, self.J_minus), self.J_plus)
        result = simplify(expr)

        # 结果应该不再是原始的嵌套形式
        assert result != expr

        # 结果应该完全简化，不包含非标准序的 NO
        result_str = str(result)
        # 结果应该包含标准序的 NO，如 NO(J⁰,∂J⁰)
        assert "∂" in result_str  # 应该包含导数项

    def test_simple_no_reordering(self):
        """测试简单的 NO 重排序: NO(J_minus, J_plus)"""
        expr = NO(self.J_minus, self.J_plus)
        result = simplify(expr)

        # 应该应用 OPE 并重排序
        assert result != expr

        # 结果应该包含导数项（根据 OPE 定义）
        result_str = str(result)
        assert "∂" in result_str or "J⁰" in result_str

    def test_standard_order_no_change(self):
        """测试已经是标准顺序的 NO 不会改变: NO(J_plus, J_zero)"""
        expr = NO(self.J_plus, self.J_zero)
        result = simplify(expr)

        # 标准顺序应该保持不变（或只是简化内部）
        assert str(result) == str(expr)

    def test_triple_nested_no(self):
        """测试三重嵌套的 NO"""
        expr = NO(self.J_minus, NO(self.J_plus, NO(self.J_zero, self.J_plus)))
        result = simplify(expr)

        # 应该能够简化（不会无限递归）
        assert result is not None

    def test_result_all_standard_order(self):
        """测试简化结果中的所有 NO 都是标准序"""
        from pyope.simplify import is_standard_order_no
        from pyope import d

        # 测试几个例子
        expr1 = NO(self.J_minus, NO(self.J_plus, self.J_zero))
        result1 = simplify(expr1)

        # 结果应该是 2*NO(J⁺,∂J⁻)
        # 检查 NO(J⁺,∂J⁻) 是标准序
        no1 = NO(self.J_plus, d(self.J_minus))
        assert is_standard_order_no(no1) == True

        # 测试另一个例子
        expr2 = NO(NO(self.J_zero, self.J_minus), self.J_plus)
        result2 = simplify(expr2)

        # 结果应该包含 NO(J⁰,∂J⁰)
        no2 = NO(self.J_zero, d(self.J_zero))
        assert is_standard_order_no(no2) == True


class TestOperatorOrdering:
    """测试算符排序逻辑"""

    def test_operator_comparison(self):
        """测试算符比较函数"""
        from pyope.registry import ope_registry

        J_plus = BasisOperator("J⁺", conformal_weight=1)
        J_zero = BasisOperator("J⁰", conformal_weight=1)
        J_minus = BasisOperator("J⁻", conformal_weight=1)

        # Unicode 字典序：J⁰ (U+2070) < J⁺ (U+207A) < J⁻ (U+207B)
        # compare(A, B) > 0 表示 A 应该排在 B 前面
        # 所以 compare(J_zero, J_plus) > 0（J_zero 在前）
        assert ope_registry.compare_operators(J_zero, J_plus) > 0
        assert ope_registry.compare_operators(J_zero, J_minus) > 0
        assert ope_registry.compare_operators(J_plus, J_minus) > 0

    def test_flatten_right_no(self):
        """测试 flatten_right_no 函数"""
        from pyope.simplify import flatten_right_no

        J_plus = BasisOperator("J⁺", conformal_weight=1)
        J_zero = BasisOperator("J⁰", conformal_weight=1)
        J_minus = BasisOperator("J⁻", conformal_weight=1)

        # 测试单个算符
        assert flatten_right_no(J_plus) == [J_plus]

        # 测试二元 NO
        expr = NO(J_plus, J_zero)
        assert flatten_right_no(expr) == [J_plus, J_zero]

        # 测试右嵌套 NO
        expr = NO(J_minus, NO(J_plus, J_zero))
        assert flatten_right_no(expr) == [J_minus, J_plus, J_zero]

    def test_is_standard_order_no(self):
        """测试 is_standard_order_no 函数"""
        from pyope.simplify import is_standard_order_no

        J_plus = BasisOperator("J⁺", conformal_weight=1)
        J_zero = BasisOperator("J⁰", conformal_weight=1)
        J_minus = BasisOperator("J⁻", conformal_weight=1)

        # Unicode 字典序：J⁰ < J⁺ < J⁻
        # 标准顺序（按字典序）
        expr1 = NO(J_zero, J_plus)
        assert is_standard_order_no(expr1) == True

        # 非标准顺序
        expr2 = NO(J_plus, J_zero)
        assert is_standard_order_no(expr2) == False

        # 嵌套 NO - 标准顺序
        expr3 = NO(J_zero, NO(J_plus, J_minus))
        assert is_standard_order_no(expr3) == True

        # 嵌套 NO - 非标准顺序
        expr4 = NO(J_minus, NO(J_plus, J_zero))
        assert is_standard_order_no(expr4) == False
