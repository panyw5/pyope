"""
测试 simplify 函数的正规排序功能

验证两种正规排序：
1. 算符之间按标准顺序排列：NO(C, NO(A, B)) -> NO(A, NO(B, C)) + additional terms
2. NO 化成右嵌套：NO(NO(...), ...) -> NO(..., NO(..., NO(...)))
"""

import pytest
import sympy as sp
from pyope import BasisOperator, NO, OPE, simplify, d
from pyope.ope_data import OPEData
from pyope.registry import ope_registry


@pytest.fixture(autouse=True)
def reset_registry():
    """每个测试前重置 OPE 注册表"""
    ope_registry.clear()
    yield
    ope_registry.clear()


class TestCanonicalOrdering:
    """测试正规排序功能"""

    def test_right_nested_ordering_simple(self):
        """测试右嵌套的简单排序：NO(C, NO(A, B)) 其中 A < B < C"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, B] = OPEData({})
        OPE[A, C] = OPEData({})
        OPE[B, C] = OPEData({})

        expr = NO(C, NO(A, B))
        result = simplify(expr)

        expected = NO(A, NO(B, C))
        assert result == expected

    def test_right_nested_ordering_with_ope(self):
        """测试右嵌套排序（有非零 OPE）"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, B] = OPEData({1: C})
        OPE[A, C] = OPEData({})
        OPE[B, C] = OPEData({})

        expr = NO(C, NO(A, B))
        result = simplify(expr)

        assert isinstance(result, (sp.Add, type(NO(A, B))))

    def test_left_nested_to_right_nested(self):
        """测试左嵌套转换为右嵌套：NO(NO(A, B), C) -> NO(A, NO(B, C))"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, B] = OPEData({})
        OPE[A, C] = OPEData({})
        OPE[B, C] = OPEData({})

        expr = NO(NO(A, B), C)
        result = simplify(expr)

        expected = NO(A, NO(B, C))
        assert result == expected

    def test_left_nested_to_right_nested_with_ope(self):
        """测试左嵌套转换（有非零 OPE）"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, C] = OPEData({1: B})
        OPE[B, C] = OPEData({})

        expr = NO(NO(A, B), C)
        result = simplify(expr)

        assert isinstance(result, (sp.Add, type(NO(A, B))))

    def test_triple_nested_ordering(self):
        """测试三重嵌套的排序"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)
        D = BasisOperator("D", conformal_weight=1)

        OPE[A, B] = OPEData({})
        OPE[A, C] = OPEData({})
        OPE[A, D] = OPEData({})
        OPE[B, C] = OPEData({})
        OPE[B, D] = OPEData({})
        OPE[C, D] = OPEData({})

        expr = NO(NO(NO(A, B), C), D)
        result = simplify(expr)

        expected = NO(A, NO(B, NO(C, D)))
        assert result == expected

    def test_preserve_nested_structure_flag(self):
        """测试 preserve_nested_structure 参数"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, B] = OPEData({})
        OPE[A, C] = OPEData({})
        OPE[B, C] = OPEData({})

        expr = NO(NO(A, B), C)

        result_default = simplify(expr, preserve_nested_structure=False)
        assert result_default == NO(A, NO(B, C))

        result_preserve = simplify(expr, preserve_nested_structure=True)
        assert isinstance(result_preserve, type(NO(A, B)))

    def test_ordering_consistency(self):
        """测试排序的一致性"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, B] = OPEData({})
        OPE[A, C] = OPEData({})
        OPE[B, C] = OPEData({})

        expr1 = NO(NO(A, B), C)
        expr2 = NO(A, NO(B, C))

        result1 = simplify(expr1)
        result2 = simplify(expr2)

        assert result1 == result2


class TestInnerLayerOrdering:
    """测试内层排序的处理"""

    def test_inner_layer_needs_reordering(self):
        """测试内层需要重排序的情况"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, B] = OPEData({})
        OPE[A, C] = OPEData({})
        OPE[B, C] = OPEData({})

        expr = NO(A, NO(C, B))
        result = simplify(expr)

        expected = NO(A, NO(B, C))
        assert result == expected

    def test_both_layers_need_reordering(self):
        """测试内外层都需要重排序"""
        A = BasisOperator("A", conformal_weight=1)
        B = BasisOperator("B", conformal_weight=1)
        C = BasisOperator("C", conformal_weight=1)

        OPE[A, B] = OPEData({})
        OPE[A, C] = OPEData({})
        OPE[B, C] = OPEData({})

        expr = NO(C, NO(B, A))
        result = simplify(expr)

        expected = NO(A, NO(B, C))
        assert result == expected


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
