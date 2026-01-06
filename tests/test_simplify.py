"""
测试 simplify 模块
"""

import pytest
import sympy as sp
from sympy import Symbol

from pyope import (
    BasisOperator,
    d,
    One,
    Zero,
    OPE,
    NO,
    MakeOPE,
    Bosonic,
    simplify,
    canonicalize,
    collect_normal_ordered_terms,
)


@pytest.fixture(autouse=True)
def clear_registry():
    """每个测试前清空注册表"""
    from pyope.registry import ope_registry
    ope_registry.clear()
    yield
    ope_registry.clear()


class TestSimplifyBasic:
    """测试基本化简功能"""

    def test_simplify_zero(self):
        """测试零的化简"""
        assert simplify(0) == 0
        assert simplify(Zero) == 0

    def test_simplify_one(self):
        """测试单位的化简"""
        assert simplify(One) == One

    def test_simplify_basis_operator(self):
        """测试基本算符的化简"""
        T = BasisOperator("T", bosonic=True)
        Bosonic(T)

        assert simplify(T) == T

    def test_simplify_derivative(self):
        """测试导数算符的化简"""
        T = BasisOperator("T", bosonic=True)
        Bosonic(T)

        dT = d(T)
        assert simplify(dT) == dT

    def test_simplify_scalar_multiplication(self):
        """测试标量乘法的化简"""
        T = BasisOperator("T", bosonic=True)
        Bosonic(T)

        expr = 2 * T
        assert simplify(expr) == 2 * T

    def test_simplify_addition(self):
        """测试加法的化简"""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        expr = T + J
        result = simplify(expr)
        # 应该保持原样
        assert result == T + J


class TestSimplifyNormalOrdered:
    """测试正规序化简"""

    def test_simplify_simple_no(self):
        """测试简单正规序的化简"""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        no_TJ = NO(T, J)
        result = simplify(no_TJ)

        # 简单 NO 保持不变
        assert result == NO(T, J)

    def test_simplify_no_with_derivatives(self):
        """测试包含导数的正规序"""
        T = BasisOperator("T", bosonic=True)
        Bosonic(T)

        no_TdT = NO(T, d(T))
        result = simplify(no_TdT)

        # 应该保持原样
        assert result == NO(T, d(T))

    def test_simplify_no_with_scalar(self):
        """测试标量乘以正规序"""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        expr = 3 * NO(T, J)
        result = simplify(expr)

        assert result == 3 * NO(T, J)

    def test_simplify_sum_of_no(self):
        """测试正规序的和"""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        expr = NO(T, J) + NO(J, T)
        result = simplify(expr)

        # 应该保持为和的形式
        assert isinstance(result, sp.Add)


class TestSimplifyOPEData:
    """测试 OPEData 的化简"""

    def test_simplify_ope_data(self):
        """测试 OPEData 对象的化简"""
        from pyope import OPEData

        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        # 创建 OPEData
        ope_data = OPEData({
            2: 2 * T + J,
            1: d(T)
        })

        result = simplify(ope_data)

        # 应该返回化简后的 OPEData
        assert isinstance(result, OPEData)
        assert result.max_pole == 2


class TestCollectNormalOrderedTerms:
    """测试正规序项的收集"""

    def test_collect_simple_terms(self):
        """测试收集简单项"""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        expr = 2 * NO(T, J) + 3 * NO(T, J)
        terms = collect_normal_ordered_terms(expr)

        # 应该合并同类项
        assert len(terms) == 1
        key = ('NO', ('Basis', 'T', True), ('Basis', 'J', True))
        assert terms[key] == 5

    def test_collect_different_terms(self):
        """测试收集不同的项"""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        expr = NO(T, J) + NO(J, T) + NO(J, J)
        terms = collect_normal_ordered_terms(expr)

        # 应该有三个不同的键
        assert len(terms) == 3


class TestSimplifyIntegration:
    """测试与 OPE 计算的集成"""

    def test_simplify_ope_result(self):
        """测试化简 OPE 计算结果"""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        Bosonic(T, J)

        # 定义 OPE
        OPE[T, J] = MakeOPE([J, d(J)])

        # 计算 OPE
        result = OPE(T, J)

        # 化简
        simplified = simplify(result)

        # 应该保持结构
        assert simplified.max_pole == 2
        assert simplified.pole(2) == J
        assert simplified.pole(1) == d(J)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
