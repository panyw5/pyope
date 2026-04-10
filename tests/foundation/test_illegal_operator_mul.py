"""
测试非法算符乘积的检测

本测试文件验证所有入口点都能正确检测并拒绝非法的算符乘积。
"""

import pytest
import sympy as sp

from pyope import *
from pyope.exceptions import IllegalOperatorProductError


class TestOperatorMultiplication:
    """测试 Operator.__mul__ 和 __rmul__ 的禁止行为"""

    def test_operator_times_operator_raises(self):
        T = BasicOperator("T")
        J = BasicOperator("J")

        with pytest.raises(
            IllegalOperatorProductError, match="Illegal operator product"
        ):
            T * J

    def test_operator_times_derivative_raises(self):
        T = BasicOperator("T")

        with pytest.raises(
            IllegalOperatorProductError, match="Illegal operator product"
        ):
            T * d(T)

    def test_operator_times_normal_ordered_raises(self):
        T = BasicOperator("T")
        J = BasicOperator("J")

        with pytest.raises(
            IllegalOperatorProductError, match="Illegal operator product"
        ):
            T * NO(J, T)

    def test_scalar_times_operator_ok(self):
        T = BasicOperator("T")
        result = 2 * T
        assert result == 2 * T

    def test_symbol_times_operator_ok(self):
        T = BasicOperator("T")
        c = sp.Symbol("c")
        result = c * T
        assert result == c * T

    def test_scalar_times_operator_sum_ok(self):
        T = BasicOperator("T")
        J = BasicOperator("J")
        result = 2 * (T + J)
        assert result == 2 * T + 2 * J


class TestNOFunction:
    """测试 NO() 函数的输入验证"""

    def test_NO_rejects_mul_operator_operator_left(self):
        T = BasicOperator("T")
        J = BasicOperator("J")
        W = BasicOperator("W")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError):
            NO(mul_expr, W)

    def test_NO_rejects_mul_operator_operator_right(self):
        T = BasicOperator("T")
        J = BasicOperator("J")
        W = BasicOperator("W")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError):
            NO(W, mul_expr)

    def test_NO_accepts_scalar_mul(self):
        T = BasicOperator("T")
        J = BasicOperator("J")

        result = NO(2 * T, J)
        assert result == 2 * NO(T, J)


class TestOPEFunction:
    """测试 OPE() 函数的输入验证"""

    def test_OPE_rejects_mul_operator_operator_left(self):
        T = BasicOperator("T")
        J = BasicOperator("J")
        W = BasicOperator("W")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError, match="OPE\\(left\\)"):
            OPE(mul_expr, W)

    def test_OPE_rejects_mul_operator_operator_right(self):
        T = BasicOperator("T")
        J = BasicOperator("J")
        W = BasicOperator("W")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError, match="OPE\\(right\\)"):
            OPE(W, mul_expr)

    def test_OPE_accepts_scalar_mul(self):
        T = BasicOperator("T")
        OPE[T, T] = MakeOPE([2 * T])

        result = OPE(2 * T, T)
        assert result.max_pole == 1


class TestBracketFunction:
    """测试 bracket() 函数的输入验证"""

    def test_bracket_rejects_mul_operator_operator_left(self):
        T = BasicOperator("T")
        J = BasicOperator("J")
        W = BasicOperator("W")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError, match="bracket\\(left\\)"):
            bracket(mul_expr, W, 1)

    def test_bracket_rejects_mul_operator_operator_right(self):
        T = BasicOperator("T")
        J = BasicOperator("J")
        W = BasicOperator("W")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError, match="bracket\\(right\\)"):
            bracket(W, mul_expr, 1)


class TestMakeOPE:
    """测试 MakeOPE() 和 OPE.make() 的输入验证"""

    def test_MakeOPE_rejects_mul_operator_operator(self):
        T = BasicOperator("T")
        J = BasicOperator("J")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(
            IllegalOperatorProductError, match="MakeOPE\\(data\\[0\\]\\)"
        ):
            MakeOPE([mul_expr])

    def test_OPE_make_rejects_mul_operator_operator(self):
        T = BasicOperator("T")
        J = BasicOperator("J")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(
            IllegalOperatorProductError, match="OPE\\.make\\(data\\[0\\]\\)"
        ):
            OPE.make([mul_expr])

    def test_MakeOPE_accepts_scalar_mul(self):
        T = BasicOperator("T")

        result = MakeOPE([2 * T])
        assert result.pole(1) == 2 * T


class TestExtractScalarOperator:
    """测试 extract_scalar_operator() 的行为"""

    def test_extract_rejects_mul_operator_operator(self):
        from src.pyope.local_operator import extract_scalar_operator

        T = BasicOperator("T")
        J = BasicOperator("J")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(
            IllegalOperatorProductError, match="extract_scalar_operator"
        ):
            extract_scalar_operator(mul_expr)

    def test_extract_accepts_scalar_mul(self):
        from src.pyope.local_operator import extract_scalar_operator

        T = BasicOperator("T")
        coeff, op = extract_scalar_operator(2 * T)

        assert coeff == 2
        assert op == T

    def test_extract_handles_operator_sum_correctly(self):
        from src.pyope.local_operator import extract_scalar_operator

        T = BasicOperator("T")
        J = BasicOperator("J")

        expr = sp.Mul(2, sp.Add(T, J, evaluate=False), evaluate=False)
        coeff, op = extract_scalar_operator(expr)
        assert coeff == 2
        assert isinstance(op, sp.Add)
        assert set(op.args) == {T, J}


class TestAssertNoIllegalOperatorMul:
    """测试 assert_no_illegal_operator_mul() 函数"""

    def test_detects_mul_operator_operator(self):
        from src.pyope.local_operator import assert_no_illegal_operator_mul

        T = BasicOperator("T")
        J = BasicOperator("J")

        mul_expr = sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError):
            assert_no_illegal_operator_mul(mul_expr, context="test")

    def test_detects_nested_mul(self):
        from src.pyope.local_operator import assert_no_illegal_operator_mul

        T = BasicOperator("T")
        J = BasicOperator("J")

        nested = 2 * sp.Mul(T, J, evaluate=False)
        with pytest.raises(IllegalOperatorProductError):
            assert_no_illegal_operator_mul(nested, context="test")

    def test_accepts_scalar_mul(self):
        from src.pyope.local_operator import assert_no_illegal_operator_mul

        T = BasicOperator("T")

        assert_no_illegal_operator_mul(2 * T, context="test")

    def test_accepts_operator_sum(self):
        from src.pyope.local_operator import assert_no_illegal_operator_mul

        T = BasicOperator("T")
        J = BasicOperator("J")

        assert_no_illegal_operator_mul(T + J, context="test")
