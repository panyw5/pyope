"""
工具函数模块

本模块提供各种辅助工具函数。
"""

import sympy as sp
from typing import Any


def _ascending_pochhammer(value: Any, order: int) -> sp.Expr:
    """Compute ascending Pochhammer symbol (value)_order = value*(value+1)*...*(value+order-1)."""
    result = sp.Integer(1)
    base = sp.sympify(value)
    for k in range(order):
        result *= base + k
    return sp.simplify(result)


def sympify_coefficient(coeff: Any) -> Any:
    """
    将系数自动转换为 SymPy 类型

    这个函数确保用户输入的数字（如 3/2, 13.0 等）被转换为 SymPy 的有理数或整数，
    而不是 Python 的浮点数。这样可以保持精确的符号计算。

    转换规则：
    - Python int -> sympy.Integer
    - Python float -> sympy.Rational (如果可以精确表示) 或 sympy.Float
    - Python Fraction -> sympy.Rational
    - SymPy 表达式中的 Float -> Rational (递归转换)
    - 算符对象 -> 保持不变

    Args:
        coeff: 要转换的系数

    Returns:
        转换后的系数

    Examples:
        >>> sympify_coefficient(3)
        3
        >>> sympify_coefficient(3.0)
        3
        >>> sympify_coefficient(0.5)
        1/2
        >>> sympify_coefficient(13.0)
        13
        >>> sympify_coefficient(13.0 * One)
        13*One
    """
    from .operators import Operator
    from .constants import ConstantOperator

    if isinstance(coeff, (Operator, ConstantOperator)):
        return coeff

    if isinstance(coeff, sp.Expr):
        if isinstance(coeff, sp.Float):
            rational = sp.nsimplify(coeff)
            if isinstance(rational, (sp.Rational, sp.Integer)):
                return rational
            return coeff

        if isinstance(coeff, (sp.Mul, sp.Add)):
            simplified = sp.nsimplify(coeff)
            return simplified

        return coeff

    # 如果是 Python 的 int，转换为 SymPy Integer
    if isinstance(coeff, int):
        return sp.Integer(coeff)

    # 如果是 Python 的 float，尝试转换为 Rational
    if isinstance(coeff, float):
        rational = sp.nsimplify(coeff)
        if isinstance(rational, (sp.Rational, sp.Integer)):
            return rational
        return sp.Float(coeff)

    # 如果是 Python 的 Fraction，转换为 SymPy Rational
    try:
        from fractions import Fraction

        if isinstance(coeff, Fraction):
            return sp.Rational(coeff.numerator, coeff.denominator)
    except ImportError:
        pass

    # 其他情况，尝试使用 sympify
    try:
        return sp.sympify(coeff)
    except (TypeError, ValueError, sp.SympifyError):
        return coeff

def _multiply_bracket_operator(bracket: Any, operator: Any) -> Any:
    """
    计算 bracket * operator，并移除结果中的 One

    当 bracket 是一个包含 One 的 sympy 表达式（如 -One = sp.Mul(-1, One)）时，
    bracket * operator 会产生 sp.Mul(-1, One, operator)，这是不正确的。

    这个函数会移除结果中的 One，返回正确的表达式。

    Args:
        bracket: bracket 系数（可能包含 One）
        operator: 算符

    Returns:
        简化后的表达式
    """
    from .constants import One

    # 如果 bracket 是 One，直接返回 operator
    if bracket is One or bracket == One:
        return operator

    # 计算乘积
    result = bracket * operator

    # 如果结果是 Mul，检查是否包含 One
    if isinstance(result, sp.Mul):
        # 提取所有不是 One 的因子
        new_args = []
        for arg in result.args:
            if arg is not One and arg != One:
                new_args.append(arg)

        # 如果移除了 One，重新构造表达式
        if len(new_args) < len(result.args):
            if len(new_args) == 0:
                return operator
            elif len(new_args) == 1:
                return new_args[0]
            else:
                return sp.Mul(*new_args)

    return result
