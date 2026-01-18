"""
工具函数模块

本模块提供各种辅助工具函数。
"""

import sympy as sp
from typing import Any


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
    # 导入算符类型（避免循环导入）
    from .operators import Operator
    from .constants import ConstantOperator

    # 如果已经是算符，直接返回
    if isinstance(coeff, (Operator, ConstantOperator)):
        return coeff

    # 如果已经是 SymPy 表达式，递归转换其中的浮点数
    if isinstance(coeff, sp.Expr):
        # 如果是 SymPy 的 Float，尝试转换为 Rational
        if isinstance(coeff, sp.Float):
            rational = sp.nsimplify(coeff)
            if isinstance(rational, (sp.Rational, sp.Integer)):
                return rational
            return coeff

        # 如果是 Mul 或 Add，递归转换其参数
        if isinstance(coeff, (sp.Mul, sp.Add)):
            # 使用 nsimplify 转换表达式中的所有浮点数
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
    except:
        return coeff

    # 如果已经是 SymPy 表达式，检查是否需要转换
    if isinstance(coeff, sp.Expr):
        # 如果是 SymPy 的 Float，尝试转换为 Rational
        if isinstance(coeff, sp.Float):
            # 尝试将浮点数转换为有理数
            rational = sp.nsimplify(coeff)
            # 如果转换后是 Rational 或 Integer，使用转换结果
            if isinstance(rational, (sp.Rational, sp.Integer)):
                return rational
        return coeff

    # 如果是 Python 的 int，转换为 SymPy Integer
    if isinstance(coeff, int):
        return sp.Integer(coeff)

    # 如果是 Python 的 float，尝试转换为 Rational
    if isinstance(coeff, float):
        # 使用 nsimplify 尝试找到精确的有理数表示
        rational = sp.nsimplify(coeff)
        # 如果转换成功（得到 Rational 或 Integer），使用转换结果
        if isinstance(rational, (sp.Rational, sp.Integer)):
            return rational
        # 否则使用 Float（保持精度）
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
    except:
        # 如果 sympify 失败，返回原值
        return coeff
