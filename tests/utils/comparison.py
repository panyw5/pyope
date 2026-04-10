"""
表达式比较工具模块

提供用于比较 VOA 表达式的工具函数，支持符号比较和数值比较。
"""

from typing import Any, Dict, Optional, Union
import math

import sympy as sp
from sympy import Add, Mul, Number, Symbol

from pyope import (
    simplify,
    Zero,
    One,
    Operator,
    BasicOperator,
    DerivativeOperator,
    NormalOrderedOperator,
)
from pyope.local_operator import is_local_operator, extract_scalar_operator
from pyope.ope_data import OPEData


def canonicalize(expr: Any, *, expand: bool = True) -> Any:
    """
    将表达式规范化为标准形式

    Args:
        expr: 要规范化的表达式
        expand: 是否展开表达式（默认 True）

    Returns:
        规范化后的表达式

    Examples:
        >>> T = BasicOperator("T")
        >>> expr = NO(T, T) + 2*T - T
        >>> canonicalize(expr)  # 返回 NO(T, T) + T
    """
    # 处理零
    if expr == 0 or expr == Zero:
        return Zero

    # 处理单位
    if expr == 1 or expr == One:
        return One

    # 处理 One * expr 的情况（将 One*c 转换为 c）
    if isinstance(expr, Mul):
        args = list(expr.args)
        # 移除 One
        args = [arg for arg in args if arg != One and arg != 1]
        if len(args) == 0:
            return One
        elif len(args) == 1:
            expr = args[0]
        else:
            expr = Mul(*args)

    # 使用 pyope 的 simplify 函数
    if expand:
        try:
            result = simplify(expr)
            # 再次处理 One * expr
            if isinstance(result, Mul):
                args = list(result.args)
                args = [arg for arg in args if arg != One and arg != 1]
                if len(args) == 0:
                    return One
                elif len(args) == 1:
                    return args[0]
                else:
                    return Mul(*args)
            return result
        except Exception:
            # 如果 simplify 失败，返回原表达式
            pass

    # 对于 sympy 表达式，尝试展开和化简
    if isinstance(expr, (Add, Mul)):
        try:
            expanded = sp.expand(expr)
            # 处理 One
            if isinstance(expanded, Mul):
                args = list(expanded.args)
                args = [arg for arg in args if arg != One and arg != 1]
                if len(args) == 0:
                    return One
                elif len(args) == 1:
                    return args[0]
                else:
                    return Mul(*args)
            return expanded
        except Exception:
            pass

    return expr


def _extract_terms(expr: Any) -> list:
    """
    提取表达式中的所有项（用于加法）

    Args:
        expr: 表达式

    Returns:
        项的列表
    """
    if expr == Zero or expr == 0:
        return []

    if isinstance(expr, Add):
        return list(expr.args)

    return [expr]


def _term_key(term: Any) -> tuple:
    """
    生成项的排序键

    用于稳定排序，减少项顺序差异导致的比较失败。

    Args:
        term: 表达式项

    Returns:
        排序键元组
    """
    # 提取标量系数和算符部分
    if is_local_operator(term):
        scalar, op = extract_scalar_operator(term)
    else:
        scalar = term
        op = None

    # 构造排序键
    key_parts = []

    # 1. 算符类型优先级
    if op is None:
        key_parts.append(0)  # 纯标量
    elif isinstance(op, BasicOperator):
        key_parts.append(1)
    elif isinstance(op, DerivativeOperator):
        key_parts.append(2)
    elif isinstance(op, NormalOrderedOperator):
        key_parts.append(3)
    else:
        key_parts.append(4)

    # 2. 算符的字符串表示
    if op is not None:
        key_parts.append(str(op))
    else:
        key_parts.append("")

    # 3. 标量系数（转换为浮点数用于排序）
    try:
        if isinstance(scalar, (int, float)):
            key_parts.append(float(scalar))
        elif isinstance(scalar, Number):
            key_parts.append(float(scalar))
        else:
            key_parts.append(0.0)
    except (TypeError, ValueError):
        key_parts.append(0.0)

    return tuple(key_parts)


def compare_expressions(
    expr1: Any, expr2: Any, *, canonicalize_first: bool = True
) -> bool:
    """
    比较两个表达式是否相等（结构比较）

    Args:
        expr1: 第一个表达式
        expr2: 第二个表达式
        canonicalize_first: 是否先规范化（默认 True）

    Returns:
        是否相等

    Examples:
        >>> T = BasicOperator("T")
        >>> expr1 = T + 2*T
        >>> expr2 = 3*T
        >>> compare_expressions(expr1, expr2)  # True
    """
    # 规范化
    if canonicalize_first:
        expr1 = canonicalize(expr1)
        expr2 = canonicalize(expr2)

    # 直接比较
    if expr1 == expr2:
        return True

    # 处理零的不同表示
    if (expr1 == Zero or expr1 == 0) and (expr2 == Zero or expr2 == 0):
        return True

    # 处理单位的不同表示
    if (expr1 == One or expr1 == 1) and (expr2 == One or expr2 == 1):
        return True

    # 提取项并排序比较
    terms1 = _extract_terms(expr1)
    terms2 = _extract_terms(expr2)

    if len(terms1) != len(terms2):
        return False

    # 排序后逐项比较
    sorted_terms1 = sorted(terms1, key=_term_key)
    sorted_terms2 = sorted(terms2, key=_term_key)

    for t1, t2 in zip(sorted_terms1, sorted_terms2):
        if t1 != t2:
            # 尝试符号化简后比较
            try:
                diff = sp.simplify(t1 - t2)
                if diff != 0 and diff != Zero:
                    return False
            except Exception:
                return False

    return True


def assert_voa_equal(
    actual: Any,
    expected: Any,
    *,
    msg: Optional[str] = None,
    canonicalize_first: bool = True,
) -> None:
    """
    断言两个 VOA 表达式相等（符号比较）

    Args:
        actual: 实际值
        expected: 期望值
        msg: 自定义错误消息
        canonicalize_first: 是否先规范化（默认 True）

    Raises:
        AssertionError: 如果表达式不相等

    Examples:
        >>> T = BasicOperator("T")
        >>> actual = simplify(NO(T, T))
        >>> expected = NO(T, T)
        >>> assert_voa_equal(actual, expected)
    """
    if canonicalize_first:
        actual = canonicalize(actual)
        expected = canonicalize(expected)

    if not compare_expressions(actual, expected, canonicalize_first=False):
        error_msg = f"\nExpected: {expected}\nActual:   {actual}"
        if msg:
            error_msg = f"{msg}\n{error_msg}"
        raise AssertionError(error_msg)


def _evaluate_numeric(
    expr: Any, subs: Dict[Symbol, Union[int, float]]
) -> Union[float, complex]:
    """
    将表达式数值化

    Args:
        expr: 表达式
        subs: 符号替换字典

    Returns:
        数值结果
    """
    # 处理零
    if expr == Zero or expr == 0:
        return 0.0

    # 处理单位
    if expr == One or expr == 1:
        return 1.0

    # 处理 OPEData
    if isinstance(expr, OPEData):
        # 对每个极点进行数值化
        result = {}
        for pole, coeff in expr.poles.items():
            result[pole] = _evaluate_numeric(coeff, subs)
        return result

    # 处理算符（不能直接数值化）
    if isinstance(expr, Operator):
        raise ValueError(f"Cannot numerically evaluate operator: {expr}")

    # 使用 sympy 的 subs 和 evalf
    try:
        if hasattr(expr, "subs"):
            substituted = expr.subs(subs)
            if hasattr(substituted, "evalf"):
                return complex(substituted.evalf())
            return float(substituted)
        return float(expr)
    except (TypeError, ValueError, AttributeError) as e:
        raise ValueError(f"Cannot numerically evaluate expression: {expr}") from e


def assert_voa_numeric_equal(
    actual: Any,
    expected: Any,
    *,
    subs: Dict[Symbol, Union[int, float]],
    tol: float = 1e-12,
    msg: Optional[str] = None,
) -> None:
    """
    断言两个 VOA 表达式在数值代入后相等（数值比较）

    Args:
        actual: 实际值
        expected: 期望值
        subs: 符号替换字典（例如 {c: 100, beta: 0.1}）
        tol: 数值容差（默认 1e-12）
        msg: 自定义错误消息

    Raises:
        AssertionError: 如果数值不相等

    Examples:
        >>> import sympy as sp
        >>> c = sp.Symbol('c')
        >>> T = BasicOperator("T")
        >>> actual = c * T
        >>> expected = 100 * T
        >>> assert_voa_numeric_equal(actual, expected, subs={c: 100})
    """
    # 先规范化
    actual = canonicalize(actual)
    expected = canonicalize(expected)

    # 处理 OPEData
    if isinstance(actual, OPEData) and isinstance(expected, OPEData):
        # 比较极点数量
        if set(actual.poles.keys()) != set(expected.poles.keys()):
            error_msg = f"\nPole mismatch:\nExpected poles: {set(expected.poles.keys())}\nActual poles:   {set(actual.poles.keys())}"
            if msg:
                error_msg = f"{msg}\n{error_msg}"
            raise AssertionError(error_msg)

        # 逐个极点比较
        for pole in actual.poles.keys():
            actual_coeff = actual.poles[pole]
            expected_coeff = expected.poles[pole]

            try:
                assert_voa_numeric_equal(
                    actual_coeff,
                    expected_coeff,
                    subs=subs,
                    tol=tol,
                    msg=f"Pole {pole} mismatch",
                )
            except AssertionError as e:
                if msg:
                    raise AssertionError(f"{msg}\n{str(e)}") from e
                raise
        return

    # 数值化
    try:
        actual_val = _evaluate_numeric(actual, subs)
        expected_val = _evaluate_numeric(expected, subs)
    except ValueError as e:
        # 如果不能数值化，回退到符号比较
        try:
            assert_voa_equal(actual, expected, msg=msg, canonicalize_first=False)
            return
        except AssertionError:
            raise AssertionError(f"Cannot numerically evaluate: {e}") from e

    # 比较数值
    if isinstance(actual_val, dict) and isinstance(expected_val, dict):
        # 比较字典（用于 OPEData）
        if set(actual_val.keys()) != set(expected_val.keys()):
            error_msg = f"\nKey mismatch:\nExpected: {set(expected_val.keys())}\nActual:   {set(actual_val.keys())}"
            if msg:
                error_msg = f"{msg}\n{error_msg}"
            raise AssertionError(error_msg)

        for key in actual_val.keys():
            if not math.isclose(
                actual_val[key], expected_val[key], rel_tol=tol, abs_tol=tol
            ):
                error_msg = f"\nKey {key} mismatch:\nExpected: {expected_val[key]}\nActual:   {actual_val[key]}"
                if msg:
                    error_msg = f"{msg}\n{error_msg}"
                raise AssertionError(error_msg)
    else:
        # 比较标量
        # 处理复数
        if isinstance(actual_val, complex) or isinstance(expected_val, complex):
            actual_val = complex(actual_val)
            expected_val = complex(expected_val)
            # 分别比较实部和虚部
            real_close = math.isclose(
                actual_val.real, expected_val.real, rel_tol=tol, abs_tol=tol
            )
            imag_close = math.isclose(
                actual_val.imag, expected_val.imag, rel_tol=tol, abs_tol=tol
            )
            if not (real_close and imag_close):
                error_msg = f"\nExpected: {expected_val}\nActual:   {actual_val}\nDifference: {abs(actual_val - expected_val)}"
                if msg:
                    error_msg = f"{msg}\n{error_msg}"
                raise AssertionError(error_msg)
        else:
            # 比较实数
            if not math.isclose(actual_val, expected_val, rel_tol=tol, abs_tol=tol):
                error_msg = f"\nExpected: {expected_val}\nActual:   {actual_val}\nDifference: {abs(actual_val - expected_val)}"
                if msg:
                    error_msg = f"{msg}\n{error_msg}"
                raise AssertionError(error_msg)
