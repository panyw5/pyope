"""
局域算符模块

本模块定义了 VOA 中的局域算符类型系统：
- LocalOperator: 局域算符基类/协议
- OperatorSum: 算符和
- OperatorProduct: 算符积
- is_local_operator: 判断是否为局域算符
- extract_scalar_operator: 从表达式中提取标量和算符
"""

from typing import Union, Tuple, Any
import sympy as sp
from sympy import Add, Mul, Number, Symbol

from .operators import (
    Operator,
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
)
from .exceptions import IllegalOperatorProductError


# 类型别名
LocalOperatorType = Union[Operator, sp.Expr]


class LocalOperator:
    """
    局域算符基类/协议

    这是一个标记类，用于表示 VOA 中的局域算符。
    实际的局域算符可以是：
    - BasisOperator 实例
    - DerivativeOperator 实例
    - NormalOrderedOperator 实例
    - sympy 表达式（包含上述算符的线性组合）
    """

    pass


class OperatorSum(Add):
    """
    算符和类

    表示多个局域算符的线性组合。
    继承自 sympy.Add 以便与 sympy 的符号计算系统集成。

    注意：在实际使用中，sympy 会自动将算符的加法转换为 Add 表达式，
    所以这个类主要用于类型标记和文档目的。
    """

    pass


class OperatorProduct(Mul):
    """
    算符积类

    表示算符与标量的乘积。
    继承自 sympy.Mul 以便与 sympy 的符号计算系统集成。

    注意：在实际使用中，sympy 会自动将算符的乘法转换为 Mul 表达式，
    所以这个类主要用于类型标记和文档目的。
    """

    pass


def is_local_operator(expr: Any) -> bool:
    """
    判断一个表达式是否为局域算符

    Args:
        expr: 要检查的表达式

    Returns:
        True 如果表达式是局域算符，否则 False

    Examples:
        >>> T = BasisOperator("T")
        >>> is_local_operator(T)  # True
        >>> is_local_operator(T + T)  # True
        >>> is_local_operator(2 * T)  # True
        >>> is_local_operator(sp.Symbol('x'))  # False
    """
    # 直接是 Operator 实例
    if isinstance(expr, Operator):
        return True

    # 是 sympy 表达式
    if isinstance(expr, sp.Expr):
        # 检查是否所有的原子符号都是 Operator
        atoms = expr.atoms(Symbol)
        if not atoms:
            # 没有符号，可能是纯数字
            return False
        # 所有符号都必须是 Operator
        return all(isinstance(atom, Operator) for atom in atoms)

    return False


def assert_no_illegal_operator_mul(expr: Any, context: str) -> None:
    """
    检测表达式中是否存在非法的算符乘积

    非法的算符乘积是指 Mul 中包含两个或更多含 Operator 的因子。
    例如：Mul(T, J)、Mul(T, Add(J, W))、Mul(NO(T,J), W) 等都是非法的。

    合法的情况：
    - Mul(2, T)：标量乘法
    - Mul(c, T)：符号标量乘法（c 不含 Operator）
    - Mul(2, Add(T, J))：标量乘以算符和

    Args:
        expr: 要检查的表达式
        context: 错误发生的上下文（用于错误消息）

    Raises:
        IllegalOperatorProductError: 如果检测到非法的算符乘积
    """
    if not isinstance(expr, sp.Expr):
        return

    # 单个 Operator 节点不可能直接包含非法的 Mul(operator, operator)。
    # 这条快速路径可以避免在 NO/simplify 热路径里反复做整树扫描。
    if isinstance(expr, Operator):
        return

    for m in expr.atoms(Mul):
        op_like_factors = [
            a for a in m.args if isinstance(a, sp.Expr) and a.has(Operator)
        ]
        if len(op_like_factors) >= 2:
            raise IllegalOperatorProductError(
                m,
                context=context,
                hint="Use NO(A, B) to represent operator products (normal ordering).",
            )


def extract_scalar_operator(expr: sp.Expr) -> Tuple[sp.Expr, Union[Operator, sp.Expr]]:
    """
    从表达式中提取标量系数和算符部分

    将形如 c * A 的表达式分解为 (c, A)，其中 c 是标量，A 是算符。

    Args:
        expr: 要分解的表达式

    Returns:
        (scalar, operator) 元组

    Examples:
        >>> T = BasisOperator("T")
        >>> extract_scalar_operator(2 * T)
        (2, T)
        >>> extract_scalar_operator(T)
        (1, T)
    """
    # 如果直接是算符，系数为 1
    if isinstance(expr, Operator):
        return (sp.Integer(1), expr)

    # 如果是乘法表达式
    if isinstance(expr, Mul):
        from .constants import One, ConstantOperator

        args = expr.args

        # 热路径快速分支：大多数情况是二元乘法 coeff * Operator。
        if len(args) == 2:
            left, right = args
            if isinstance(left, Operator) and not isinstance(right, Operator):
                if right is One or (
                    isinstance(right, ConstantOperator) and right == One
                ):
                    return (sp.Integer(1), left)
                return (right, left)
            if isinstance(right, Operator) and not isinstance(left, Operator):
                if left is One or (isinstance(left, ConstantOperator) and left == One):
                    return (sp.Integer(1), right)
                return (left, right)

        has_one = False
        scalar_parts = []
        operator_parts = []

        for arg in args:
            if arg is One or (isinstance(arg, ConstantOperator) and arg == One):
                has_one = True
                continue

            if isinstance(arg, sp.Expr) and arg.has(Operator):
                operator_parts.append(arg)
            else:
                scalar_parts.append(arg)

        if len(scalar_parts) == 0:
            coeff = sp.Integer(1)
        elif len(scalar_parts) == 1:
            coeff = scalar_parts[0]
        else:
            coeff = Mul(*scalar_parts)

        if len(operator_parts) == 0 and has_one:
            operator_parts.append(One)

        if len(operator_parts) == 0:
            return (coeff, sp.Integer(1))
        elif len(operator_parts) == 1:
            return (coeff, operator_parts[0])
        else:
            raise IllegalOperatorProductError(
                expr,
                context="extract_scalar_operator",
                hint="Operator products must be represented explicitly via NO(A, B).",
            )

    # 如果是加法表达式，不能简单分解
    if isinstance(expr, Add):
        # 返回系数 1 和原表达式
        return (sp.Integer(1), expr)

    # 其他情况
    return (sp.Integer(1), expr)


def get_operator_parity(expr: LocalOperatorType) -> int:
    """
    获取局域算符的 parity

    Args:
        expr: 局域算符表达式

    Returns:
        parity 值（0 或 1）

    Raises:
        ValueError: 如果表达式不是有效的局域算符或 parity 不一致
    """
    # 直接是 Operator 实例
    if isinstance(expr, Operator):
        return expr.parity

    # 是 sympy 表达式
    if isinstance(expr, sp.Expr):
        # 对于加法，所有项必须有相同的 parity
        if isinstance(expr, Add):
            parities = set()
            for arg in expr.args:
                parities.add(get_operator_parity(arg))
            if len(parities) > 1:
                raise ValueError("Operators in sum have inconsistent parities")
            return parities.pop() if parities else 0

        # 对于乘法，提取算符部分
        _, operator = extract_scalar_operator(expr)
        if isinstance(operator, Operator):
            return operator.parity
        elif isinstance(operator, sp.Expr):
            return get_operator_parity(operator)

    raise ValueError(f"Cannot determine parity of {expr}")


# 为了方便，添加一些辅助函数


def simplify_with_sympy(expr: LocalOperatorType) -> LocalOperatorType:
    """
    使用 SymPy 简化算符表达式

    Args:
        expr: 算符表达式

    Returns:
        简化后的表达式
    """
    if isinstance(expr, sp.Expr):
        return sp.simplify(expr)
    return expr


def collect_operator_terms(expr: LocalOperatorType) -> dict:
    """
    收集算符表达式中的同类项

    Args:
        expr: 算符表达式

    Returns:
        字典，键为算符，值为系数
    """
    terms: dict[Any, sp.Expr] = {}

    def add_term(operator: Any, coeff: sp.Expr) -> None:
        coeff = sp.sympify(coeff)
        if coeff == 0:
            return
        terms[operator] = sp.sympify(terms.get(operator, 0) + coeff)
        if terms[operator] == 0:
            del terms[operator]

    def visit(node: LocalOperatorType, scalar: sp.Expr = sp.Integer(1)) -> None:
        scalar = sp.sympify(scalar)
        if scalar == 0:
            return

        if isinstance(node, Operator):
            add_term(node, scalar)
            return

        if isinstance(node, Add):
            for arg in node.args:
                visit(arg, scalar)
            return

        if isinstance(node, Mul):
            coeff, op = extract_scalar_operator(node)
            visit(op, scalar * sp.sympify(coeff))
            return

        if isinstance(node, sp.Expr) and not node.has(Operator):
            add_term(sp.Integer(1), scalar * node)
            return

        add_term(node, scalar)

    visit(expr)
    return terms
