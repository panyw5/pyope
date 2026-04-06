"""
常数算符模块

本模块定义了 VOA 计算中使用的常数算符：
- ConstantOperator: 常数算符基类
- One: 单位算符
- Zero: 零算符
- Delta: Kronecker delta 函数
"""

from typing import Any, Union
import sympy as sp
from sympy import Function, KroneckerDelta

from .operators import Operator


class ConstantOperator(Operator):
    """
    常数算符类

    表示特殊的常数算符，如单位算符 One、零算符 Zero 等。
    常数算符总是玻色子（parity = 0）。

    Attributes:
        name: 算符名称
    """

    _name: str  # 类型注解，告诉类型检查器这个属性存在

    def __new__(cls, name: str, **assumptions):
        """
        创建常数算符

        Args:
            name: 算符名称
        """
        obj = Operator.__new__(cls, name, **assumptions)
        obj._name = name  # type: ignore[attr-defined]
        return obj

    @property
    def is_bosonic(self) -> bool:
        """常数算符总是玻色子"""
        return True

    @property
    def is_fermionic(self) -> bool:
        """常数算符不是费米子"""
        return False

    @property
    def parity(self) -> int:
        """常数算符的 parity 总是 0（玻色子）"""
        return 0

    def __eq__(self, other):
        if isinstance(other, ConstantOperator):
            return self._name == other._name  # type: ignore[attr-defined]
        if self._name == "Zero":  # type: ignore[attr-defined]
            return other == 0 or other == sp.Integer(0)
        return False

    def __bool__(self):
        if self._name == "Zero":  # type: ignore[attr-defined]
            return False
        return True

    def __mul__(self, other) -> sp.Expr:
        """乘法运算"""
        if self._name == "Zero":  # type: ignore[attr-defined]
            return self  # Zero * anything = Zero
        if self._name == "One":  # type: ignore[attr-defined]
            return other  # One * anything = anything
        return sp.Mul(self, other)

    def __add__(self, other) -> sp.Expr:
        """加法运算"""
        if self._name == "Zero":  # type: ignore[attr-defined]
            other_expr = sp.sympify(other)
            if other_expr == 0 or other_expr == sp.Integer(0):
                return self
            return other_expr
        return sp.Add(self, other)

    def __radd__(self, other) -> sp.Expr:
        """右加法运算"""
        if self._name == "Zero":  # type: ignore[attr-defined]
            other_expr = sp.sympify(other)
            if other_expr == 0 or other_expr == sp.Integer(0):
                return self
            return other_expr
        return sp.Add(other, self)

    def __sub__(self, other) -> sp.Expr:
        """减法运算"""
        if self._name == "Zero":  # type: ignore[attr-defined]
            other_expr = sp.sympify(other)
            if other_expr == 0 or other_expr == sp.Integer(0):
                return self
            return -other_expr
        return sp.Add(self, -sp.sympify(other))

    def __rsub__(self, other) -> sp.Expr:
        """右减法运算"""
        if self._name == "Zero":  # type: ignore[attr-defined]
            other_expr = sp.sympify(other)
            if other_expr == 0 or other_expr == sp.Integer(0):
                return self
            return other_expr
        return sp.Add(sp.sympify(other), -self)

    def __rmul__(self, other) -> sp.Expr:
        """右乘法运算"""
        if self._name == "Zero":  # type: ignore[attr-defined]
            return self  # anything * Zero = Zero

        if self._name == "One":  # type: ignore[attr-defined]
            other_expr = sp.sympify(other)

            # 1 * One -> One
            if other_expr == 1 or other_expr == sp.Integer(1):
                return self

            # (operator expression) * One -> operator expression
            # Keep symbolic scalars (like c) as c*One.
            if other_expr.has(Operator):
                return other_expr

            return sp.Mul(other_expr, self)

        return sp.Mul(other, self)

    def __hash__(self):
        """哈希值"""
        return hash(self._name)  # type: ignore[attr-defined]

    def __repr__(self):
        """字符串表示"""
        return f"ConstantOperator('{self._name}')"  # type: ignore[attr-defined]

    def _latex(self, printer=None):
        """
        LaTeX 渲染

        返回 LaTeX 格式的字符串，用于 sympy 的 latex() 函数。

        Returns:
            LaTeX 格式的字符串
        """
        # One 渲染为 1，Zero 渲染为 0
        if self._name == "One":  # type: ignore[attr-defined]
            return "1"
        elif self._name == "Zero":  # type: ignore[attr-defined]
            return "0"
        else:
            return self._name  # type: ignore[attr-defined]


# 预定义的常数算符实例

# 单位算符 One
One = ConstantOperator("One")
"""
单位算符

在 VOA 中，One 表示单位元素，满足：
- NO(One, A) = A
- NO(A, One) = A
- OPE(One, A) 只有 0 阶极点，系数为 A
"""

# 零算符 Zero
Zero = ConstantOperator("Zero")
"""
零算符

在 VOA 中，Zero 表示零元素，满足：
- A + Zero = A
- c * Zero = Zero (对任意标量 c)
- NO(Zero, A) = Zero
- NO(A, Zero) = Zero
"""


# Delta 函数


class DeltaFunction(Function):
    """
    Kronecker Delta 函数

    Delta(i, j) 返回：
    - 1 如果 i == j
    - 0 如果 i != j
    - 符号表达式如果无法确定

    这是对 sympy.KroneckerDelta 的包装，提供更简洁的接口。
    """

    @classmethod
    def eval(cls, i, j):  # type: ignore[override]
        """
        计算 Delta(i, j)

        Args:
            i: 第一个索引
            j: 第二个索引

        Returns:
            1, 0, 或符号表达式
        """
        # 如果两个参数相同
        if i == j:
            return sp.Integer(1)

        # 如果都是数值且不相等
        if i.is_number and j.is_number:
            return sp.Integer(0)

        # 否则返回 KroneckerDelta
        return KroneckerDelta(i, j)


def Delta(i: Any, j: Any):
    """
    Kronecker Delta 函数

    计算 δ_ij，即：
    - 返回 1 如果 i == j
    - 返回 0 如果 i != j（且都是数值）
    - 返回符号表达式如果无法确定

    Args:
        i: 第一个索引
        j: 第二个索引

    Returns:
        1, 0, 或 sympy 表达式

    Examples:
        >>> Delta(1, 1)
        1
        >>> Delta(1, 2)
        0
        >>> i = sp.Symbol('i')
        >>> Delta(i, i)
        1
        >>> j = sp.Symbol('j')
        >>> Delta(i, j)  # 返回符号表达式
        KroneckerDelta(i, j)
    """
    # 转换为 sympy 对象
    i_sym = sp.sympify(i)
    j_sym = sp.sympify(j)

    # 使用 DeltaFunction 进行计算
    result = DeltaFunction(i_sym, j_sym)

    # 如果结果是 DeltaFunction 实例，返回 KroneckerDelta
    if isinstance(result, DeltaFunction):
        return KroneckerDelta(i_sym, j_sym)

    return result
