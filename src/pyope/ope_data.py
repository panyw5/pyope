"""
OPE 数据结构模块

本模块定义了用于存储和管理 OPE（算符积展开）数据的类：
- OPEData: 存储 OPE 的极点信息
"""

from typing import Dict, List, Optional, Union, Callable, Any
import sympy as sp
from sympy import sympify, simplify


class OPEData:
    """
    OPE 数据类

    存储算符积展开 A(z)B(w) 的极点信息。
    OPE 的形式为：
        A(z)B(w) = Σ_n {AB}_n(w) / (z-w)^n

    其中 {AB}_n 是第 n 阶极点的系数（局域算符）。

    Attributes:
        poles: 字典，键为极点阶数 n，值为对应的系数
    """

    def __init__(self, poles: Optional[Dict[int, Any]] = None):
        """
        创建 OPEData 实例

        Args:
            poles: 极点字典，键为阶数，值为系数（可以是算符或 sympy 表达式）
        """
        if poles is None:
            self._poles = {}
        else:
            # 过滤掉零系数的极点
            self._poles = {n: coeff for n, coeff in poles.items()
                          if not self._is_zero(coeff)}

    @staticmethod
    def _is_zero(expr: Any) -> bool:
        """
        判断表达式是否为零

        Args:
            expr: 要判断的表达式

        Returns:
            True 如果表达式为零
        """
        if expr == 0:
            return True
        if isinstance(expr, sp.Expr):
            # 使用简单的相等性检查，避免触发复杂的 sympy 简化
            return expr == 0 or expr == sp.Integer(0)
        return False

    @classmethod
    def from_list(cls, pole_list: List[Any]) -> 'OPEData':
        """
        从列表创建 OPEData（Mathematica 风格）

        列表按照从高阶到低阶排列：[pole_n, pole_{n-1}, ..., pole_2, pole_1]
        其中 n 是列表长度。

        Args:
            pole_list: 极点系数列表

        Returns:
            OPEData 实例

        Examples:
            >>> pole_list = [c/2, 0, 2*T, dT]  # [pole_4, pole_3, pole_2, pole_1]
            >>> ope = OPEData.from_list(pole_list)
        """
        n = len(pole_list)
        poles = {}
        for i, coeff in enumerate(pole_list):
            pole_order = n - i
            if not cls._is_zero(coeff):
                poles[pole_order] = coeff
        return cls(poles)

    @property
    def poles(self) -> Dict[int, Any]:
        """返回极点字典的副本"""
        return self._poles.copy()

    @property
    def max_pole(self) -> int:
        """
        返回最高阶极点

        Returns:
            最高阶极点的阶数，如果没有极点则返回 0
        """
        if not self._poles:
            return 0
        return max(self._poles.keys())

    def pole(self, n: int) -> Any:
        """
        获取第 n 阶极点的系数

        Args:
            n: 极点阶数

        Returns:
            第 n 阶极点的系数，如果不存在则返回 0
        """
        return self._poles.get(n, 0)

    def set_pole(self, n: int, coeff: Any):
        """
        设置第 n 阶极点的系数

        Args:
            n: 极点阶数
            coeff: 系数
        """
        if self._is_zero(coeff):
            # 如果系数为零，删除该极点
            self._poles.pop(n, None)
        else:
            self._poles[n] = coeff

    def is_zero(self) -> bool:
        """
        判断 OPE 是否为零（没有任何极点）

        Returns:
            True 如果没有极点
        """
        return len(self._poles) == 0

    def __add__(self, other: 'OPEData') -> 'OPEData':
        """
        OPE 加法

        Args:
            other: 另一个 OPEData

        Returns:
            新的 OPEData，包含两个 OPE 的和
        """
        if not isinstance(other, OPEData):
            return NotImplemented

        # 合并极点
        result_poles = self._poles.copy()
        for n, coeff in other._poles.items():
            if n in result_poles:
                # 同阶极点相加
                new_coeff = result_poles[n] + coeff
                if self._is_zero(new_coeff):
                    del result_poles[n]
                else:
                    result_poles[n] = new_coeff
            else:
                result_poles[n] = coeff

        return OPEData(result_poles)

    def __radd__(self, other):
        """
        右加法（支持 sum() 函数）

        Args:
            other: 左操作数

        Returns:
            加法结果
        """
        if other == 0:
            # 支持 sum([ope1, ope2, ...])
            return self
        return self.__add__(other)

    def __mul__(self, scalar: Any) -> 'OPEData':
        """
        标量乘法（右乘）

        Args:
            scalar: 标量

        Returns:
            新的 OPEData
        """
        result_poles = {n: scalar * coeff for n, coeff in self._poles.items()}
        return OPEData(result_poles)

    def __rmul__(self, scalar: Any) -> 'OPEData':
        """
        标量乘法（左乘）

        Args:
            scalar: 标量

        Returns:
            新的 OPEData
        """
        return self.__mul__(scalar)

    def __neg__(self) -> 'OPEData':
        """
        取负

        Returns:
            新的 OPEData
        """
        return self.__mul__(-1)

    def __sub__(self, other: 'OPEData') -> 'OPEData':
        """
        OPE 减法

        Args:
            other: 另一个 OPEData

        Returns:
            新的 OPEData
        """
        return self.__add__(-other)

    def __eq__(self, other: Any) -> bool:
        """
        相等性比较

        Args:
            other: 另一个对象

        Returns:
            True 如果两个 OPE 相等
        """
        if not isinstance(other, OPEData):
            return False

        # 比较极点集合
        if set(self._poles.keys()) != set(other._poles.keys()):
            return False

        # 比较每个极点的系数
        for n in self._poles.keys():
            if self._poles[n] != other._poles[n]:
                # 尝试使用 sympy 的 equals 方法
                if isinstance(self._poles[n], sp.Expr) and isinstance(other._poles[n], sp.Expr):
                    if not self._poles[n].equals(other._poles[n]):
                        return False
                else:
                    return False

        return True

    def simplify(self, simplify_func: Optional[Callable] = None) -> 'OPEData':
        """
        简化 OPE 中的表达式

        Args:
            simplify_func: 简化函数（默认使用 sympy.simplify）

        Returns:
            简化后的新 OPEData
        """
        if simplify_func is None:
            simplify_func = sp.simplify

        result_poles = {}
        for n, coeff in self._poles.items():
            if isinstance(coeff, sp.Expr):
                result_poles[n] = simplify_func(coeff)
            else:
                result_poles[n] = coeff

        return OPEData(result_poles)

    def __repr__(self) -> str:
        """字符串表示"""
        if self.is_zero():
            return "OPEData({})"

        pole_strs = [f"{n}: {coeff}" for n, coeff in sorted(self._poles.items(), reverse=True)]
        return f"OPEData({{{', '.join(pole_strs)}}})"

    def __str__(self) -> str:
        """用户友好的字符串表示"""
        if self.is_zero():
            return "0"

        terms = []
        for n in sorted(self._poles.keys(), reverse=True):
            coeff = self._poles[n]
            if n == 1:
                terms.append(f"({coeff})/(z-w)")
            else:
                terms.append(f"({coeff})/(z-w)^{n}")

        return " + ".join(terms)
