"""
OPE 注册表模块

本模块实现了 OPE 定义的存储和管理机制：
- OPERegistry: 注册表类，存储算符的 parity 和 OPE 定义
- ope_registry: 全局注册表实例
- Bosonic, Fermionic: 声明算符类型的辅助函数
"""

from typing import Dict, Tuple, Optional, Any, Union
import sympy as sp

from .operators import Operator, BasisOperator
from .ope_data import OPEData


class OPERegistry:
    """
    OPE 注册表类

    存储和管理算符的 OPE 定义。功能包括：
    - 注册算符的 parity（玻色/费米）
    - 存储基本算符的 OPE 定义
    - 查询已定义的 OPE
    - 管理算符的排序位置

    Attributes:
        _opes: OPE 定义字典，键为 (A, B) 元组
        _parities: 算符 parity 字典
        _positions: 算符位置字典（用于排序）
        _position_counter: 位置计数器
    """

    def __init__(self):
        """初始化空的注册表"""
        self._opes: Dict[Tuple[Any, Any], OPEData] = {}
        self._parities: Dict[Any, int] = {}
        self._positions: Dict[Any, int] = {}
        self._position_counter: int = 0

    def register_operator(self, operator: Any, parity: int) -> None:
        """
        注册算符及其 parity

        Args:
            operator: 要注册的算符（可以是 Operator 实例或模式）
            parity: 算符的 parity（0 表示玻色子，1 表示费米子）

        Raises:
            ValueError: 如果 parity 不是 0 或 1
        """
        if parity not in (0, 1):
            raise ValueError(f"Parity must be 0 (bosonic) or 1 (fermionic), got {parity}")

        # 如果算符已经注册，发出警告
        if self.is_registered(operator):
            import warnings
            warnings.warn(f"Operator {operator} is already registered, overwriting parity")

        # 注册 parity
        self._parities[operator] = parity

        # 分配位置（用于排序）
        if operator not in self._positions:
            self._positions[operator] = self._position_counter
            self._position_counter += 1

    def is_registered(self, operator: Any) -> bool:
        """
        检查算符是否已注册

        Args:
            operator: 要检查的算符

        Returns:
            True 如果算符已注册
        """
        return operator in self._parities

    def get_parity(self, operator: Any) -> Optional[int]:
        """
        获取算符的 parity

        Args:
            operator: 算符

        Returns:
            parity 值（0 或 1），如果未注册则返回 None
        """
        return self._parities.get(operator)

    def get_position(self, operator: Any) -> Optional[int]:
        """
        获取算符的位置（用于排序）

        Args:
            operator: 算符

        Returns:
            位置值，如果未注册则返回 None
        """
        return self._positions.get(operator)

    def define_ope(self, left: Any, right: Any, ope_data: OPEData) -> None:
        """
        定义两个算符的 OPE

        Args:
            left: 左侧算符
            right: 右侧算符
            ope_data: OPE 数据

        Examples:
            >>> T = BasisOperator("T", bosonic=True)
            >>> registry.define_ope(T, T, OPEData({2: 2*T, 1: d(T)}))
        """
        # 创建规范化的键（使用算符的字符串表示）
        key = self._make_key(left, right)
        self._opes[key] = ope_data

    def get_ope(self, left: Any, right: Any) -> Optional[OPEData]:
        """
        查询两个算符的 OPE

        Args:
            left: 左侧算符
            right: 右侧算符

        Returns:
            OPEData 如果已定义，否则返回 None
        """
        key = self._make_key(left, right)
        return self._opes.get(key)

    def has_ope(self, left: Any, right: Any) -> bool:
        """
        检查是否定义了两个算符的 OPE

        Args:
            left: 左侧算符
            right: 右侧算符

        Returns:
            True 如果已定义
        """
        key = self._make_key(left, right)
        return key in self._opes

    def _make_key(self, left: Any, right: Any) -> Tuple[str, str]:
        """
        创建 OPE 字典的键

        使用算符的字符串表示作为键，以便支持模式匹配

        Args:
            left: 左侧算符
            right: 右侧算符

        Returns:
            (left_str, right_str) 元组
        """
        # 对于 Operator 实例，使用 name 属性
        if isinstance(left, Operator):
            left_key = left.name
        else:
            left_key = str(left)

        if isinstance(right, Operator):
            right_key = right.name
        else:
            right_key = str(right)

        return (left_key, right_key)

    def clear(self) -> None:
        """清空注册表（主要用于测试）"""
        self._opes.clear()
        self._parities.clear()
        self._positions.clear()
        self._position_counter = 0

    def __repr__(self) -> str:
        """字符串表示"""
        return f"OPERegistry(operators={len(self._parities)}, opes={len(self._opes)})"


# 全局注册表实例
ope_registry = OPERegistry()
"""
全局 OPE 注册表

这是一个全局的 OPERegistry 实例，用于存储所有的算符和 OPE 定义。
在大多数情况下，应该使用这个全局实例而不是创建新的注册表。
"""


# 辅助函数

def Bosonic(*operators) -> None:
    """
    声明算符为玻色算符

    将一个或多个算符注册为玻色算符（parity = 0）。

    Args:
        *operators: 要声明的算符

    Examples:
        >>> T = BasisOperator("T")
        >>> J = BasisOperator("J")
        >>> Bosonic(T, J)
    """
    for op in operators:
        ope_registry.register_operator(op, parity=0)


def Fermionic(*operators) -> None:
    """
    声明算符为费米算符

    将一个或多个算符注册为费米算符（parity = 1）。

    Args:
        *operators: 要声明的算符

    Examples:
        >>> psi = BasisOperator("ψ")
        >>> chi = BasisOperator("χ")
        >>> Fermionic(psi, chi)
    """
    for op in operators:
        ope_registry.register_operator(op, parity=1)


class OPEDefiner:
    """
    OPE 定义器类

    提供类似 Mathematica 的 OPE[A, B] = ... 语法。

    使用方式：
        >>> OPE = OPEDefiner()
        >>> OPE[T, T] = OPEData({2: 2*T, 1: d(T)})
        >>> ope_data = OPE(T, T)  # 或 OPE[T, T]
    """

    def __init__(self, registry: Optional[OPERegistry] = None):
        """
        创建 OPE 定义器

        Args:
            registry: 使用的注册表，默认使用全局注册表
        """
        self.registry = registry if registry is not None else ope_registry

    def __setitem__(self, key: Tuple[Any, Any], value: OPEData) -> None:
        """
        定义 OPE: OPE[A, B] = OPEData(...)

        Args:
            key: (left, right) 算符元组
            value: OPE 数据
        """
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("Key must be a tuple of two operators: OPE[A, B]")

        left, right = key

        if not isinstance(value, OPEData):
            raise TypeError(f"Value must be OPEData, got {type(value)}")

        self.registry.define_ope(left, right, value)

    def __getitem__(self, key: Tuple[Any, Any]) -> Optional[OPEData]:
        """
        查询 OPE: OPE[A, B]

        Args:
            key: (left, right) 算符元组

        Returns:
            OPEData 如果已定义，否则返回 None
        """
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("Key must be a tuple of two operators: OPE[A, B]")

        left, right = key
        return self.registry.get_ope(left, right)

    def __call__(self, left: Any, right: Any):
        """
        计算 OPE: OPE(A, B)

        这个方法会被 api.py 中的 OPEComputer 类覆盖。

        Args:
            left: 左侧算符
            right: 右侧算符

        Returns:
            OPEData 实例
        """
        raise NotImplementedError("OPE computation is implemented in api.py")

    @staticmethod
    def make(data):
        """
        创建 OPEData（静态方法）

        这个方法会被 api.py 中的 MakeOPE 函数覆盖。

        Args:
            data: 极点列表或 OPEData 实例

        Returns:
            OPEData 实例
        """
        raise NotImplementedError("MakeOPE is implemented in api.py")
