"""
算符类定义模块

本模块定义了 VOA 计算中使用的各种算符类型：
- Operator: 算符基类
- BasisOperator: 基本算符（用户定义的基础局域算符）
- DerivativeOperator: 导数算符
- NormalOrderedOperator: 正规序算符
"""

from typing import Optional, Tuple, Any
import sympy as sp
from sympy.core.symbol import Symbol


class Operator(Symbol):
    """
    算符基类

    继承自 sympy.Symbol 以便与 sympy 的符号计算系统集成。
    所有算符都具有 parity（奇偶性）属性。
    """

    def __new__(cls, name: str, **assumptions):
        """创建新的算符实例"""
        obj = Symbol.__new__(cls, name, **assumptions)
        return obj

    @property
    def parity(self) -> int:
        """
        算符的 parity（奇偶性）

        Returns:
            0 表示玻色子（bosonic）
            1 表示费米子（fermionic）
        """
        raise NotImplementedError("Subclasses must implement parity property")


class BasisOperator(Operator):
    """
    基本算符类

    表示用户定义的基础局域算符，具有特定的 conformal weight 和 parity。

    Attributes:
        name: 算符名称
        bosonic: 是否为玻色算符（True）或费米算符（False）
        indexed: 是否支持索引（如 J[i]）
        conformal_weight: 共形权重（可选）
        indices: 索引元组（仅当 indexed=True 时使用）
    """

    def __new__(cls, name: str, bosonic: bool = True, indexed: bool = False,
                conformal_weight: Optional[float] = None, indices: Optional[Tuple] = None,
                base_name: Optional[str] = None, **assumptions):
        """
        创建基本算符

        Args:
            name: 算符名称
            bosonic: True 表示玻色子，False 表示费米子
            indexed: 是否支持索引
            conformal_weight: 共形权重
            indices: 索引元组（用于索引算符）
            base_name: 基础名称（用于索引算符保持原始名称）
        """
        obj = Operator.__new__(cls, name, **assumptions)
        obj._bosonic = bosonic
        obj._indexed = indexed
        obj._conformal_weight = conformal_weight
        obj._indices = indices if indices is not None else ()
        obj._base_name = base_name if base_name is not None else name
        return obj

    @property
    def is_bosonic(self) -> bool:
        """是否为玻色算符"""
        return self._bosonic

    @property
    def is_fermionic(self) -> bool:
        """是否为费米算符"""
        return not self._bosonic

    @property
    def parity(self) -> int:
        """
        算符的 parity

        Returns:
            0 表示玻色子
            1 表示费米子
        """
        return 0 if self._bosonic else 1

    @property
    def base_name(self) -> str:
        """基础名称（不包含索引）"""
        return self._base_name

    @property
    def indices(self) -> Tuple:
        """索引元组"""
        return self._indices

    @property
    def conformal_weight(self) -> Optional[float]:
        """共形权重"""
        return self._conformal_weight

    def __getitem__(self, index):
        """
        支持索引操作，如 J[i]

        Args:
            index: 索引（可以是 sympy Symbol 或其他类型）

        Returns:
            带索引的新 BasisOperator 实例
        """
        if not self._indexed:
            raise TypeError(f"Operator {self.name} is not indexed")

        # 确保索引是元组
        if not isinstance(index, tuple):
            index = (index,)

        # 创建新的索引算符
        indexed_name = f"{self._base_name}_{index}"
        return BasisOperator(
            indexed_name,
            bosonic=self._bosonic,
            indexed=True,
            conformal_weight=self._conformal_weight,
            indices=index,
            base_name=self._base_name  # 保持原始的 base_name
        )

    def __eq__(self, other):
        """相等性比较"""
        if not isinstance(other, BasisOperator):
            return False
        return (self._base_name == other._base_name and
                self._bosonic == other._bosonic and
                self._indices == other._indices)

    def __hash__(self):
        """哈希值"""
        return hash((self._base_name, self._bosonic, self._indices))

    def __repr__(self):
        """字符串表示"""
        if self._indices:
            indices_str = ",".join(str(i) for i in self._indices)
            return f"BasisOperator('{self._base_name}'[{indices_str}], bosonic={self._bosonic})"
        return f"BasisOperator('{self.name}', bosonic={self._bosonic})"


class DerivativeOperator(Operator):
    """
    导数算符类

    表示对算符的导数，如 ∂A(z) 或 ∂^n A(z)。

    Attributes:
        base: 被求导的基础算符
        order: 导数阶数
    """

    def __new__(cls, base: Operator, order: int = 1, **assumptions):
        """
        创建导数算符

        Args:
            base: 被求导的算符
            order: 导数阶数（默认为 1）
        """
        # 生成导数算符的名称
        if order == 1:
            name = f"∂{base.name}"
        else:
            name = f"∂^{order}{base.name}"

        obj = Operator.__new__(cls, name, **assumptions)
        obj._base = base
        obj._order = order
        return obj

    @property
    def base(self) -> Operator:
        """被求导的基础算符"""
        return self._base

    @property
    def order(self) -> int:
        """导数阶数"""
        return self._order

    @property
    def parity(self) -> int:
        """
        导数算符的 parity 与基础算符相同
        """
        return self._base.parity

    def __eq__(self, other):
        """相等性比较"""
        if not isinstance(other, DerivativeOperator):
            return False
        return self._base == other._base and self._order == other._order

    def __hash__(self):
        """哈希值"""
        return hash((self._base, self._order))

    def __repr__(self):
        """字符串表示"""
        if self._order == 1:
            return f"d({self._base})"
        return f"d^{self._order}({self._base})"


class NormalOrderedOperator(Operator):
    """
    正规序算符类

    表示两个算符的正规序乘积 NO(AB)。

    Attributes:
        left: 左侧算符
        right: 右侧算符
        factors: 算符元组 (left, right)
    """

    def __new__(cls, left: Operator, right: Operator, **assumptions):
        """
        创建正规序算符

        Args:
            left: 左侧算符
            right: 右侧算符
        """
        # 生成正规序算符的名称
        name = f"NO({left.name},{right.name})"

        obj = Operator.__new__(cls, name, **assumptions)
        obj._left = left
        obj._right = right
        obj._factors = (left, right)
        return obj

    @property
    def left(self) -> Operator:
        """左侧算符"""
        return self._left

    @property
    def right(self) -> Operator:
        """右侧算符"""
        return self._right

    @property
    def factors(self) -> Tuple[Operator, Operator]:
        """算符元组"""
        return self._factors

    @property
    def parity(self) -> int:
        """
        正规序算符的 parity 是两个算符 parity 之和模 2
        """
        return (self._left.parity + self._right.parity) % 2

    def __eq__(self, other):
        """相等性比较"""
        if not isinstance(other, NormalOrderedOperator):
            return False
        return self._left == other._left and self._right == other._right

    def __hash__(self):
        """哈希值"""
        return hash((self._left, self._right))

    def __repr__(self):
        """字符串表示"""
        return f"NO({self._left}, {self._right})"


# 辅助函数

def d(operator: Operator, order: int = 1) -> DerivativeOperator:
    """
    计算算符的导数

    Args:
        operator: 要求导的算符
        order: 导数阶数（默认为 1）

    Returns:
        DerivativeOperator 实例

    Examples:
        >>> T = BasisOperator("T", bosonic=True)
        >>> dT = d(T)  # 一阶导数
        >>> d2T = d(T, 2)  # 二阶导数
        >>> d(d(T))  # 等价于 d(T, 2)
    """
    # 如果 operator 已经是 DerivativeOperator，累加阶数
    if isinstance(operator, DerivativeOperator):
        return DerivativeOperator(operator.base, operator.order + order)
    return DerivativeOperator(operator, order)


def dn(order: int, operator: Operator) -> DerivativeOperator:
    """
    计算算符的 n 阶导数（参数顺序与 d 不同）

    Args:
        order: 导数阶数
        operator: 要求导的算符

    Returns:
        DerivativeOperator 实例

    Examples:
        >>> T = BasisOperator("T", bosonic=True)
        >>> d3T = dn(3, T)  # 三阶导数
        >>> dn(2, d(T))  # 等价于 dn(3, T)
    """
    # 如果 operator 已经是 DerivativeOperator，累加阶数
    if isinstance(operator, DerivativeOperator):
        return DerivativeOperator(operator.base, operator.order + order)
    return DerivativeOperator(operator, order)
